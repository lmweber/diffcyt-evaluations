##########################################################################################
# Script to run methods
# 
# - method: cydar
# - data set: BCR-XL-sim
# 
# - using all markers (does not work correctly: does not detect any DA hyperspheres)
# 
# Lukas Weber, February 2018
##########################################################################################


library(cydar)
library(flowCore)
library(ncdfFlow)
library(edgeR)
library(SummarizedExperiment)
library(Rtsne)


DIR_BENCHMARK <- "../../../../../benchmark_data/BCR_XL_sim/data/main"
DIR_CYDAR_FILES <- "../../../../cydar_files/BCR_XL_sim/supp_all_markers"
DIR_RDATA <- "../../../../RData/BCR_XL_sim/comparisons_cydar"
DIR_SESSION_INFO <- "../../../../session_info/BCR_XL_sim/comparisons_cydar"




###########################
# Load data, pre-processing
###########################

# filenames

files <- list.files(DIR_BENCHMARK, pattern = "\\.fcs$", full.names = TRUE)
files_base <- files[grep("base\\.fcs$", files)]
files_spike <- files[grep("spike\\.fcs$", files)]

files_load <- c(files_base, files_spike)
files_load

# load data as 'ncdfFlowSet' object (for cydar)

d_input <- read.ncdfFlowSet(files_load, transformation = FALSE, truncate_max_range = FALSE)

# sample information

sample_IDs <- gsub("^BCR_XL_sim_", "", 
                   gsub("\\.fcs$", "", basename(files_load)))
sample_IDs

group_IDs <- factor(gsub("^.*_", "", sample_IDs), levels = c("base", "spike"))
group_IDs

patient_IDs <- factor(gsub("_.*$", "", sample_IDs))
patient_IDs

sample_info <- data.frame(group_IDs, patient_IDs, sample_IDs)
sample_info

# marker information

# indices of all marker columns, lineage markers, and functional markers
# (10 surface markers / 14 functional markers; see Bruggner et al. 2014, Table 1)
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)

marker_names <- colnames(d_input[[1]])
marker_names <- gsub("\\(.*$", "", marker_names)

is_marker <- is_celltype_marker <- is_state_marker <- rep(FALSE, length(marker_names))

is_marker[cols_markers] <- TRUE
is_celltype_marker[cols_lineage] <- TRUE
is_state_marker[cols_func] <- TRUE

marker_info <- data.frame(marker_names, is_marker, is_celltype_marker, is_state_marker)
marker_info




######################################
# Additional pre-processing for Citrus
######################################

# --------------
# markers to use
# --------------

cols_to_use <- cols_markers




################
# cydar pipeline
################

# following steps in Bioconductor vignette


# -----------------------------
# Pre-processing of intensities
# -----------------------------

runtime_pre <- system.time({
  
  # Pooling cells together
  # note: keep all cells
  
  pool.ff <- poolCells(d_input, equalize = FALSE)
  
  
  # Estimating transformation parameters
  # note: we use 'asinh' transform with 'cofactor' = 5 (see Bendall et al. 2011, Supp. Fig. S2), 
  # instead of 'logicle' as shown in Bioconductor vignette
  
  cofactor <- 5
  
  transf_exprs <- asinh(exprs(pool.ff)[, cols_markers] / cofactor)
  
  proc.ff <- pool.ff
  exprs(proc.ff)[, cols_markers] <- transf_exprs
  
  
  # Gating out uninteresting cells
  # note: not required for this data set
  
  # keep only markers of interest
  processed.exprs <- proc.ff[, cols_to_use]
  
  
  # Normalizing intensities across batches
  # note: not required for this data set
  
})


# --------------------------------
# Counting cells into hyperspheres
# --------------------------------

runtime_count <- system.time({
  
  # note: data object must be list (or ncdfFlowSet object), with one list item per sample
  
  # convert processed data to named list (this step is not shown in Bioconductor vignette)
  
  # check sample order
  sampleNames(d_input)
  sample_IDs
  
  n_cells <- as.vector(fsApply(d_input, nrow))
  sample_IDs_rep <- rep(sample_IDs, n_cells)
  sample_IDs_rep <- factor(sample_IDs_rep, levels = unique(sample_IDs_rep))
  
  d_list <- split.data.frame(exprs(processed.exprs), sample_IDs_rep)
  
  
  set.seed(1234)
  
  # prepare 'CyData' object
  cd <- prepareCellData(d_list)
  
  # assign cells to hyperspheres
  cd <- countCells(cd, tol = 0.5)
  
  # check hypersphere counts
  head(assay(cd))
  dim(assay(cd))
  
  # check hypersphere positions (median intensities)
  head(intensities(cd))
  
})


# --------------------------------------------------------------------------
# Test for differentially abundant (DA) hyperspheres and control spatial FDR
# --------------------------------------------------------------------------

# using 'edgeR' for testing

runtime_test <- system.time({
  
  # create object for edgeR
  y <- DGEList(assay(cd), lib.size = cd$totals)
  
  # filtering
  keep <- aveLogCPM(y) >= aveLogCPM(1, mean(cd$totals))
  cd <- cd[keep, ]
  y <- y[keep, ]
  
  # design matrix
  # note: including 'patient_IDs' for paired design
  # note: design matrix specification with intercept term
  design <- model.matrix(~ group_IDs + patient_IDs)
  
  # estimate dispersions and fit models
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  
  # set up contrast
  contr_string <- "group_IDsspike"
  contrast <- makeContrasts(contr_string, levels = design)
  
  # calculate differential tests
  res_cydar <- glmQLFTest(fit, contrast = contrast)
  
  # raw p-values
  p_vals <- res_cydar$table$PValue
  
  # controlling the spatial false discovery rate (FDR)
  q_vals <- spatialFDR(intensities(cd), res_cydar$table$PValue)
  
  # significant hyperspheres
  is.sig <- q_vals <= 0.1
  print(table(is.sig))
  
})


# ---------------------------------------------------------------------
# Visualizing and interpreting the results (from Bioconductor vignette)
# ---------------------------------------------------------------------

# 2D PCA plots
if (sum(is.sig) > 0) {
  pdf(file.path(DIR_CYDAR_FILES, "cydar_populations_PCA_supp_all_markers.pdf"))
  
  sig.coords <- intensities(cd)[is.sig, ]
  sig.res <- res_cydar$table[is.sig, ]
  coords <- prcomp(sig.coords)
  
  plotCellLogFC(coords$x[, 1], coords$x[, 2], sig.res$logFC)
  
  dev.off()
}


# 2D tSNE plots (does not work)
# if (sum(is.sig) > 0) {
#   pdf(file.path(DIR_CYDAR_FILES, "cydar_populations_tSNE_supp_all_markers.pdf"))
#   
#   sig.coords <- intensities(cd)[is.sig, ]
#   sig.res <- res_cydar$table[is.sig, ]
#   
#   set.seed(123)
#   out_tsne <- Rtsne(sig.coords, pca = FALSE, verbose = TRUE)
#   tsne_coords <- as.data.frame(out_tsne$Y)
#   colnames(tsne_coords) <- c("tSNE_1", "tSNE_2")
#   
#   plotCellLogFC(tsne_coords[, 1], tsne_coords[, 2], sig.res$logFC)
#   
#   dev.off()
# }


# Median marker intensities of each hypersphere
if (sum(is.sig) > 0) {
  pdf(file.path(DIR_CYDAR_FILES, "cydar_medians_supp_all_markers.pdf"), width = 8, height = 3.5)
  
  par(mfrow = c(2, 6), mar = c(2.1, 1.1, 3.1, 1.1))
  limits <- intensityRanges(cd, p = 0.05)
  all.markers <- rownames(markerData(cd))
  for (i in order(all.markers)) { 
    plotCellIntensity(coords$x[, 1], coords$x[, 2], sig.coords[, i], 
                      irange = limits[, i], main = all.markers[i])
  }
  par(mfrow = c(1, 1))
  
  dev.off()
}


# -------
# Runtime
# -------

runtime_total <- runtime_pre[["elapsed"]] + runtime_count[["elapsed"]] + runtime_test[["elapsed"]]
print(runtime_total)

runtime_cydar_supp_all_markers <- runtime_total




##############################
# Return results at cell level
##############################

# Note: cydar returns q-values at the hypersphere level. Since hyperspheres overlap, the
# q-values are not unique for each cell. To evaluate performance at the cell level, we
# assign a unique q-value to each cell, by selecting the smallest q-value for any
# hypersphere containing that cell.


# number of cells per sample (including spike-in cells)
n_cells <- as.vector(fsApply(d_input, nrow))

# identify B cells (these contain the true differential signal; from both 'base' and
# 'spike' conditions)
is_B_cell <- exprs(pool.ff)[, "B_cell"]
stopifnot(length(is_B_cell) == sum(n_cells))


# get smallest q-value for each cell, across all hyperspheres

stopifnot(length(p_vals) == length(q_vals))

# cell assignments
cells <- cellAssignments(cd)
# 'unpack' indices (e.g. '142183 -142188' means all values from 142183 to 142188)
cells <- unpackIndices(cells)

length(cells)
stopifnot(length(cells) == length(p_vals))

# repeat p/q-values for each hypersphere
cells_rep <- unlist(cells)
p_vals_rep <- rep(p_vals, sapply(cells, length))
q_vals_rep <- rep(q_vals, sapply(cells, length))
stopifnot(length(p_vals_rep) == length(cells_rep), 
          length(q_vals_rep) == length(cells_rep))
# split by cell indices
p_vals_rep_split <- split(p_vals_rep, cells_rep)
q_vals_rep_split <- split(q_vals_rep, cells_rep)
# get minimum p/q-value for each unique cell
cells_p_vals <- sapply(p_vals_rep_split, function(v) min(v))
cells_q_vals <- sapply(q_vals_rep_split, function(v) min(v))

# fill in NAs for missing cells
p_vals_all <- q_vals_all <- rep(NA, sum(n_cells))
names(p_vals_all) <- names(q_vals_all) <- 1:sum(n_cells)
p_vals_all[names(cells_p_vals)] <- cells_p_vals
q_vals_all[names(cells_q_vals)] <- cells_q_vals

stopifnot(length(p_vals_all) == length(q_vals_all))


# set up data frame with results (for marker pS6) and true B-cell status at cell level

res_p_vals <- p_vals_all
res_p_adj <- q_vals_all

# replace NAs to ensure same cells are returned for all methods
res_p_vals[is.na(res_p_vals)] <- 1
res_p_adj[is.na(res_p_adj)] <- 1

stopifnot(length(res_p_vals) == length(res_p_adj), 
          length(res_p_vals) == length(is_B_cell))

res <- data.frame(p_vals = res_p_vals, 
                  q_vals = res_p_adj, 
                  B_cell = is_B_cell)

# store results
out_cydar_supp_all_markers <- res




#####################
# Save output objects
#####################

save(out_cydar_supp_all_markers, runtime_cydar_supp_all_markers, 
     file = file.path(DIR_RDATA, "outputs_BCR_XL_sim_cydar_supp_all_markers.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_BCR_XL_sim_cydar_supp_all_markers.txt"))
sessionInfo()
sink()



