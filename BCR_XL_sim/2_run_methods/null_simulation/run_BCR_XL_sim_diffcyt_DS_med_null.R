##########################################################################################
# Script to run methods
# 
# - method: diffcyt-DS-med
# - data set: BCR-XL-sim
# 
# - null simulation
# 
# Lukas Weber, November 2017
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)


DIR_BENCHMARK <- "../../../../../benchmark_data/BCR_XL_sim/data/null_simulation"
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/null_simulation_diagnostic/diffcyt_DS_med"
DIR_RDATA <- "../../../../RData/BCR_XL_sim/null_simulation"
DIR_SESSION_INFO <- "../../../../session_info/BCR_XL_sim/null_simulation"




#############
# Preliminary
#############

# contrast (to compare 'spike' vs. 'base')
# note: include zeros for patient_IDs fixed effects
contrasts_list <- list(spike = c(0, 1, 0, 0, 0, 0, 0, 0, 0))




###########################
# Load data, pre-processing
###########################

# filenames
files <- list.files(DIR_BENCHMARK, pattern = "\\.fcs$", full.names = TRUE)
files_base <- files[grep("base\\.fcs$", files)]
files_spike <- files[grep("spike\\.fcs$", files)]

# load data
files_load <- c(files_base, files_spike)
files_load

d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)

# sample IDs, group IDs, patient IDs
sample_IDs <- gsub("^BCR_XL_sim_", "", 
                   gsub("\\.fcs$", "", basename(files_load)))
sample_IDs

group_IDs <- factor(gsub("^.*_", "", sample_IDs), levels = c("base", "spike"))
group_IDs

patient_IDs <- factor(gsub("_.*$", "", sample_IDs))
patient_IDs

# check
data.frame(sample_IDs, group_IDs, patient_IDs)

# indices of all marker columns, lineage markers, and functional markers
# (10 surface markers / 14 functional markers; see Bruggner et al. 2014, Table 1)
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)




##################
# diffcyt pipeline
##################

# --------------------
# pre-processing steps
# --------------------

runtime_preprocessing <- system.time({
  
  # prepare data into required format
  d_se <- prepareData(d_input, sample_IDs, group_IDs, 
                      cols_markers, cols_lineage, cols_func)
  
  colnames(d_se)[cols_lineage]
  colnames(d_se)[cols_func]
  
  # transform data
  d_se <- transformData(d_se, cofactor = 5)
  
  # clustering
  # (runtime: ~5 sec with xdim = 10, ydim = 10)
  seed <- 100
  d_se <- generateClusters(d_se, xdim = 10, ydim = 10, seed = seed)
  
  length(table(rowData(d_se)$cluster))  # number of clusters
  nrow(rowData(d_se))                   # number of cells
  sum(table(rowData(d_se)$cluster))
  min(table(rowData(d_se)$cluster))     # size of smallest cluster
  max(table(rowData(d_se)$cluster))     # size of largest cluster
  
  # calculate cluster cell counts
  d_counts <- calcCounts(d_se)
  
  dim(d_counts)
  rowData(d_counts)
  length(assays(d_counts))
  
  # calculate cluster medians by sample
  d_medians <- calcMedians(d_se)
  
  dim(d_medians)
  rowData(d_medians)
  length(assays(d_medians))
  names(assays(d_medians))
  
  # calculate cluster medians across all samples
  d_medians_all <- calcMediansAll(d_se)
  
  dim(d_medians_all)
  length(assays(d_medians_all))
  
  # calculate ECDFs (runtime: ~10 sec)
  d_ecdfs <- calcECDFs(d_se)
  
  dim(d_ecdfs)
  length(assays(d_ecdfs))
  
  # subset marker expression values
  d_vals <- calcSubsetVals(d_se)
  
  dim(d_vals)
  length(assays(d_vals))
  
})


# ---------------------------------
# store data objects (for plotting)
# ---------------------------------

out_objects_diffcyt_DS_med_null <- list(d_se = d_se, 
                                        d_counts = d_counts, 
                                        d_medians = d_medians, 
                                        d_medians_all = d_medians_all)


# -------------------------------------------------------
# test for differential functional states within clusters
# -------------------------------------------------------

runtime_tests <- system.time({
  
  # set up design matrix
  # note: include 'patient_IDs' as fixed effects ('block_IDs' argument)
  design <- createDesignMatrix(group_IDs, block_IDs = patient_IDs)
  design
  
  # set up contrast matrix
  contrast <- createContrast(group_IDs, contrast = contrasts_list$spike)
  contrast
  
  # run tests
  # note: including 'patient_IDs' as fixed effects ('block_IDs' argument)
  res <- testDS_med(d_counts, d_medians, design, contrast, 
                    block_IDs = patient_IDs, path = DIR_PLOTS)
  
})

# show results
rowData(res)

# sort to show top (most highly significant) cluster-marker combinations first
res_sorted <- rowData(res)[order(rowData(res)$adj.P.Val), ]
print(head(res_sorted, 10))
#View(as.data.frame(res_sorted))

# number of significant tests (note: one test per cluster-marker combination)
print(table(res_sorted$adj.P.Val <= 0.05))

# runtime (~1 min on laptop)
runtime_total <- runtime_preprocessing[["elapsed"]] + runtime_tests[["elapsed"]]
print(runtime_total)

runtime_diffcyt_DS_med_null <- runtime_total


# ---------------------------------------------
# store results at cluster level (for plotting)
# ---------------------------------------------

res_clusters <- as.data.frame(rowData(res))

out_clusters_diffcyt_DS_med_null <- res_clusters





##############################
# Return results at cell level
##############################

# Note: diffcyt methods return results for each cluster-marker combination. To enable
# performance comparisons between methods at the cell level, we assign the same p-values
# to all cells within a given cluster-marker combination.

# Note: return cell-level results for marker pS6 only, since the comparative evaluations
# (ROC curves etc) are based on pS6 in B cells only.


# number of cells per sample (including spike-in cells)
n_cells <- sapply(d_input, nrow)

# spike-in status for each cell
is_spikein <- unlist(sapply(d_input, function(d) exprs(d)[, "spikein"]))
stopifnot(length(is_spikein) == sum(n_cells))


# match cluster-level p-values for marker pS6 to individual cells

stopifnot(nrow(rowData(res)) == nlevels(rowData(d_se)$cluster) * length(cols_func), 
          all(levels(rowData(res)$cluster) == levels(rowData(d_se)$cluster)), 
          all(levels(rowData(res)$cluster) %in% rowData(res)$cluster))

# select results for pS6
res_pS6 <- res[rowData(res)$marker == "pS6(Yb172)Dd", ]

# match cells to clusters
ix_match <- match(rowData(d_se)$cluster, rowData(res_pS6)$cluster)

p_vals_clusters <- rowData(res_pS6)$P.Value
p_adj_clusters <- rowData(res_pS6)$adj.P.Val

p_vals_cells <- p_vals_clusters[ix_match]
p_adj_cells <- p_adj_clusters[ix_match]


# set up data frame with results (for marker pS6) and true spike-in status at cell level

# select results from spike-in samples only
which_keep <- rowData(d_se)$group == "spike"

is_spikein_keep <- is_spikein[which_keep]

res_p_vals <- p_vals_cells[which_keep]
res_p_adj <- p_adj_cells[which_keep]

# replace NAs (due to filtering) to ensure same cells are returned for all methods
res_p_vals[is.na(res_p_vals)] <- 1
res_p_adj[is.na(res_p_adj)] <- 1

stopifnot(length(res_p_vals) == length(res_p_adj), 
          length(res_p_vals) == length(is_spikein_keep))

res <- data.frame(p_vals = res_p_vals, 
                  p_adj = res_p_adj, 
                  spikein = is_spikein_keep)

# store results
out_diffcyt_DS_med_null <- res




#####################
# Save output objects
#####################

save(out_diffcyt_DS_med_null, runtime_diffcyt_DS_med_null, 
     file = file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_med_null.RData"))

save(out_clusters_diffcyt_DS_med_null, 
     file = file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DS_med_null.RData"))

save(out_objects_diffcyt_DS_med_null, 
     file = file.path(DIR_RDATA, "out_objects_BCR_XL_sim_diffcyt_DS_med_null.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_BCR_XL_sim_diffcyt_DS_med_null.txt"))
sessionInfo()
sink()



