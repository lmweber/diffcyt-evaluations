##########################################################################################
# Script to run methods
# 
# - method: diffcyt-DA-voom
# - data set: BCR-XL-sim
# 
# - supplementary results: treating all markers as 'cell type' markers; using DA methods
# 
# Lukas Weber, November 2018
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)


DIR_BENCHMARK <- "../../../../../benchmark_data/BCR_XL_sim/data/main"
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/supp_cell_type_markers/diagnostic/diffcyt_DA_voom"
DIR_RDATA <- "../../../../RData/BCR_XL_sim/supp_cell_type_markers"
DIR_SESSION_INFO <- "../../../../session_info/BCR_XL_sim/supp_cell_type_markers"




###########################
# Load data, pre-processing
###########################

# filenames

files <- list.files(DIR_BENCHMARK, pattern = "\\.fcs$", full.names = TRUE)
files_base <- files[grep("base\\.fcs$", files)]
files_spike <- files[grep("spike\\.fcs$", files)]

files_load <- c(files_base, files_spike)
files_load

# load data

d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)

# sample information

sample_id <- gsub("^BCR_XL_sim_", "", 
                  gsub("\\.fcs$", "", basename(files_load)))
sample_id

group_id <- factor(gsub("^.*_", "", sample_id), levels = c("base", "spike"))
group_id

patient_id <- factor(gsub("_.*$", "", sample_id))
patient_id

experiment_info <- data.frame(group_id, patient_id, sample_id)
experiment_info

# marker information

# indices of all marker columns, lineage markers, and functional markers
# (10 surface markers / 14 functional markers; see Bruggner et al. 2014, Table 1)
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)

marker_name <- colnames(d_input[[1]])
marker_name <- gsub("\\(.*$", "", marker_name)

marker_class <- rep("none", length(marker_name))
marker_class[cols_lineage] <- "type"
marker_class[cols_func] <- "state"
marker_class <- factor(marker_class, levels = c("type", "state", "none"))

# exclude CD45 from clustering
marker_class[marker_name == "CD45"] <- "none"

marker_info <- data.frame(marker_name, marker_class)
marker_info


# supplementary results: treat all markers as 'cell type' markers; then use DA methods

marker_info$marker_class[marker_info$marker_class == "state"] <- "type"
marker_info




##################
# diffcyt pipeline
##################

# --------------------
# pre-processing steps
# --------------------

runtime_preprocessing <- system.time({
  
  # prepare data into required format
  d_se <- prepareData(d_input, experiment_info, marker_info)
  
  colnames(d_se)[colData(d_se)$marker_class == "type"]
  colnames(d_se)[colData(d_se)$marker_class == "state"]
  
  # transform data
  d_se <- transformData(d_se, cofactor = 5)
  
  # clustering
  # (runtime: ~5 sec with xdim = 10, ydim = 10)
  seed <- 1234
  d_se <- generateClusters(d_se, xdim = 10, ydim = 10, seed_clustering = seed)
  
  length(table(rowData(d_se)$cluster_id))  # number of clusters
  nrow(rowData(d_se))                      # number of cells
  sum(table(rowData(d_se)$cluster_id))
  min(table(rowData(d_se)$cluster_id))     # size of smallest cluster
  max(table(rowData(d_se)$cluster_id))     # size of largest cluster
  
  # calculate cluster cell counts
  d_counts <- calcCounts(d_se)
  
  dim(d_counts)
  rowData(d_counts)
  length(assays(d_counts))
  
  # calculate cluster medians
  d_medians <- calcMedians(d_se)
  
  dim(d_medians)
  rowData(d_medians)
  length(assays(d_medians))
  names(assays(d_medians))
  
  # calculate medians by cluster and marker
  d_medians_by_cluster_marker <- calcMediansByClusterMarker(d_se)
  
  dim(d_medians_by_cluster_marker)
  length(assays(d_medians_by_cluster_marker))
  
  # calculate medians by sample and marker
  d_medians_by_sample_marker <- calcMediansBySampleMarker(d_se)
  
  dim(d_medians_by_sample_marker)
  length(assays(d_medians_by_sample_marker))
  
})


# ---------------------------------
# store data objects (for plotting)
# ---------------------------------

out_objects_diffcyt_DA_voom_supp_cell_type_markers <- list(
  d_se = d_se, 
  d_counts = d_counts, 
  d_medians = d_medians, 
  d_medians_by_cluster_marker = d_medians_by_cluster_marker, 
  d_medians_by_sample_marker = d_medians_by_sample_marker
)


# -----------------------------------------
# test for differentially abundant clusters
# -----------------------------------------

# contrast (to compare 'spike' vs. 'base')
# note: include fixed effects for 'patient_id'
contrast_vec <- c(0, 1, 0, 0, 0, 0, 0, 0, 0)

runtime_tests <- system.time({
  
  # set up design matrix
  # note: include fixed effects for 'patient_id'
  design <- createDesignMatrix(experiment_info, cols_design = 1:2)
  design
  
  # set up contrast matrix
  contrast <- createContrast(contrast_vec)
  contrast
  
  # run tests
  path <- DIR_PLOTS
  res <- testDA_voom(d_counts, design, contrast, path = path)
  
})

# show results
rowData(res)

# sort to show top (most highly significant) clusters first
res_sorted <- rowData(res)[order(rowData(res)$p_adj), ]
print(head(res_sorted, 10))
#View(as.data.frame(res_sorted))

# number of significant tests (note: one test per cluster)
print(table(res_sorted$p_adj <= 0.1))

# runtime
runtime_total <- runtime_preprocessing[["elapsed"]] + runtime_tests[["elapsed"]]
print(runtime_total)

runtime_diffcyt_DA_voom_supp_cell_type_markers <- runtime_total


# ---------------------------------------------
# store results at cluster level (for plotting)
# ---------------------------------------------

res_clusters <- as.data.frame(rowData(res))

out_clusters_diffcyt_DA_voom_supp_cell_type_markers <- res_clusters




##############################
# Return results at cell level
##############################

# Note: diffcyt methods return results for each cluster. To enable performance
# comparisons between methods at the cell level, we assign the cluster-level p-values
# to all cells within each cluster.


# number of cells per sample (including spike-in cells)
n_cells <- sapply(d_input, nrow)

# identify B cells (these contain the true differential signal; from both 'base' and
# 'spike' conditions)
is_B_cell <- unlist(sapply(d_input, function(d) exprs(d)[, "B_cell"]))
stopifnot(length(is_B_cell) == sum(n_cells))


# match cluster-level p-values to individual cells

stopifnot(nrow(rowData(res)) == length(levels(rowData(d_se)$cluster_id)), 
          all(rowData(res)$cluster_id == levels(rowData(d_se)$cluster_id)))

rowData(res)$cluster_id <- factor(rowData(res)$cluster_id, levels = levels(rowData(d_se)$cluster_id))

# match cells to clusters
ix_match <- match(rowData(d_se)$cluster_id, rowData(res)$cluster_id)

p_vals_clusters <- rowData(res)$p_val
p_adj_clusters <- rowData(res)$p_adj

p_vals_cells <- p_vals_clusters[ix_match]
p_adj_cells <- p_adj_clusters[ix_match]


# set up data frame with results and true B-cell status at cell level

res_p_vals <- p_vals_cells
res_p_adj <- p_adj_cells

# replace NAs (due to filtering) to ensure same cells are returned for all methods
res_p_vals[is.na(res_p_vals)] <- 1
res_p_adj[is.na(res_p_adj)] <- 1

# return values for this condition and healthy
res <- data.frame(p_val = res_p_vals, 
                  p_adj = res_p_adj, 
                  B_cell = is_B_cell)

# store results
out_diffcyt_DA_voom_supp_cell_type_markers <- res




#####################
# Save output objects
#####################

save(out_diffcyt_DA_voom_supp_cell_type_markers, runtime_diffcyt_DA_voom_supp_cell_type_markers, 
     file = file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DA_voom_supp_cell_type_markers.RData"))

save(out_clusters_diffcyt_DA_voom_supp_cell_type_markers, 
     file = file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DA_voom_supp_cell_type_markers.RData"))

save(out_objects_diffcyt_DA_voom_supp_cell_type_markers, 
     file = file.path(DIR_RDATA, "out_objects_BCR_XL_sim_diffcyt_DA_voom_supp_cell_type_markers.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_BCR_XL_sim_diffcyt_DA_voom_supp_cell_type_markers.txt"))
sessionInfo()
sink()



