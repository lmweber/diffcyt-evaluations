##########################################################################################
# Script to run methods
# 
# - method: diffcyt-DS-LMM
# - data set: BCR-XL
# 
# - main results
# 
# Lukas Weber, May 2018
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)


DIR_BENCHMARK <- "../../../../../benchmark_data/BCR_XL"
DIR_RDATA <- "../../../../RData/BCR_XL/main"
DIR_SESSION_INFO <- "../../../../session_info/BCR_XL/main"




#########################
# Load data and meta-data
#########################

# ---------------
# Load .fcs files
# ---------------

DIR_FCS_FILES <- "Bodenmiller_BCR_XL_fcs_files"

files_load_fcs <- list.files(file.path(DIR_BENCHMARK, DIR_FCS_FILES), pattern = "\\.fcs$", full.names = TRUE)

# load as flowSet
d_flowSet <- read.flowSet(files_load_fcs, transformation = FALSE, truncate_max_range = FALSE)


# ----------------------
# Load population labels
# ----------------------

DIR_POP_IDS <- "Bodenmiller_BCR_XL_population_IDs"

files_load_pop <- list.files(file.path(DIR_BENCHMARK, DIR_POP_IDS), pattern = "\\.csv$", full.names = TRUE)

d_population_IDs <- lapply(files_load_pop, read.csv)

# check numbers of cells match
stopifnot(all(sapply(seq_along(files_load_fcs), function(i) {
  nrow(d_population_IDs[[i]]) == nrow(d_flowSet[[i]])
})))


# ----------------
# Create meta-data
# ----------------

# sample information

# check sample order
stopifnot(all(pData(d_flowSet)$name == basename(files_load_fcs)))

# sample information
sample_IDs <- gsub("^PBMC8_30min_", "", gsub("\\.fcs$", "", basename(files_load_fcs)))
group_IDs <- factor(gsub("^patient[0-9]+_", "", sample_IDs), levels = c("Reference", "BCR-XL"))
patient_IDs <- factor(gsub("_.*$", "", sample_IDs))

experiment_info <- data.frame(group_id = group_IDs, patient_id = patient_IDs, sample_id = sample_IDs, 
                              stringsAsFactors = FALSE)
experiment_info


# marker information

# indices of all marker columns, lineage markers, and functional markers
# (10 surface markers / 14 functional markers; see Bruggner et al. 2014, Table 1)
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)

channel_name <- colnames(d_flowSet)

marker_name <- gsub("\\(.*$", "", channel_name)

marker_class <- rep("none", ncol(d_flowSet[[1]]))
marker_class[cols_lineage] <- "type"
marker_class[cols_func] <- "state"
marker_class <- factor(marker_class, levels = c("type", "state", "none"))

# exclude CD45 from clustering
marker_class[marker_name == "CD45"] <- "none"

marker_info <- data.frame(marker_name, marker_class, stringsAsFactors = FALSE)
marker_info




##################
# diffcyt pipeline
##################

# --------------------
# pre-processing steps
# --------------------

runtime_preprocessing <- system.time({
  
  # prepare data into required format
  d_se <- prepareData(d_flowSet, experiment_info, marker_info)
  
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

out_objects_diffcyt_DS_LMM_main <- list(
  d_se = d_se, 
  d_counts = d_counts, 
  d_medians = d_medians, 
  d_medians_by_cluster_marker = d_medians_by_cluster_marker, 
  d_medians_by_sample_marker = d_medians_by_sample_marker
)


# --------------------------------------------
# test for differential states within clusters
# --------------------------------------------

# contrast (to compare 'BCR-XL' vs. 'Reference')
# note: include random effects for 'patient_id'
contrast_vec <- c(0, 1)

runtime_tests <- system.time({
  
  # set up model formula
  # note: include random effects for 'patient_id'
  formula <- createFormula(experiment_info, cols_fixed = 1, cols_random = 2)
  formula
  
  # set up contrast matrix
  contrast <- createContrast(contrast_vec)
  contrast
  
  # run tests
  res <- testDS_LMM(d_counts, d_medians, formula, contrast)
  
})

# show results
rowData(res)

# sort to show top (most highly significant) cluster-marker combinations first
res_sorted <- rowData(res)[order(rowData(res)$p_adj), ]
print(head(res_sorted, 10))
#View(as.data.frame(res_sorted))

# number of significant tests (note: one test per cluster-marker combination)
print(table(res_sorted$p_adj <= 0.1))

# runtime
runtime_total <- runtime_preprocessing[["elapsed"]] + runtime_tests[["elapsed"]]
print(runtime_total)

runtime_diffcyt_DS_LMM_main <- runtime_total


# ---------------------------------------------
# store results at cluster level (for plotting)
# ---------------------------------------------

res_clusters <- as.data.frame(rowData(res))

out_clusters_diffcyt_DS_LMM_main <- res_clusters




######################
# Store population IDs
######################

# add population IDs to 'rowData' of 'd_se' object

n_cells <- sapply(as(d_flowSet, "list"), nrow)
n_cells

# check sample order
stopifnot(all(names(n_cells) == basename(files_load_fcs)))
rowData(d_se)[(n_cells[1] - 1):(n_cells[1] + 2), ]

population_IDs_rep <- unlist(do.call("rbind", d_population_IDs))

# add to 'rowData'
rowData(d_se)$population_id <- population_IDs_rep

# store object
out_objects_diffcyt_DS_LMM_main$d_se <- d_se




#####################
# Save output objects
#####################

save(runtime_diffcyt_DS_LMM_main, 
     file = file.path(DIR_RDATA, "outputs_BCR_XL_diffcyt_DS_LMM_main.RData"))

save(out_clusters_diffcyt_DS_LMM_main, 
     file = file.path(DIR_RDATA, "out_clusters_BCR_XL_diffcyt_DS_LMM_main.RData"))

save(out_objects_diffcyt_DS_LMM_main, 
     file = file.path(DIR_RDATA, "out_objects_BCR_XL_diffcyt_DS_LMM_main.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_BCR_XL_diffcyt_DS_LMM_main.txt"))
sessionInfo()
sink()



