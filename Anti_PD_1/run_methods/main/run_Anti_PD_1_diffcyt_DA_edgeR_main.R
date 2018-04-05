##########################################################################################
# Script to run methods
# 
# - method: diffcyt-DA-edgeR
# - data set: Anti-PD-1
# 
# - main results
# 
# Lukas Weber, April 2018
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)
library(readxl)


DIR_BENCHMARK <- "../../../../../benchmark_data/Anti_PD_1"
DIR_RDATA <- "../../../../RData/Anti_PD_1/main"
DIR_SESSION_INFO <- "../../../../session_info/Anti_PD_1/main"




#########################
# Load data and meta-data
#########################

# note: combine data from batches '23' and '29' (see analysis for Krieg et al., 2018,
# using CellCnn: https://github.com/lmweber/PD1_analysis_CellCnn)


fn_metadata_23 <- file.path(DIR_BENCHMARK, "CK_metadata", "metadata_23_03all.xlsx")
fn_metadata_29 <- file.path(DIR_BENCHMARK, "CK_metadata", "metadata_29_03all3.xlsx")

path_23 <- file.path(DIR_BENCHMARK, "CK_2016-06-23_03all", "010_cleanfcs")
path_29 <- file.path(DIR_BENCHMARK, "CK_2016-06-29_03all3", "010_cleanfcs")

fn_panel <- file.path(DIR_BENCHMARK, "CK_panels", "panel3_v3.xlsx")


# -------------
# load metadata
# -------------

# load metadata spreadsheets for each data set ("data 23" and "data 29")
metadata_23 <- read_excel(fn_metadata_23)
metadata_29 <- read_excel(fn_metadata_29)

#View(metadata_23)
#View(metadata_29)

ix_keep <- 6:15

paths <- c(rep(path_23, length(ix_keep)), rep(path_29, length(ix_keep)))

files <- c(metadata_23$filename[ix_keep], metadata_29$filename[ix_keep])


# -----------
# sample info
# -----------

# group IDs
group_IDs <- factor(gsub("^base_", "", c(metadata_23$condition[ix_keep], metadata_29$condition[ix_keep])), 
                    levels = c("NR", "R"))
group_IDs

# batch (data set) IDs
batch_IDs <- factor(c(rep("batch23", length(ix_keep)), rep("batch29", length(ix_keep))), 
                    levels = c("batch23", "batch29"))
batch_IDs

# sample IDs
sample_IDs <- factor(gsub("^base_", "", c(metadata_23$shortname[ix_keep], metadata_29$shortname[ix_keep])), 
                     levels = c(paste0("NR", 1:9), paste0("R", 1:11)))
sample_IDs

sample_info <- data.frame(group = group_IDs, batch = batch_IDs, sample = sample_IDs)
sample_info


# ----------------
# load data into R
# ----------------

stopifnot(length(paths) == length(files))

fn <- file.path(paths, files)

d_input <- lapply(fn, read.FCS, transformation = FALSE, truncate_max_range = FALSE)

# check column names match (note: exclude columns 58 and 59, which contain "beadDist" and "Time")
ix <- 1:57
check_cols <- lapply(d_input, function(d) pData(parameters(d))$name)
all(sapply(check_cols, function(ch) all(ch[ix] == check_cols[[1]][ix])))


# -----------
# marker info
# -----------

# load panel details spreadsheet
panel <- as.data.frame(read_excel(fn_panel))
panel

# update panel to include CD45
panel$transform[panel$fcs_colname == "Y89Di"] <- 1
panel

# replace NAs
panel$Antigen[is.na(panel$Antigen)] <- panel$fcs_colname[is.na(panel$Antigen)]
panel


# subset and rearrange data; so order of columns matches rows in panel

markers_ix <- match(panel$fcs_colname, pData(parameters(d_input[[1]]))$name)

stopifnot(length(markers_ix) == nrow(panel))

d_input <- lapply(d_input, function(d) {
  e <- exprs(d)[, markers_ix]
  colnames(e) <- panel$Antigen
  flowFrame(e)
})


# marker information
# (note: all surface markers in this panel)

is_marker <- as.logical(panel$transform)

marker_type <- rep("none", nrow(panel))
marker_type[is_marker] <- "cell_type"
marker_type <- factor(marker_type, levels = c("cell_type", "cell_state", "none"))

marker_name <- panel$Antigen

marker_info <- data.frame(marker_name, is_marker, marker_type)
marker_info




##################
# diffcyt pipeline
##################

# --------------------
# pre-processing steps
# --------------------

# random seed
seed <- 10000


runtime_preprocessing <- system.time({
  
  # prepare data into required format
  d_se <- prepareData(d_input, sample_info, marker_info)
  
  colnames(d_se)[colData(d_se)$marker_type == "cell_type"]
  colnames(d_se)[colData(d_se)$marker_type == "cell_state"]
  
  # transform data
  d_se <- transformData(d_se, cofactor = 5)
  
  # clustering
  # (runtime: ~20 sec with xdim = 20, ydim = 20)
  d_se <- generateClusters(d_se, xdim = 20, ydim = 20, seed = seed)
  
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

out_objects_diffcyt_DA_edgeR_main <- list(
  d_se = d_se, 
  d_counts = d_counts, 
  d_medians = d_medians, 
  d_medians_by_cluster_marker = d_medians_by_cluster_marker, 
  d_medians_by_sample_marker = d_medians_by_sample_marker
)


# -----------------------------------------
# test for differentially abundant clusters
# -----------------------------------------

runtime_test <- system.time({
  
  # set up design matrix
  # note: include fixed effects for 'batch'
  design <- createDesignMatrix(sample_info, cols_include = 1:2)
  design
  
  # set up contrast matrix
  contrast_vec <- c(0, 1, 0)
  contrast <- createContrast(contrast_vec)
  contrast
  
  # run tests
  # note: adjust filtering parameters
  min_cells <- 3
  min_samples <- min(table(sample_info$group))
  res <- testDA_edgeR(d_counts, design, contrast, 
                      min_cells = min_cells, min_samples = min_samples)
  
})

# show results
rowData(res)

# sort to show top (most highly significant) clusters first
res_sorted <- rowData(res)[order(rowData(res)$FDR), ]
print(head(res_sorted, 10))
#View(as.data.frame(res_sorted))

hist(res_sorted$PValue)
hist(res_sorted$FDR)

# number of significant tests (note: one test per cluster)
print(table(res_sorted$FDR <= 0.2))

# runtime
runtime_total <- runtime_preprocessing[["elapsed"]] + runtime_test[["elapsed"]]
print(runtime_total)

runtime_diffcyt_DA_edgeR_main <- runtime_total


# -----------------------------
# save results at cluster level
# -----------------------------

# note: do not need results at cell level, since this is not simulated data

res_clusters <- as.data.frame(rowData(res))

out_clusters_diffcyt_DA_edgeR_main <- res_clusters




#####################
# Save output objects
#####################

# note: do not need results at cell level, since this is not simulated data

save(runtime_diffcyt_DA_edgeR_main, 
     file = file.path(DIR_RDATA, "outputs_Anti_PD_1_diffcyt_DA_edgeR_main.RData"))

save(out_clusters_diffcyt_DA_edgeR_main, 
     file = file.path(DIR_RDATA, "out_clusters_Anti_PD_1_diffcyt_DA_edgeR_main.RData"))

save(out_objects_diffcyt_DA_edgeR_main, 
     file = file.path(DIR_RDATA, "out_objects_Anti_PD_1_diffcyt_DA_edgeR_main.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_Anti_PD_1_diffcyt_DA_edgeR_main.txt"))
sessionInfo()
sink()



