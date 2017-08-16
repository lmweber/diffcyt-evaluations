##########################################################################################
# Script to run methods
# 
# - method: diffcyt-DE-med
# - data set: BCR-XL-sim
# 
# Lukas Weber, August 2017
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)


DIR_BENCHMARK <- "../../../../benchmark_data/BCR_XL_sim/data"
DIR_PLOTS <- "../../../plots/BCR_XL_sim/diffcyt_DE_med"
DIR_RDATA <- "../../../RData/BCR_XL_sim"
DIR_SESSION_INFO <- "../../../session_info/BCR_XL_sim"




#######################
# Loop over simulations
#######################

# to do: add loop

sim <- "sim_full"

DIR_SIM <- file.path(DIR_BENCHMARK, sim)

# set up contrast (to compare 'spike' vs. 'base')
contrasts_list <- list(spike = c(0, 1))




###########
# Load data
###########

# filenames
files <- list.files(DIR_SIM, pattern = "\\.fcs$", full.names = TRUE)
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


# ------------------------------------
# choose markers to use for clustering
# ------------------------------------

cols_clustering <- cols_lineage




##################
# diffcyt pipeline
##################

# --------------------
# pre-processing steps
# --------------------

# prepare data into required format
d_se <- prepareData(d_input, sample_IDs, group_IDs, 
                    cols_markers, cols_clustering, cols_func)

colnames(d_se)[cols_clustering]
colnames(d_se)[cols_func]

# transform data
d_se <- transformData(d_se, cofactor = 5)

# clustering
# (note: clustering all samples together)
seed <- 123
runtime_clustering <- system.time(
  d_se <- generateClusters(d_se, xdim = 30, ydim = 30, seed = seed)
)

runtime_clustering  # ~30 sec (30x30 clusters)

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
## to do: possibly change this function so it returns marker columns only

# calculate ECDFs (slow)
# d_ecdfs <- calcECDFs(d_se)
# 
# dim(d_ecdfs)

# subset marker expression values
# d_vals <- calcSubsetVals(d_se)
# 
# dim(d_vals)


# ---------------------------------------------------------------------------
# test for differential expression (DE) of functional markers within clusters
# ---------------------------------------------------------------------------

# set up design matrix
design <- createDesignMatrix(group_IDs)
design

# set up contrast matrix
contrast <- createContrast(group_IDs, contrast = contrasts_list[[1]])
contrast

# run tests
# - note: include 'patient_IDs' as random effects using limma 'duplicateCorrelation' methodology
path <- file.path(DIR_PLOTS, sim)
runtime <- system.time(
  res <- testDE_med(d_se, d_counts, d_medians, 
                    design, contrast, block_IDs = patient_IDs, 
                    path = path)
)

print(runtime)

# show results
rowData(res)

# sort to show top (most highly significant) clusters first
res_sorted <- rowData(res)[order(rowData(res)$adj.P.Val), ]
print(head(res_sorted, 10))
#View(as.data.frame(res_sorted))

# number of significant DE clusters
print(table(res_sorted$adj.P.Val <= 0.05))




##############################
# Return results at cell level
##############################

# Note: diffcyt methods return results at cluster level (e.g. 900 small clusters), for
# each functional marker. To enable performance comparisons between methods at the cell
# level, we assign the cluster-level p-values to all cells within each cluster, for each
# functional marker.


# number of cells per sample
n_cells <- sapply(d_input, nrow)

# spike-in status for each cell
is_spikein <- unlist(sapply(d_input, function(d) exprs(d)[, "is_spikein"]))
stopifnot(length(is_spikein) == sum(n_cells))


# match cluster-level p-values to individual cells

stopifnot(nrow(rowData(res)) == nlevels(rowData(res)$cluster) * nlevels(rowData(res)$marker))
stopifnot(all(levels(rowData(d_se)$cluster) == levels(rowData(res)$cluster)))

# split by marker
row_data_res_split <- split(rowData(res), rowData(res)$marker)

# check markers are in correct order after split
stopifnot(names(row_data_res_split) == unique(rowData(res)$marker))

p_vals_cells <- p_adj_cells <- vector("list", length(row_data_res_split))
names(p_vals_cells) <- names(p_adj_cells) <- names(row_data_res_split)

for (i in seq_along(row_data_res_split)) {
  ix_match <- match(rowData(d_se)$cluster, row_data_res_split[[i]]$cluster)
  stopifnot(length(ix_match) == nrow(rowData(d_se)))
  
  p_vals_clusters <- row_data_res_split[[i]]$P.Value
  p_adj_clusters <- row_data_res_split[[i]]$adj.P.Val
  
  p_vals_cells[[i]] <- p_vals_clusters[ix_match]
  p_adj_cells[[i]] <- p_adj_clusters[ix_match]
  
  # replace NAs (due to filtering) to ensure same cells are returned for all methods
  p_vals_cells[[i]][is.na(p_vals_cells[[i]])] <- 1
  p_adj_cells[[i]][is.na(p_adj_cells[[i]])] <- 1
}


# set up data frames with results and true spike-in status at cell level

res <- vector("list", length(row_data_res_split))
names(res) <- names(row_data_res_split)

for (i in seq_along(res)) {
  stopifnot(length(p_vals_cells[[i]]) == length(is_spikein), 
            length(p_adj_cells[[i]]) == length(is_spikein))
  
  res[[i]] <- data.frame(p_vals = p_vals_cells[[i]], 
                         p_adj = p_adj_cells[[i]], 
                         spikein = is_spikein)
}

# store results
out_diffcyt_DE_med <- res




#####################
# Save output objects
#####################

save(out_diffcyt_DE_med, file = file.path(DIR_RDATA, "/outputs_BCR_XL_sim_diffcyt_DE_med.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "/session_info_BCR_XL_sim_diffcyt_DE_med.txt"))
sessionInfo()
sink()



