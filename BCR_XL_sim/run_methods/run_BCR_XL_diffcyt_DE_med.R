##########################################################################################
# Script to run methods
# 
# - method: diffcyt-DE-med
# - data set: BCR-XL
# 
# Lukas Weber, August 2017
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)
library(Rtsne)
library(magrittr)
library(dplyr)
library(reshape2)
library(ggplot2)


DIR_BENCHMARK <- "../../../../benchmark_data/BCR_XL/data"
DIR_PLOTS <- ""
DIR_RDATA <- ""
DIR_SESSION_INFO <- ""




######################################
# Load data, pre-processing, transform
######################################

# ---------
# load data
# ---------

# filenames
files <- list.files(DIR_BENCHMARK, pattern = "\\.fcs$", full.names = TRUE)
files_BCRXL <- files[grep("BCR-XL", files)]
files_ref <- files[grep("Reference", files)]

# load data
files_load <- c(files_BCRXL, files_ref)
files_load

d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)

# sample IDs, group IDs, patient IDs
sample_IDs <- gsub("^PBMC8_30min_", "", gsub("\\.fcs$", "", basename(files_load)))
sample_IDs

group_IDs <- gsub("^patient[0-9]_", "", sample_IDs)
group_IDs

patient_IDs <- factor(gsub("_.*$", "", 
                           gsub("^patient", "", sample_IDs)))
patient_IDs

# set group_IDs reference level (for differential tests)
group_IDs <- factor(group_IDs, levels = c("Reference", "BCR-XL"))
group_IDs

# check all match correctly
data.frame(sample_IDs, group_IDs, patient_IDs)

# indices of all marker columns, lineage markers, and functional markers
# (10 surface markers / 14 functional markers; see Bruggner et al. 2014, Table 1)
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)


# ------------------------------------
# choose markers to use for clustering
# ------------------------------------

cols_to_use <- cols_lineage




##################
# diffcyt pipeline
##################

# --------------------
# pre-processing steps
# --------------------

# prepare data into required format
d_se <- prepareData(d_input, sample_IDs, group_IDs, 
                    cols_markers, cols_to_use, cols_func)

colnames(d_se)[cols_to_use]
colnames(d_se)[cols_func]

# transform data
d_se <- transformData(d_se, cofactor = 5)

# clustering
# (note: clustering all samples together)
seed <- 123
runtime_clustering <- system.time(
  d_se <- generateClusters(d_se, xdim = 30, ydim = 30, seed = seed)
)

runtime_clustering  # ~60 sec (30x30 clusters)

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
colData(d_medians)
length(assays(d_medians))
names(assays(d_medians))

# calculate cluster medians across all samples
d_medians_all <- calcMediansAll(d_se)

dim(d_medians_all)

# calculate ECDFs
d_ecdfs <- calcECDFs(d_se)

dim(d_ecdfs)

# subset marker expression values
d_vals <- calcSubsetVals(d_se)

dim(d_vals)


# ------------------
# visualize clusters
# ------------------

# t-SNE plot

n_cells <- rowData(d_counts)$n_cells

d_plot <- data.frame(n_cells = n_cells)

# data for Rtsne
d_tsne <- assay(d_medians_all)[, colData(d_medians_all)$is_clustering_col]
d_tsne <- as.matrix(d_tsne)

# remove any duplicate rows (required by Rtsne)
dups <- duplicated(d_tsne)
d_tsne <- d_tsne[!dups, ]

# also remove duplicated rows from plotting data
d_plot <- d_plot[!dups, ]

# run Rtsne
# (note: initial PCA step not required, since we do not have too many dimensions)
set.seed(1234)
out_tsne <- Rtsne(d_tsne, pca = FALSE, verbose = TRUE)

tsne_coords <- as.data.frame(out_tsne$Y)
colnames(tsne_coords) <- c("tSNE_1", "tSNE_2")

d_plot <- cbind(d_plot, tsne_coords)

# plot
ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, size = n_cells)) + 
  geom_point(alpha = 0.75) + 
  scale_size_continuous(range = c(0.25, 4)) + 
  coord_fixed() + 
  ggtitle(paste0("t-SNE map: BCR-XL data set")) + 
  theme_bw()

#path <- "../../../plots/BCR_XL/diffcyt_DE_med"
#filename <- file.path(path, paste0("results_BCR_XL_diffcyt_clustering_tSNE.pdf"))

#ggsave(filename, width = 9, height = 9)



## use 20 clusters + merging strategy from Gosia's paper instead

## follow workflow up to about page 25 to get merged clusters and cell population IDs

## then vary the amount of pS6 signal between BCR-XL and Reference groups in B cells and monocytes only
## (and also demonstrate that the detected small clusters are indeed B cells and monocytes)

## or: spike in activated B cells from the BCR-XL samples into the Reference samples
## (and split the Reference samples into two halves: one with spike-in and one without)


## then: for the third data set (Stephane cancer data), just show some qualitative results


## for each data set: first show qualitative results, then simulations (and no simulations for third data set)




# ---------------------------------------------------------------------------
# test for differential expression (DE) of functional markers within clusters
# ---------------------------------------------------------------------------



