##########################################################################################
# Script to run 'diffcyt' methods for data set 'BCR-XL-all'
#
# Lukas Weber, May 2017
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)



##############################
# Load data and pre-processing
##############################

# ---------
# load data
# ---------

# filenames
files <- list.files("../../../../benchmark_data/BCR_XL/data", pattern = "\\.fcs$", full.names = TRUE)
files_BCRXL <- files[grep("BCR-XL", files)]
files_ref <- files[grep("Reference", files)]

# load data
files_load <- c(files_BCRXL, files_ref)
files_load

d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)

# sample IDs and group IDs
sample_IDs <- gsub("^PBMC8_30min_", "", gsub("\\.fcs$", "", basename(files_load)))
sample_IDs
group_IDs <- gsub("^patient[0-9]_", "", sample_IDs)
group_IDs

# set group reference level for differential testing
group_IDs <- factor(group_IDs, levels = c("Reference", "BCR-XL"))
group_IDs

# indices of all marker columns, lineage markers, and functional markers
# (see Table 1 in Bruggner et al. 2014)
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)

# prepare data into required format
# (note: using lineage markers for clustering, and functional markers for DE testing)
d_se <- prepareData(d_input, sample_IDs, cols_markers, cols_lineage, cols_func)

# check markers (see Table 1, Bruggner et al. 2014)
colnames(d_se)[cols_lineage]
colnames(d_se)[cols_func]


# --------------
# transform data
# --------------

# transform all marker columns
# (using asinh with cofactor = 5; see Bendall et al. 2011, Supp. Fig. S2)
d_se <- transformData(d_se, cofactor = 5)


# ----------
# clustering
# ----------

# generate mini-clusters
seed <- 123
runtime_clustering <- system.time(
  d_se <- generateClusters(d_se, xdim = 20, ydim = 20, seed = seed)
)

# check
nrow(rowData(d_se))                   # number of cells
sum(table(rowData(d_se)$cluster))
length(table(rowData(d_se)$cluster))  # number of clusters


# ---------------------------------------------
# calculate cluster cell counts, medians, ECDFs
# ---------------------------------------------

# calculate cluster cell counts
d_counts <- calcCounts(d_se)
dim(d_counts)

# calculate cluster medians
d_medians <- calcMedians(d_se)
dim(d_medians)

# calculate ECDFs
d_ecdfs <- calcECDFs(d_se)
dim(d_ecdfs)


# --------------------------
# block IDs for paired tests
# --------------------------

# create block IDs for paired tests (one block per patient, since this is a paired data set)
patient_IDs <- factor(gsub("_(BCR-XL|Reference)$", "", sample_IDs))
patient_IDs <- as.numeric(patient_IDs)
patient_IDs

# check sample IDs, group IDs, and block IDs match correctly
sample_IDs
group_IDs
patient_IDs




################################################
# Test for differentially abundant (DA) clusters
################################################

# test for differentially abundant (DA) clusters
runtime_DA <- system.time(
  res_DA <- testDA(d_counts, group_IDs, paired = TRUE, block_IDs = patient_IDs, 
                   plot = TRUE, path = "../../../plots/diffcyt")
)

# show results
rowData(res_DA)

# sort to show top (most highly significant) clusters first
res_DA_sorted <- rowData(res_DA)[order(rowData(res_DA)$adj.P.Val), ]

head(res_DA_sorted, 10)
#View(res_DA_sorted)


# --------
# analysis
# --------

# number of significant DA clusters
table(res_DA_sorted$adj.P.Val < 0.05)


# plots:
# - counts for top DA cluster
# - MST with point color/size for p-value
# - expression profiles of significant clusters, with hierarchical clustering (possibly
# consensus clustering) to group




#############################################################################
# Test for differential expression (DE) of functional markers within clusters
# (method 'diffcyt-med')
#############################################################################

# test for differential expression (DE) of functional markers within clusters
runtime_DE_med <- system.time(
  res_DE_med <- testDE_med(d_counts, d_medians, group_IDs, paired = TRUE, block_IDs = patient_IDs, 
                           plot = TRUE, path = "../../../plots/diffcyt")
)

# show results
rowData(res_DE_med)

# sort to show top (most highly significant) cluster-marker combinations first
res_DE_med_sorted <- rowData(res_DE_med)[order(rowData(res_DE_med)$adj.P.Val), ]

head(res_DE_med_sorted, 10)
#View(res_DE_med_sorted)


# --------
# analysis
# --------

# number of significant DE cluster-marker combinations
table(res_DE_med_sorted$adj.P.Val < 0.05)


# plots:
# - boxplots (with points) of medians for top cluster-marker combination
# - one MST plot for each functional marker, with point color = significance, point size =
# no. of cells
# - heatmap showing up/down for difference in median marker expression (14 functional
# markers) for each cluster




#############################################################################
# Test for differential expression (DE) of functional markers within clusters
# (method 'diffcyt-FDA')
#############################################################################

# ----------
# unweighted
# ----------

# set seed (for permutation tests)
set.seed(123)

# test for differential expression (DE) of functional markers within clusters
runtime_DE_FDA_unwtd <- system.time(
  res_DE_FDA_unwtd <- testDE_FDA(d_counts, d_medians, d_ecdfs, group_IDs, weighted = FALSE, 
                                 paired = TRUE, block_IDs = patient_IDs, 
                                 n_perm = 1000, n_cores = 24)
)

# show results
rowData(res_DE_FDA_unwtd)

# sort to show top (most highly significant) cluster-marker combinations first
res_DE_FDA_unwtd_sorted <- rowData(res_DE_FDA_unwtd)[order(rowData(res_DE_FDA_unwtd)$p_adj), ]

head(res_DE_FDA_unwtd_sorted, 10)
#View(res_DE_FDA_unwtd_sorted)


# --------
# weighted
# --------

# set seed (for permutation tests)
set.seed(123)

# test for differential expression (DE) of functional markers within clusters
runtime_DE_FDA_wtd <- system.time(
  res_DE_FDA_wtd <- testDE_FDA(d_counts, d_medians, d_ecdfs, group_IDs, weighted = TRUE, 
                               paired = TRUE, block_IDs = patient_IDs, 
                               n_perm = 1000, n_cores = 24)
)

# show results
rowData(res_DE_FDA_wtd)

# sort to show top (most highly significant) cluster-marker combinations first
res_DE_FDA_wtd_sorted <- rowData(res_DE_FDA_wtd)[order(rowData(res_DE_FDA_wtd)$p_adj), ]

head(res_DE_FDA_wtd_sorted, 10)
#View(res_DE_FDA_wtd_sorted)


# --------
# analysis
# --------

# number of significant DE cluster-marker combinations

# unweighted
table(res_DE_FDA_unwtd_sorted$adj.P.Val < 0.05)

# weighted
table(res_DE_FDA_wtd_sorted$adj.P.Val < 0.05)


# plots (for both unweighted and weighted):
# - ECDF curves for top cluster-marker combination
# - one MST plot for each functional marker, with point color = significance, point size =
# no. of cells
# - heatmap showing up/down for difference in *median* marker expression (14 functional
# markers) for each cluster

# additional plots:
# - heatmap of correlation matrix between unweighted and weighted
# - runtime for each method (diffcyt-med, diffcyt-FDA-unwtd, diffcyt-FDA-wtd)




#####################
# Session information
#####################

sink("../../../session_info/session_info_BCR_XL_all.txt")
sessionInfo()
sink()



