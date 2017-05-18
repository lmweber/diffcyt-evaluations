###################################################################
# Script to run 'diffcyt' analysis pipeline for data set 'BCR-XL' #
# Lukas Weber, March 2017                                         #
###################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)
library(limma)



####################################
# (1) Load data and pre-processing #
####################################

# ---------
# load data
# ---------

# filenames
files <- list.files("../../../../benchmark_data/Citrus_paper_data/experiment_15713_files", 
                    pattern = "\\.fcs$", full.names = TRUE)

files_BCRXL <- files[grep("patient[1-8]_BCR-XL\\.fcs$", files)]
files_ref <- files[grep("patient[1-8]_Reference\\.fcs$", files)]

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

# indices of all marker columns, lineage markers, and functional markers (see Table 1 in
# Bruggner et al. 2014)
marker_cols <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
lineage_cols <- c(3:4, 9, 11,12,14, 21, 29, 31, 33)
functional_cols <- setdiff(marker_cols, lineage_cols)

# prepare data into required format
# (note: using lineage markers for clustering, and functional markers for DE testing)
d_se <- prepareData(d_input, sample_IDs, marker_cols, lineage_cols, functional_cols)

# check markers (see Table 1, Bruggner et al. 2014)
colnames(d_se)[lineage_cols]
colnames(d_se)[functional_cols]


# ---------
# transform
# ---------

# transform all marker columns (using asinh with cofactor = 5; see Bendall et al. 2011, 
# Supp. Fig. S2)
d_se <- transformData(d_se, cofactor = 5)


# ----------
# clustering
# ----------

# generate micro-clusters
d_se <- generateClusters(d_se, xdim = 20, ydim = 20, seed = 123)

# check
length(rowData(d_se)$cluster)         # number of cells
length(table(rowData(d_se)$cluster))  # number of clusters


# --------------------------------------------
# calculate cluster counts, medians, and ECDFs
# --------------------------------------------

# calculate cluster cell counts
d_counts <- calcCounts(d_se)
dim(d_counts)

# calculate cluster medians (median marker expression)
d_medians <- calcMedians(d_se)
dim(d_medians)

# calculate ECDFs
d_ecdfs <- calcECDFs(d_se)
dim(d_ecdfs)



######################################################
# (2) Test for differentially abundant (DA) clusters #
######################################################

# check group IDs factor has correct level as reference
group_IDs

# check order of group IDs factor matches with order of samples (columns) in data objects
colnames(d_counts)
colnames(d_medians)
colnames(d_ecdfs)
group_IDs

# create patient IDs (block IDs) for paired tests (for this data set: 1 block per patient)
patient_IDs <- factor(gsub("_(BCR-XL|Reference)$", "", sample_IDs))
patient_IDs <- as.numeric(patient_IDs)
patient_IDs

# test for differentially abundant (DA) clusters
res_DA <- testDA(d_counts, group_IDs, 
                 paired = TRUE, block_IDs = patient_IDs, 
                 plot = TRUE, path = "../../../plots")

# show results for top DA clusters
topTable(res_DA, number = 10)

### results: most micro-clusters are DA (but not all). Try to aggregate them using 
### hierarchical consensus clustering (then can classify whole groups as DA).

# [TO DO: plot counts for top DA clusters (similar to previous version)]

# [TO DO: plot MST with point color/size for p-value significance]

# [TO DO: plot expression profiles of significant clusters]

# [TO DO: implement cluster/hypothesis merging and aggregation: consensus clustering
# trees, aggregate p-values]



##########################################################################################
# (3b) Test for differential expression of markers within clusters: method "diffcyt-FDA" #
##########################################################################################

# runtime: ........... with 2 cores on Mac ................. 1 hr 20 mins with a single core
# runtime: ........ (unpaired tests)
# (paired tests)
system.time(
  res_DE_FDA <- testDE_FDA(d_ecdfs, group_IDs, weighted = TRUE, d_counts = d_counts)
)

# results for top DE clusters (across all functional markers)
head(res_DE_FDA, 20)



# plot ECDFs for top DE cluster
#plotECDFsTopDECluster(d_ecdfs, path = "../../../plots")


# [next: add a heatmap of all the top clusters, to identify them (B cells etc); do the
# results recover the same biologically meaningful results from the Citrus paper?]





#####################################################################################
# (3a) TEST FOR DIFFERENTIAL EXPRESSION OF FUNCTIONAL MARKERS: method "diffcyt-med" #
#####################################################################################

# re-level factor to use "Reference" as base level
#group <- factor(group_IDs, levels = c("Reference", "BCR-XL"))

# test for differential expression of functional markers: method "diffcyt-med"
res_DE_med <- testDE_med(d_clus, group, path = "../../../plots")

# top DE clusters across all functional markers
# (note: 'coef = 2' selects contrast of interest, i.e. BCR-XL vs. Reference)
topTable(res_DE_med, coef = 2, number = 20)








