######################################################################
# Script to run pre-processing: all steps prior to statistical tests #
# for data set 'BCR-XL'                                              #
######################################################################


library(diffcyt)
library(flowCore)
library(magrittr)



#############
# LOAD DATA #
#############

# filenames
files <- list.files("../../../../benchmark_data/Citrus_paper_data/experiment_15713_files", 
                    pattern = "\\.fcs$", full.names = TRUE)

files_BCRXL <- files[grep("patient[1-8]_BCR-XL\\.fcs$", files)]
files_ref <- files[grep("patient[1-8]_Reference\\.fcs$", files)]

# load data
files_load <- c(files_BCRXL, files_ref)
d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)

# sample IDs and group IDs
sample_IDs <- gsub("^PBMC8_30min_", "", gsub("\\.fcs$", "", basename(files_load)))
sample_IDs
group_IDs <- gsub("^patient[0-9]_", "", sample_IDs)
group_IDs

# indices of all marker columns, lineage markers, and functional markers (see Table 1 in
# Bruggner et al. 2014)
marker_cols <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
lineage_cols <- c(3:4, 9, 11,12,14, 21, 29, 31, 33)
functional_cols <- setdiff(marker_cols, lineage_cols)

# prepare data into required format
d_se <- prepareData(d_input, sample_IDs, group_IDs, marker_cols, lineage_cols, functional_cols)

# check markers (see Table 1, Bruggner et al. 2014)
colnames(d_se)[lineage_cols]
colnames(d_se)[functional_cols]



#############
# TRANSFORM #
#############

# transform all marker columns (using asinh with cofactor = 5; see Bendall et al. 2011, 
# Supp. Fig. S2)
d_se <- transformData(d_se, cofactor = 5)



##############
# CLUSTERING #
##############

# generate clusters
d_se <- generateClusters(d_se, cols_to_use = lineage_cols, xdim = 20, ydim = 20, 
                         seed = 123, plot = TRUE, path = "../../../plots")

# check
length(SummarizedExperiment::rowData(d_se)$cluster)         # number of cells
length(table(SummarizedExperiment::rowData(d_se)$cluster))  # number of clusters



###################################
# CLUSTER MEDIANS AND FREQUENCIES #
###################################

# calculate cluster medians and frequencies
d_clus <- calcMediansAndFreq(d_se)



###################
# CALCULATE ECDFS #
###################

# calculate ECDFs
d_ecdfs <- calcECDFs(d_se, resolution = 50)



