#####################################
# Script to load and transform data #
#####################################

# data from Citrus paper: Bruggner et al. (2014)


library(flowCore)

files <- list.files("../../Citrus_paper_data/experiment_15713_files", full.names = TRUE)

# BCR-XL and Reference samples only
files_BCRXL <- files[grep("patient[1-8]_BCR-XL\\.fcs$", files)]
files_ref <- files[grep("patient[1-8]_Reference\\.fcs$", files)]

fs <- c(files_BCRXL, files_ref)


# read in raw data
data_raw <- lapply(fs, read.FCS, transformation = FALSE, truncate_max_range = FALSE)
names(data_raw) <- basename(fs)


# select lineage markers (for clustering) and functional markers (for differential expression)
# see Table 1 in Bruggner et al. (2014)
all_marker_cols <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
markers_lineage <- c(3, 4, 9, 11, 12, 14, 21, 29, 31, 33)
markers_func <- setdiff(all_marker_cols, markers_lineage)

length(all_marker_cols)
length(markers_lineage)
length(markers_func)


# transform data
# using arcsinh transform with cofactor 5 (see Bendall et al. 2011, Supp. Fig. S2)
cofactor <- 5

data_transf <- lapply(data_raw, function(u) {
  e <- exprs(u)
  e[, all_marker_cols] <- asinh(e[, all_marker_cols] / cofactor)
  exprs(u) <- e
  u
})



