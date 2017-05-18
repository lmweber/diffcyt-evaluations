##########################################################################################
# Script to prepare benchmark data set 'BCR-XL'
#
# 'BCR-XL' data set is sourced from Bodenmiller et al. (2012), and was previously used for
# benchmark evaluations by Bruggner et al. (2014) (Citrus paper)
#
# Raw data downloaded from Cytobank (experiment 15713)
# - see Citrus wiki (section "PBMC Example 1"):
# https://github.com/nolanlab/citrus/wiki/PBMC-Example-1
# - direct link to Cytobank repository:
# https://community.cytobank.org/cytobank/experiments/15713/download_files
#
# Lukas Weber, May 2017
##########################################################################################


# Select files from BCR-XL and Reference conditions, and copy to 'data' folder. No 
# additional pre-processing required.

DIR_RAW_DATA <- "../../../../benchmark_data/BCR_XL/raw_data/experiment_15713_files"
DIR_DATA <- "../../../../benchmark_data/BCR_XL/data"

files <- list.files(DIR_RAW_DATA, pattern = "\\.fcs$", full.names = TRUE)

files_BCRXL <- files[grep("patient[1-8]_BCR-XL\\.fcs$", files)]
files_ref <- files[grep("patient[1-8]_Reference\\.fcs$", files)]

cmds <- paste("cp", c(files_BCRXL, files_ref), DIR_DATA)

for (cmd in cmds) {
  system(cmd)
}


