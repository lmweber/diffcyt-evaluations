##########################################################################################
# Script to prepare benchmark data set 'BCR-XL-sim'
# 
# The original 'BCR-XL' data set is sourced from Bodenmiller et al. (2012), and was
# previously used for benchmark evaluations by Bruggner et al. (2014) (Citrus paper).
# 
# Raw data downloaded from Cytobank (experiment 15713)
# - see Citrus wiki (section "PBMC Example 1"):
# https://github.com/nolanlab/citrus/wiki/PBMC-Example-1
# - direct link to Cytobank repository:
# https://community.cytobank.org/cytobank/experiments/15713/download_files
# 
# Cell population labels are reproduced from Nowicka et al. (2017), F1000Research, using a
# strategy of expert-guided manual merging of automatically generated clusters from the
# FlowSOM algorithm. Code to reproduce the cell population labels is available in the
# script "cell_population_labels_BCR_XL.R".
# 
# The simulations in this script are generated as follows:
# - select reference (unstimulated) samples from the main 'BCR-XL' data set
# - randomly split each sample into two halves
# - in one half, replace B cells with B cells from the corresponding stimulated sample
# - adjust "difficulty" of the simulation by scaling the average difference in pS6 signal
# 
# Methods can then be evaluated by their ability to detect the known strong differential
# signal in pS6 expression.
# 
# Lukas Weber, August 2017
##########################################################################################


library(flowCore)




###########
# LOAD DATA
###########

# load .fcs files
DIR_RAW_DATA <- "../../../benchmark_data/BCR_XL/raw_data/experiment_15713_files"

files <- list.files(DIR_RAW_DATA, pattern = "\\.fcs$", full.names = TRUE)

files_BCRXL <- files[grep("patient[1-8]_BCR-XL\\.fcs$", files)]
files_ref <- files[grep("patient[1-8]_Reference\\.fcs$", files)]

files_all <- c(files_BCRXL, files_ref)

data <- lapply(files_all, function(f) exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)))


# load population labels
DIR_LABELS <- "../../../benchmark_data/BCR_XL/population_IDs"
file_labels <- list.files(DIR_LABELS, pattern = "\\.csv$", full.names = TRUE)
data_labels_raw <- read.csv(file_labels)

# match order of labels to order of .fcs files
files_all
unique(data_labels_raw$sample)

d_labels <- split(data_labels_raw, data_labels_raw$sample)

# note 'split' has automatically sorted by factor levels in alphabetical order
d_labels <- do.call("rbind", d_labels)

# check
files_all
unique(d_labels$sample)
#View(d_labels)




#############################
# EXPORT FILES: NO SIMULATION
#############################

DIR_DATA_NOSIM <- "../../../benchmark_data/BCR_XL/data/nosim"

cmds <- paste("cp", files_all, DIR_DATA_NOSIM)

for (cmd in cmds) {
  system(cmd)
}




##################
# SIMULATION: FULL
##################

# 'full' simulation contains the full differential expression signal for pS6


# load data for reference samples
data_ref <- lapply(files_ref, function(f) exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)))



