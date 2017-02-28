###########################################################################
# Script to run pipeline of 'diffcyt' analysis scripts: data set 'BCR-XL' #
###########################################################################


# all steps prior to statistical tests: load data, transform, clustering, medians and frequencies
source("1_preprocessing_BCRXL.R")


# test for differentially abundant (DA) clusters
source("2a_testDA_BCRXL.R")


# test for differential expression of functional markers: method "diffcyt-med"
source("2b_testDE_med_BCRXL.R")


# next: add scripts for methods "diffcyt-FDA" and "diffcyt-KS"


