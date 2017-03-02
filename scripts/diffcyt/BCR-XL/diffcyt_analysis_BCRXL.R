###########################################################################
# Script to run pipeline of 'diffcyt' analysis scripts: data set 'BCR-XL' #
###########################################################################


library(diffcyt)


# all steps prior to statistical tests: load data, transform data, clustering, calculate
# medians and frequencies, calculate ECDFs
source("1_preprocessing_BCRXL.R")


# test for differentially abundant (DA) clusters
source("2_testDA_BCRXL.R")


# test for differential expression of functional markers: method "diffcyt-med"
source("2a_testDE_med_BCRXL.R")


# test for differential expression of functional markers: method "diffcyt-FDA"
# (runtime: 30 mins)
source("2b_testDE_FDA_BCRXL.R")


