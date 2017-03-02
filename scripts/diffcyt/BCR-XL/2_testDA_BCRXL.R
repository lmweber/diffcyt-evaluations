###############################################################################
# Script to test for differentially abundant (DA) clusters: data set 'BCR-XL' #
###############################################################################


library(limma)


# re-level factor to use "Reference" as base level
group <- factor(group_IDs, levels = c("Reference", "BCR-XL"))


# test for differentially abundant (DA) clusters
res_DA <- testDA(d_clus, group)


# top DA clusters
topTable(res_DA, number = 10)


# plot top differentially abundant (DA) clusters
plotTopDAClusters(res_DA, path = "../../../plots")


