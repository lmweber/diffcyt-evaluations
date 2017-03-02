#######################################################################################
# Script to test for differential expression of functional markers: data set 'BCR-XL' #
# method "diffcyt-FDA"                                                                #
#######################################################################################


library(fda)


# re-level factor to use "Reference" as base level
group <- factor(group_IDs, levels = c("Reference", "BCR-XL"))


# test for differential expression of functional markers: method "diffcyt-FDA"
# (runtime: 30 mins)
res_DE_FDA <- testDE_FDA(d_ecdfs, group, n_perm = 5000)


# results for top DE clusters (across all functional markers)
head(res_DE_FDA, 20)


# plot ECDFs for top DE cluster
plotECDFsTopDECluster(d_ecdfs, path = "../../../plots")


