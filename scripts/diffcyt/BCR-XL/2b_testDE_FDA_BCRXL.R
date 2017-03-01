#######################################################################################
# Script to test for differential expression of functional markers: data set 'BCR-XL' #
# method "diffcyt-FDA"                                                                #
#######################################################################################


library(fda)


# calculate ECDFs
d_ecdfs <- calcECDFs(d_se, resolution = 30)


# test for differential expression of functional markers: method "diffcyt-FDA"
res_DE_FDA <- testDE_FDA(d_clus, group, path = "../../../plots")


# plot ECDFs for top DE cluster
plotECDFsTopDECluster(d_ecdfs, path = "../../../plots")


# top DE clusters across all functional markers
# (note: 'coef = 2' selects contrast of interest, i.e. BCR-XL vs. Reference)
topTable(res_DE_FDA, coef = 2, number = 20)


