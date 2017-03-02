#######################################################################################
# Script to test for differential expression of functional markers: data set 'BCR-XL' #
# method "diffcyt-med"                                                                #
#######################################################################################


library(limma)


# re-level factor to use "Reference" as base level
group <- factor(group_IDs, levels = c("Reference", "BCR-XL"))


# test for differential expression of functional markers: method "diffcyt-med"
res_DE_med <- testDE_med(d_clus, group, path = "../../../plots")


# top DE clusters across all functional markers
# (note: 'coef = 2' selects contrast of interest, i.e. BCR-XL vs. Reference)
topTable(res_DE_med, coef = 2, number = 20)


