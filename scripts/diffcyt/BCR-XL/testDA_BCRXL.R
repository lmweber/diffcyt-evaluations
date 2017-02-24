#######################################################################################
# Script to calculate tests for differential abundance of clusters: data set 'BCR-XL' #
#######################################################################################


# run after script "clustering_BCRXL.R"


library(limma)


# group membership for each sample (use 'ref' as reference level in factor)
smp <- rownames(flowCore::phenoData(data_BCRXL))
group <- gsub("^patient[0-9]_", "", smp)
group <- factor(group, levels = rev(unique(group)))  # flips reference level

# check
smp
group
as.numeric(group)


# test for differentially abundant (DA) clusters
f1_BCRXL <- testDA(data_BCRXL, clus_BCRXL, group)


# display results for top DA clusters (top 6 clusters)
topTable(f1_BCRXL, coef = 2, number = 6)


# plot for top DA cluster (one cluster)
plotTopDACluster(f1_BCRXL, filename = "../../../plots/plot_top_DA_cluster.pdf")


