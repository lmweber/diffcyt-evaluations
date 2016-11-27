############################
# Script to run clustering #
############################

# run clustering to automatically detect cell populations
# using FlowSOM for clustering; lineage markers only

# run after script "load_and_transform_data.R"
#source("load_and_transform_data.R")


library(FlowSOM)

data_flowSet <- as(data_transf, "flowSet")


# set parameters for FlowSOM
xdim <- 20
ydim <- 20
#rlen <- 20

# set seed (required for consistent cluster IDs)
set.seed(100)

# run FlowSOM
system.time({
  fsom <- FlowSOM::ReadInput(data_flowSet, transform = FALSE, scale = FALSE)
  fsom <- FlowSOM::BuildSOM(fsom, colsToUse = markers_lineage, xdim = xdim, ydim = ydim)
  fsom <- FlowSOM::BuildMST(fsom)
})

# extract cluster labels
clus <- fsom$map$mapping[, 1]


# set additional parameters for FlowSOM meta-clustering (if required)
k <- 40

# run FlowSOM meta-clustering (if required; runtime ~1 min)
seed <- 100
system.time({
  fsom_mc <- ConsensusClusterPlus::ConsensusClusterPlus(t(fsom$map$codes), maxK = k, seed = seed)
  fsom_mc <- fsom_mc[[k]]$consensusClass
})

# extract cluster labels
clus_mc <- fsom_mc[clus]



# plot minimum spanning tree (MST) (if required)
# with point sizes scaled to cluster size

pdf("../plots/plot_MST.pdf")
plot(fsom$MST$l, cex = (fsom$MST$size - 7.5) / 3, pch = 19)
dev.off()



