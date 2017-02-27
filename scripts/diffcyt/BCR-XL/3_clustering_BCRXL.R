###############################################
# Script to run clustering: data set 'BCR-XL' #
###############################################


suppressPackageStartupMessages(library(FlowSOM))


clus_BCRXL <- generateClusters(data_BCRXL, 
                               cols_to_use = lineage_cols, 
                               xdim = 20, ydim = 20, seed = 123, 
                               plot = FALSE)


# check

length(clus_BCRXL)               # number of samples
names(clus_BCRXL)                # sample names

sapply(clus_BCRXL, length)       # number of cells per sample
sum(sapply(clus_BCRXL, length))  # total number of cells

clus_BCRXL_all <- do.call("c", clus_BCRXL)
length(clus_BCRXL_all)

table(clus_BCRXL_all)          # cluster abundances
length(table(clus_BCRXL_all))  # total 400 clusters
sum(is.na(table(clus_BCRXL_all)))
sum(table(clus_BCRXL_all))


