###############################################
# Script to run clustering: data set 'BCR-XL' #
###############################################


# run after script "transform_BCRXL.R"


suppressPackageStartupMessages(library(FlowSOM))


clus_BCRXL <- generateClusters(data_BCRXL, 
                               cols_to_use = lineage_cols, 
                               xdim = 20, ydim = 20, seed = 123, 
                               plot = FALSE)

# check cluster labels
table(clus_BCRXL)              # cluster abundances
length(table(clus_BCRXL))      # total 400 clusters
sum(is.na(table(clus_BCRXL)))
sum(table(clus_BCRXL))         # total number of cells


