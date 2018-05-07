##########################################################################################
# Generate summary tables
# 
# - data set: BCR-XL-sim
# - tables: summary of filtering
# - method: diffcyt methods
# 
# - main results
# 
# Lukas Weber, May 2018
##########################################################################################


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL_sim/main"

load(file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DS_limma_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DS_LMM_main.RData"))




#########################
# Generate summary tables
#########################

# loop over methods (m)

# methods
method_names <- c("diffcyt_DS_limma", "diffcyt_DS_LMM")


# total number of clusters for this data set
n_clus_total <- 100
# total number of cell state markers
n_state_markers <- 14


# combined list of objects
out_clusters <- list(
  out_clusters_diffcyt_DS_limma_main, 
  out_clusters_diffcyt_DS_LMM_main
)
names(out_clusters) <- method_names


# store in lists
summary_kept <- summary_filtered <- vector("list", length(method_names))
names(summary_kept) <- names(summary_filtered) <- method_names


for (m in 1:length(method_names)) {
  
  # numbers of clusters kept
  n_clus_kept <- nrow(out_clusters[[m]]) / n_state_markers
  
  # numbers of clusters filtered
  n_clus_filtered <- n_clus_kept - n_clus_total
  
  # store in lists
  summary_kept[[m]] <- n_clus_kept
  summary_filtered[[m]] <- n_clus_filtered
  
}




#####################
# Show summary tables
#####################

# numbers of clusters kept
print(summary_kept)

# numbers of clusters filtered
print(summary_filtered)



