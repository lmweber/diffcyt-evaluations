##########################################################################################
# Generate summary tables
# 
# - data set: AML-sim
# - tables: summary of filtering
# - method: diffcyt methods
# 
# - main results
# 
# Lukas Weber, May 2018
##########################################################################################


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/main"

load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_voom_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_GLMM_main.RData"))




#########################
# Generate summary tables
#########################

# loop over methods (m), thresholds (th), conditions (j)

# methods
method_names <- c("diffcyt_DA_edgeR", "diffcyt_DA_voom", "diffcyt_DA_GLMM")

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# condition names
cond_names <- c("CN", "CBF")


# total number of clusters for this data set
n_clus_total <- 400


# combined list of objects
out_clusters <- list(
  out_clusters_diffcyt_DA_edgeR_main, 
  out_clusters_diffcyt_DA_voom_main, 
  out_clusters_diffcyt_DA_GLMM_main
)
names(out_clusters) <- method_names


# store in lists
summary_kept <- summary_filtered <- vector("list", length(method_names))
names(summary_kept) <- names(summary_filtered) <- method_names


for (m in 1:length(method_names)) {
  
  # store numbers of clusters in matrix
  n_clus_kept <- matrix(as.numeric(NA), nrow = length(cond_names), ncol = length(thresholds))
  rownames(n_clus_kept) <- cond_names
  colnames(n_clus_kept) <- thresholds
  
  
  for (th in 1:length(thresholds)) {
    
    for (j in 1:length(cond_names)) {
      
      # numbers of clusters kept
      n_clus <- nrow(out_clusters[[m]][[th]][[j]])
      n_clus_kept[j, th] <- n_clus
      
    }
  }
  
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



