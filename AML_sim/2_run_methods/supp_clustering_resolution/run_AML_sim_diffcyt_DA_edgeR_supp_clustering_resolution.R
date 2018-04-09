##########################################################################################
# Script to run methods
# 
# - method: diffcyt-DA-edgeR
# - data set: AML-sim
# 
# - supplementary results: varying clustering resolution
# 
# Lukas Weber, April 2018
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)


DIR_BENCHMARK <- "../../../../../benchmark_data/AML_sim/data/main"
DIR_RDATA <- "../../../../RData/AML_sim/supp_clustering_resolution"
DIR_SESSION_INFO <- "../../../../session_info/AML_sim/supp_clustering_resolution"




################################
# Loop to run for each threshold
################################

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# condition names
cond_names <- c("CN", "CBF")

# varying clustering resolution: grid size for FlowSOM (e.g. 10x10 grid)
resolution <- c(3, 5, 7, 10, 14, 20, 30, 40)
resolution_sq <- resolution^2

# contrasts (to compare each of 'CN' and 'CBF' vs. 'healthy')
# note: include fixed effects for 'patient_id'
contrasts_list <- list(CN = c(0, 1, 0, 0, 0, 0, 0), CBF = c(0, 0, 1, 0, 0, 0, 0))

# lists to store objects and runtime
out_diffcyt_DA_edgeR_supp_clustering_resolution <- runtime_diffcyt_DA_edgeR_supp_clustering_resolution <- 
  out_clusters_diffcyt_DA_edgeR_supp_clustering_resolution <- out_objects_diffcyt_DA_edgeR_supp_clustering_resolution <- 
  vector("list", length(resolution))
names(out_diffcyt_DA_edgeR_supp_clustering_resolution) <- names(runtime_diffcyt_DA_edgeR_supp_clustering_resolution) <- 
  names(out_clusters_diffcyt_DA_edgeR_supp_clustering_resolution) <- names(out_objects_diffcyt_DA_edgeR_supp_clustering_resolution) <- 
  paste("k", resolution_sq, sep = "_")




for (k in 1:length(resolution)) {
  
  out_diffcyt_DA_edgeR_supp_clustering_resolution[[k]] <- runtime_diffcyt_DA_edgeR_supp_clustering_resolution[[k]] <- 
    out_clusters_diffcyt_DA_edgeR_supp_clustering_resolution[[k]] <- out_objects_diffcyt_DA_edgeR_supp_clustering_resolution[[k]] <- 
    vector("list", length(thresholds))
  names(out_diffcyt_DA_edgeR_supp_clustering_resolution[[k]]) <- names(runtime_diffcyt_DA_edgeR_supp_clustering_resolution[[k]]) <- 
    names(out_clusters_diffcyt_DA_edgeR_supp_clustering_resolution[[k]]) <- names(out_objects_diffcyt_DA_edgeR_supp_clustering_resolution[[k]]) <- 
    thresholds
  
  
  for (th in 1:length(thresholds)) {
    
    ###########################
    # Load data, pre-processing
    ###########################
    
    # filenames
    
    files_healthy <- list.files(file.path(DIR_BENCHMARK, "healthy"), 
                                pattern = "\\.fcs$", full.names = TRUE)
    files_CN <- list.files(file.path(DIR_BENCHMARK, "CN"), 
                           pattern = paste0("_", thresholds[th], "\\.fcs$"), full.names = TRUE)
    files_CBF <- list.files(file.path(DIR_BENCHMARK, "CBF"), 
                            pattern = paste0("_", thresholds[th], "\\.fcs$"), full.names = TRUE)
    
    files_load <- c(files_healthy, files_CN, files_CBF)
    files_load
    
    # load data
    
    d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)
    
    # sample IDs, group IDs, patient IDs
    sample_id <- gsub("(_[0-9]+pc$)|(_0\\.[0-9]+pc$)", "", 
                      gsub("^AML_sim_", "", 
                           gsub("\\.fcs$", "", basename(files_load))))
    sample_id
    
    group_id <- factor(gsub("_.*$", "", sample_id), levels = c("healthy", "CN", "CBF"))
    group_id
    
    patient_id <- factor(gsub("^.*_", "", sample_id))
    patient_id
    
    experiment_info <- data.frame(group_id, patient_id, sample_id)
    experiment_info
    
    # marker information
    
    # indices of all marker columns, lineage markers, and functional markers
    # (16 surface markers / 15 functional markers; see Levine et al. 2015, Supplemental 
    # Information, p. 4)
    cols_markers <- 11:41
    cols_lineage <- c(35, 29, 14, 30, 12, 26, 17, 33, 41, 32, 22, 40, 27, 37, 23, 39)
    cols_func <- setdiff(cols_markers, cols_lineage)
    
    stopifnot(all(sapply(seq_along(d_input), function(i) all(colnames(d_input[[i]]) == colnames(d_input[[1]])))))
    
    marker_name <- colnames(d_input[[1]])
    marker_name <- gsub("\\(.*$", "", marker_name)
    
    marker_class <- rep("none", length(marker_name))
    marker_class[cols_lineage] <- "cell_type"
    marker_class[cols_func] <- "cell_state"
    marker_class <- factor(marker_class, levels = c("cell_type", "cell_state", "none"))
    
    marker_info <- data.frame(marker_name, marker_class)
    marker_info
    
    
    
    
    ##################
    # diffcyt pipeline
    ##################
    
    # --------------------
    # pre-processing steps
    # --------------------
    
    runtime_preprocessing <- system.time({
      
      # prepare data into required format
      d_se <- prepareData(d_input, experiment_info, marker_info)
      
      colnames(d_se)[colData(d_se)$marker_class == "cell_type"]
      colnames(d_se)[colData(d_se)$marker_class == "cell_state"]
      
      # transform data
      d_se <- transformData(d_se, cofactor = 5)
      
      # clustering
      # (runtime: ~30 sec with xdim = 20, ydim = 20)
      # note: varying clustering resolution
      seed <- 123
      d_se <- generateClusters(d_se, xdim = resolution[k], ydim = resolution[k], seed = seed)
      
      length(table(rowData(d_se)$cluster_id))  # number of clusters
      nrow(rowData(d_se))                      # number of cells
      sum(table(rowData(d_se)$cluster_id))
      min(table(rowData(d_se)$cluster_id))     # size of smallest cluster
      max(table(rowData(d_se)$cluster_id))     # size of largest cluster
      
      # calculate cluster cell counts
      d_counts <- calcCounts(d_se)
      
      dim(d_counts)
      rowData(d_counts)
      length(assays(d_counts))
      
      # calculate cluster medians
      d_medians <- calcMedians(d_se)
      
      dim(d_medians)
      rowData(d_medians)
      length(assays(d_medians))
      names(assays(d_medians))
      
      # calculate medians by cluster and marker
      d_medians_by_cluster_marker <- calcMediansByClusterMarker(d_se)
      
      dim(d_medians_by_cluster_marker)
      length(assays(d_medians_by_cluster_marker))
      
      # calculate medians by sample and marker
      d_medians_by_sample_marker <- calcMediansBySampleMarker(d_se)
      
      dim(d_medians_by_sample_marker)
      length(assays(d_medians_by_sample_marker))
      
    })
    
    
    # ---------------------------------
    # store data objects (for plotting)
    # ---------------------------------
    
    out_objects_diffcyt_DA_edgeR_supp_clustering_resolution[[k]][[th]] <- list(
      d_se = d_se, 
      d_counts = d_counts, 
      d_medians = d_medians, 
      d_medians_by_cluster_marker = d_medians_by_cluster_marker, 
      d_medians_by_sample_marker = d_medians_by_sample_marker
    )
    
    
    # -----------------------------------------
    # test for differentially abundant clusters
    # -----------------------------------------
    
    # note: test separately for each condition: CN vs. healthy, CBF vs. healthy
    
    out_diffcyt_DA_edgeR_supp_clustering_resolution[[k]][[th]] <- runtime_diffcyt_DA_edgeR_supp_clustering_resolution[[k]][[th]] <- 
      out_clusters_diffcyt_DA_edgeR_supp_clustering_resolution[[k]][[th]] <- vector("list", length(cond_names))
    names(out_diffcyt_DA_edgeR_supp_clustering_resolution[[k]][[th]]) <- names(runtime_diffcyt_DA_edgeR_supp_clustering_resolution[[k]][[th]]) <- 
      names(out_clusters_diffcyt_DA_edgeR_supp_clustering_resolution[[k]][[th]]) <- cond_names
    
    
    for (j in 1:length(cond_names)) {
      
      runtime_j <- system.time({
        
        # set up design matrix
        # note: include fixed effects for 'patient_id'
        design <- createDesignMatrix(experiment_info, cols_design = 1:2)
        design
        
        # set up contrast matrix
        contrast <- createContrast(contrasts_list[[j]])
        contrast
        
        # run tests
        # note: adjust filtering parameter 'min_samples' (since there are 3 conditions)
        res <- testDA_edgeR(d_counts, design, contrast, 
                            min_cells = 3, min_samples = nrow(experiment_info) / 3)
        
      })
      
      # show results
      rowData(res)
      
      # sort to show top (most highly significant) clusters first
      res_sorted <- rowData(res)[order(rowData(res)$FDR), ]
      print(head(res_sorted, 10))
      #View(as.data.frame(res_sorted))
      
      # number of significant tests (note: one test per cluster)
      print(table(res_sorted$FDR <= 0.1))
      
      # runtime
      runtime_total <- runtime_preprocessing[["elapsed"]] + runtime_j[["elapsed"]]
      print(runtime_total)
      
      runtime_diffcyt_DA_edgeR_supp_clustering_resolution[[k]][[th]][[j]] <- runtime_total
      
      
      # ---------------------------------------------
      # store results at cluster level (for plotting)
      # ---------------------------------------------
      
      res_clusters <- as.data.frame(rowData(res))
      
      out_clusters_diffcyt_DA_edgeR_supp_clustering_resolution[[k]][[th]][[j]] <- res_clusters
      
      
      
      
      ##############################
      # Return results at cell level
      ##############################
      
      # Note: diffcyt methods return results for each cluster. To enable performance
      # comparisons between methods at the cell level, we assign the cluster-level p-values
      # to all cells within each cluster.
      
      
      # number of cells per sample (including spike-in cells)
      n_cells <- sapply(d_input, nrow)
      
      # spike-in status for each cell
      is_spikein <- unlist(sapply(d_input, function(d) exprs(d)[, "spikein"]))
      stopifnot(length(is_spikein) == sum(n_cells))
      
      # select samples for this condition and healthy
      ix_keep_cnd <- group_id %in% c("healthy", cond_names[j])
      
      
      # match cluster-level p-values to individual cells
      
      stopifnot(nrow(rowData(res)) == length(levels(rowData(d_se)$cluster_id)), 
                all(rowData(res)$cluster_id == levels(rowData(d_se)$cluster_id)))
      
      rowData(res)$cluster_id <- factor(rowData(res)$cluster_id, levels = levels(rowData(d_se)$cluster_id))
      
      # match cells to clusters
      ix_match <- match(rowData(d_se)$cluster_id, rowData(res)$cluster_id)
      
      p_vals_clusters <- rowData(res)$PValue
      p_adj_clusters <- rowData(res)$FDR
      
      p_vals_cells <- p_vals_clusters[ix_match]
      p_adj_cells <- p_adj_clusters[ix_match]
      
      
      # set up data frame with results and true spike-in status at cell level
      
      which_cnd <- rep(ix_keep_cnd, n_cells)
      is_spikein_cnd <- is_spikein[which_cnd]
      
      stopifnot(length(p_vals_cells[which_cnd]) == length(is_spikein_cnd), 
                length(p_adj_cells[which_cnd]) == length(is_spikein_cnd))
      
      res_p_vals <- p_vals_cells[which_cnd]
      res_p_adj <- p_adj_cells[which_cnd]
      
      # replace NAs (due to filtering) to ensure same cells are returned for all methods
      res_p_vals[is.na(res_p_vals)] <- 1
      res_p_adj[is.na(res_p_adj)] <- 1
      
      # return values for this condition and healthy
      res <- data.frame(p_vals = res_p_vals, 
                        p_adj = res_p_adj, 
                        spikein = is_spikein_cnd)
      
      # store results
      out_diffcyt_DA_edgeR_supp_clustering_resolution[[k]][[th]][[j]] <- res
      
    }
  }
}




#####################
# Save output objects
#####################

save(out_diffcyt_DA_edgeR_supp_clustering_resolution, runtime_diffcyt_DA_edgeR_supp_clustering_resolution, 
     file = file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_edgeR_supp_clustering_resolution.RData"))

save(out_clusters_diffcyt_DA_edgeR_supp_clustering_resolution, 
     file = file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_edgeR_supp_clustering_resolution.RData"))

save(out_objects_diffcyt_DA_edgeR_supp_clustering_resolution, 
     file = file.path(DIR_RDATA, "out_objects_AML_sim_diffcyt_DA_edgeR_supp_clustering_resolution.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_AML_sim_diffcyt_DA_edgeR_supp_clustering_resolution.txt"))
sessionInfo()
sink()



