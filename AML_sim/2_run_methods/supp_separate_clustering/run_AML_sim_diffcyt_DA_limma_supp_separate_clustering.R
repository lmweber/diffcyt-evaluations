##########################################################################################
# Script to run methods
# 
# - method: diffcyt-DA-limma
# - data set: AML-sim
# 
# - supplementary results: separate clustering for each condition
# 
# Lukas Weber, April 2018
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)


DIR_BENCHMARK <- "../../../../../benchmark_data/AML_sim/data/main"
DIR_PLOTS <- "../../../../plots/AML_sim/supp_separate_clustering/diagnostic/diffcyt_DA_limma"
DIR_RDATA <- "../../../../RData/AML_sim/supp_separate_clustering"
DIR_SESSION_INFO <- "../../../../session_info/AML_sim/supp_separate_clustering"




################################
# Loop to run for each threshold
################################

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# condition names
cond_names <- c("CN", "CBF")

# contrasts (to compare each of 'CN' and 'CBF' vs. 'healthy')
# note: include fixed effects for 'patient'
# note: running whole pipeline separately for each condition, so CN or CBF term is in same column
contrast_vec <- c(0, 1, 0, 0, 0, 0)

# lists to store objects and runtime
out_diffcyt_DA_limma_supp_separate_clustering <- runtime_diffcyt_DA_limma_supp_separate_clustering <- 
  out_clusters_diffcyt_DA_limma_supp_separate_clustering <- out_objects_diffcyt_DA_limma_supp_separate_clustering <- 
  vector("list", length(thresholds))
names(out_diffcyt_DA_limma_supp_separate_clustering) <- names(runtime_diffcyt_DA_limma_supp_separate_clustering) <- 
  names(out_clusters_diffcyt_DA_limma_supp_separate_clustering) <- names(out_objects_diffcyt_DA_limma_supp_separate_clustering) <- 
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
  sample_IDs <- gsub("(_[0-9]+pc$)|(_0\\.[0-9]+pc$)", "", 
                     gsub("^AML_sim_", "", 
                          gsub("\\.fcs$", "", basename(files_load))))
  sample_IDs
  
  group_IDs <- factor(gsub("_.*$", "", sample_IDs), levels = c("healthy", "CN", "CBF"))
  group_IDs
  
  patient_IDs <- factor(gsub("^.*_", "", sample_IDs))
  patient_IDs
  
  sample_info <- data.frame(group = group_IDs, patient = patient_IDs, sample = sample_IDs)
  sample_info
  
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
  
  is_marker <- rep(FALSE, length(marker_name))
  is_marker[cols_markers] <- TRUE
  
  marker_type <- rep("none", length(marker_name))
  marker_type[cols_lineage] <- "cell_type"
  marker_type[cols_func] <- "cell_state"
  marker_type <- factor(marker_type, levels = c("cell_type", "cell_state", "none"))
  
  marker_info <- data.frame(marker_name, is_marker, marker_type)
  marker_info
  
  
  
  
  ##################
  # diffcyt pipeline
  ##################
  
  # run clustering separately for each condition
  
  
  out_diffcyt_DA_limma_supp_separate_clustering[[th]] <- runtime_diffcyt_DA_limma_supp_separate_clustering[[th]] <- 
    out_clusters_diffcyt_DA_limma_supp_separate_clustering[[th]] <- out_objects_diffcyt_DA_limma_supp_separate_clustering[[th]] <- 
    vector("list", length(cond_names))
  names(out_diffcyt_DA_limma_supp_separate_clustering[[th]]) <- names(runtime_diffcyt_DA_limma_supp_separate_clustering[[th]]) <- 
    names(out_clusters_diffcyt_DA_limma_supp_separate_clustering[[th]]) <- names(out_objects_diffcyt_DA_limma_supp_separate_clustering[[th]]) <- 
    cond_names
  
  
  for (j in 1:length(cond_names)) {
    
    # ------------------------------
    # subset data for each condition
    # ------------------------------
    
    ix_subset <- group_IDs %in% c("healthy", cond_names[j])
    
    d_input_sub <- d_input[ix_subset]
    sample_info_sub <- droplevels(sample_info[ix_subset, ])
    
    
    # --------------------
    # pre-processing steps
    # --------------------
    
    runtime_preprocessing <- system.time({
      
      # prepare data into required format
      d_se <- prepareData(d_input_sub, sample_info_sub, marker_info)
      
      colnames(d_se)[colData(d_se)$marker_type == "cell_type"]
      colnames(d_se)[colData(d_se)$marker_type == "cell_state"]
      
      # transform data
      d_se <- transformData(d_se, cofactor = 5)
      
      # clustering
      # (runtime: ~30 sec with xdim = 20, ydim = 20)
      seed <- 123
      d_se <- generateClusters(d_se, xdim = 20, ydim = 20, seed = seed)
      
      length(table(rowData(d_se)$cluster))  # number of clusters
      nrow(rowData(d_se))                   # number of cells
      sum(table(rowData(d_se)$cluster))
      min(table(rowData(d_se)$cluster))     # size of smallest cluster
      max(table(rowData(d_se)$cluster))     # size of largest cluster
      
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
    
    out_objects_diffcyt_DA_limma_supp_separate_clustering[[th]][[j]] <- list(
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
    
    
    runtime_j <- system.time({
      
      # set up design matrix
      # note: include fixed effects for 'patient'
      design <- createDesignMatrix(sample_info, cols_include = 1:2)
      design
      
      # set up contrast matrix
      contrast <- createContrast(contrast_vec)
      contrast
      
      # run tests
      # note: adjust filtering parameter 'min_samples' (since there are 2 conditions)
      path <- file.path(DIR_PLOTS, thresholds[th], cond_names[j])
      res <- testDA_limma(d_counts, design, contrast, 
                          min_cells = 3, min_samples = nrow(sample_info) / 2, 
                          path = path)
      
    })
    
    # show results
    rowData(res)
    
    # sort to show top (most highly significant) clusters first
    res_sorted <- rowData(res)[order(rowData(res)$adj.P.Val), ]
    print(head(res_sorted, 10))
    #View(as.data.frame(res_sorted))
    
    # number of significant tests (note: one test per cluster)
    print(table(res_sorted$adj.P.Val <= 0.1))
    
    # runtime
    runtime_total <- runtime_preprocessing[["elapsed"]] + runtime_j[["elapsed"]]
    print(runtime_total)
    
    runtime_diffcyt_DA_limma_supp_separate_clustering[[th]][[j]] <- runtime_total
    
    
    # ---------------------------------------------
    # store results at cluster level (for plotting)
    # ---------------------------------------------
    
    res_clusters <- as.data.frame(rowData(res))
    
    out_clusters_diffcyt_DA_limma_supp_separate_clustering[[th]][[j]] <- res_clusters
    
    
    
    
    ##############################
    # Return results at cell level
    ##############################
    
    # Note: diffcyt methods return results for each cluster. To enable performance
    # comparisons between methods at the cell level, we assign the cluster-level p-values
    # to all cells within each cluster.
    
    
    # number of cells per sample (including spike-in cells)
    n_cells <- sapply(d_input_sub, nrow)
    
    # spike-in status for each cell
    is_spikein <- unlist(sapply(d_input_sub, function(d) exprs(d)[, "spikein"]))
    stopifnot(length(is_spikein) == sum(n_cells))
    
    
    # match cluster-level p-values to individual cells
    
    stopifnot(nrow(rowData(res)) == length(levels(rowData(d_se)$cluster)), 
              all(rowData(res)$cluster == levels(rowData(d_se)$cluster)))
    
    rowData(res)$cluster <- factor(rowData(res)$cluster, levels = levels(rowData(d_se)$cluster))
    
    # match cells to clusters
    ix_match <- match(rowData(d_se)$cluster, rowData(res)$cluster)
    
    p_vals_clusters <- rowData(res)$P.Value
    p_adj_clusters <- rowData(res)$adj.P.Val
    
    p_vals_cells <- p_vals_clusters[ix_match]
    p_adj_cells <- p_adj_clusters[ix_match]
    
    
    # set up data frame with results and true spike-in status at cell level
    
    stopifnot(length(p_vals_cells) == length(is_spikein), 
              length(p_adj_cells) == length(is_spikein))
    
    res_p_vals <- p_vals_cells
    res_p_adj <- p_adj_cells
    
    # replace NAs (due to filtering) to ensure same cells are returned for all methods
    res_p_vals[is.na(res_p_vals)] <- 1
    res_p_adj[is.na(res_p_adj)] <- 1
    
    # return values for this condition and healthy
    res <- data.frame(p_vals = res_p_vals, 
                      p_adj = res_p_adj, 
                      spikein = is_spikein)
    
    # store results
    out_diffcyt_DA_limma_supp_separate_clustering[[th]][[j]] <- res
    
  }
}




#####################
# Save output objects
#####################

save(out_diffcyt_DA_limma_supp_separate_clustering, runtime_diffcyt_DA_limma_supp_separate_clustering, 
     file = file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_limma_supp_separate_clustering.RData"))

save(out_clusters_diffcyt_DA_limma_supp_separate_clustering, 
     file = file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_limma_supp_separate_clustering.RData"))

save(out_objects_diffcyt_DA_limma_supp_separate_clustering, 
     file = file.path(DIR_RDATA, "out_objects_AML_sim_diffcyt_DA_limma_supp_separate_clustering.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_AML_sim_diffcyt_DA_limma_supp_separate_clustering.txt"))
sessionInfo()
sink()



