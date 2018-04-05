##########################################################################################
# Script to run methods
# 
# - method: diffcyt-DA-edgeR
# - data set: AML-sim
# 
# - null simulations
# 
# Lukas Weber, April 2018
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)


DIR_BENCHMARK <- "../../../../../benchmark_data/AML_sim/data/null_simulations"
DIR_RDATA <- "../../../../RData/AML_sim/null_simulations"
DIR_SESSION_INFO <- "../../../../session_info/AML_sim/null_simulations"




################################
# Loop to run for each threshold
################################

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# names of random seeds used
seed_names <- c("seed1", "seed2", "seed3")

# lists to store objects
out_diffcyt_DA_edgeR_null  <- 
  out_clusters_diffcyt_DA_edgeR_null <- 
  out_objects_diffcyt_DA_edgeR_null <- 
  runtime_diffcyt_DA_edgeR_null <- vector("list", length(thresholds))
names(out_diffcyt_DA_edgeR_null) <- 
  names(out_clusters_diffcyt_DA_edgeR_null) <- 
  names(out_objects_diffcyt_DA_edgeR_null) <- 
  names(runtime_diffcyt_DA_edgeR_null) <- thresholds




for (th in 1:length(thresholds)) {
  
  # lists to store objects
  out_diffcyt_DA_edgeR_null[[th]]  <- 
    out_clusters_diffcyt_DA_edgeR_null[[th]] <- 
    out_objects_diffcyt_DA_edgeR_null[[th]] <- 
    runtime_diffcyt_DA_edgeR_null[[th]] <- vector("list", length(seed_names))
  names(out_diffcyt_DA_edgeR_null[[th]]) <- 
    names(out_clusters_diffcyt_DA_edgeR_null[[th]]) <- 
    names(out_objects_diffcyt_DA_edgeR_null[[th]]) <- 
    names(runtime_diffcyt_DA_edgeR_null[[th]]) <- seed_names
  
  
  for (s in 1:length(seed_names)) {
    
    
    ###########################
    # Load data, pre-processing
    ###########################
    
    # note: load data from each random seed
    
    
    # filenames
    
    files_null1 <- list.files(file.path(DIR_BENCHMARK, seed_names[s], "null1"), pattern = "\\.fcs$", full.names = TRUE)
    files_null2 <- list.files(file.path(DIR_BENCHMARK, seed_names[s], "null2"), pattern = "\\.fcs$", full.names = TRUE)
    
    files_load <- c(files_null1, files_null2)
    files_load
    
    # load data
    
    d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)
    
    # sample information
    
    sample_IDs <- gsub("^AML_sim_", "", 
                       gsub("\\.fcs$", "", basename(files_load)))
    sample_IDs
    
    group_IDs <- factor(gsub("^.*_", "", sample_IDs), levels = c("null1", "null2"))
    group_IDs
    
    patient_IDs <- factor(gsub("_.*$", "", sample_IDs))
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
    
    # --------------------
    # pre-processing steps
    # --------------------
    
    runtime_preprocessing <- system.time({
      
      # prepare data into required format
      d_se <- prepareData(d_input, sample_info, marker_info)
      
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
    
    out_objects_diffcyt_DA_edgeR_null[[th]][[s]] <- list(
      d_se = d_se, 
      d_counts = d_counts, 
      d_medians = d_medians, 
      d_medians_by_cluster_marker = d_medians_by_cluster_marker, 
      d_medians_by_sample_marker = d_medians_by_sample_marker
    )
    
    
    # --------------------------------------------
    # test for differential states within clusters
    # --------------------------------------------
    
    # contrast (to compare 'null2' vs. 'null1')
    # note: include zeros for 'patient'
    contrast_vec <- c(0, 1, 0, 0, 0, 0)
    
    runtime_tests <- system.time({
      
      # set up design matrix
      # note: include fixed effects for 'patient'
      design <- createDesignMatrix(sample_info, cols_include = 1:2)
      design
      
      # set up contrast matrix
      contrast <- createContrast(contrast_vec)
      contrast
      
      # run tests
      # note: use default filtering parameter 'min_samples' (since there are 2 conditions null comparison)
      res <- testDA_edgeR(d_counts, design, contrast)
      
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
    runtime_total <- runtime_preprocessing[["elapsed"]] + runtime_tests[["elapsed"]]
    print(runtime_total)
    
    runtime_diffcyt_DA_edgeR_null[[th]][[s]] <- runtime_total
    
    
    # ---------------------------------------------
    # store results at cluster level (for plotting)
    # ---------------------------------------------
    
    res_clusters <- as.data.frame(rowData(res))
    
    out_clusters_diffcyt_DA_edgeR_null[[th]][[s]] <- res_clusters
    
    
    
    
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
    
    
    # match cluster-level p-values to individual cells
    
    stopifnot(nrow(rowData(res)) == length(levels(rowData(d_se)$cluster)), 
              all(rowData(res)$cluster == levels(rowData(d_se)$cluster)))
    
    rowData(res)$cluster <- factor(rowData(res)$cluster, levels = levels(rowData(d_se)$cluster))
    
    # match cells to clusters
    ix_match <- match(rowData(d_se)$cluster, rowData(res)$cluster)
    
    p_vals_clusters <- rowData(res)$PValue
    p_adj_clusters <- rowData(res)$FDR
    
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
    out_diffcyt_DA_edgeR_null[[th]][[s]] <- res
    
  }
}




#####################
# Save output objects
#####################

save(out_diffcyt_DA_edgeR_null, runtime_diffcyt_DA_edgeR_null, 
     file = file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_edgeR_null.RData"))

save(out_clusters_diffcyt_DA_edgeR_null, 
     file = file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_edgeR_null.RData"))

save(out_objects_diffcyt_DA_edgeR_null, 
     file = file.path(DIR_RDATA, "out_objects_AML_sim_diffcyt_DA_edgeR_null.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_AML_sim_diffcyt_DA_edgeR_null.txt"))
sessionInfo()
sink()



