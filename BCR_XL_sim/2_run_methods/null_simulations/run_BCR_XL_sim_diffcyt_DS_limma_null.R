##########################################################################################
# Script to run methods
# 
# - method: diffcyt-DS-limma
# - data set: BCR-XL-sim
# 
# - null simulations
# 
# Lukas Weber, May 2018
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)


DIR_BENCHMARK <- "../../../../../benchmark_data/BCR_XL_sim/data/null_simulations"
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/null_simulations/diagnostic/diffcyt_DS_limma"
DIR_RDATA <- "../../../../RData/BCR_XL_sim/null_simulations"
DIR_SESSION_INFO <- "../../../../session_info/BCR_XL_sim/null_simulations"




#############
# Preliminary
#############

# names of random seeds used
seed_names <- c("seed1", "seed2", "seed3")

# lists to store objects
out_diffcyt_DS_limma_null  <- 
  out_clusters_diffcyt_DS_limma_null <- 
  out_objects_diffcyt_DS_limma_null <- 
  runtime_diffcyt_DS_limma_null <- vector("list", length(seed_names))
names(out_diffcyt_DS_limma_null) <- 
  names(out_clusters_diffcyt_DS_limma_null) <- 
  names(out_objects_diffcyt_DS_limma_null) <- 
  names(runtime_diffcyt_DS_limma_null) <- seed_names




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
  
  sample_id <- gsub("^BCR_XL_sim_", "", 
                    gsub("\\.fcs$", "", basename(files_load)))
  sample_id
  
  group_id <- factor(gsub("^.*_", "", sample_id), levels = c("null1", "null2"))
  group_id
  
  patient_id <- factor(gsub("_.*$", "", sample_id))
  patient_id
  
  experiment_info <- data.frame(group_id, patient_id, sample_id)
  experiment_info
  
  # marker information
  
  # indices of all marker columns, lineage markers, and functional markers
  # (10 surface markers / 14 functional markers; see Bruggner et al. 2014, Table 1)
  cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
  cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
  cols_func <- setdiff(cols_markers, cols_lineage)
  
  marker_name <- colnames(d_input[[1]])
  marker_name <- gsub("\\(.*$", "", marker_name)
  
  marker_class <- rep("none", length(marker_name))
  marker_class[cols_lineage] <- "type"
  marker_class[cols_func] <- "state"
  marker_class <- factor(marker_class, levels = c("type", "state", "none"))
  
  # exclude CD45 from clustering
  marker_class[marker_name == "CD45"] <- "none"
  
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
    
    colnames(d_se)[colData(d_se)$marker_class == "type"]
    colnames(d_se)[colData(d_se)$marker_class == "state"]
    
    # transform data
    d_se <- transformData(d_se, cofactor = 5)
    
    # clustering
    # (runtime: ~5 sec with xdim = 10, ydim = 10)
    seed <- 1234
    d_se <- generateClusters(d_se, xdim = 10, ydim = 10, seed_clustering = seed)
    
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
  
  out_objects_diffcyt_DS_limma_null[[s]] <- list(
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
  # note: include fixed effects for 'patient_id'
  contrast_vec <- c(0, 1, 0, 0, 0, 0, 0, 0, 0)
  
  runtime_tests <- system.time({
    
    # set up design matrix
    # note: include fixed effects for 'patient_id'
    design <- createDesignMatrix(experiment_info, cols_design = 1:2)
    design
    
    # set up contrast matrix
    contrast <- createContrast(contrast_vec)
    contrast
    
    # run tests
    res <- testDS_limma(d_counts, d_medians, design, contrast, path = DIR_PLOTS)
    
  })
  
  # show results
  rowData(res)
  
  # sort to show top (most highly significant) cluster-marker combinations first
  res_sorted <- rowData(res)[order(rowData(res)$p_adj), ]
  print(head(res_sorted, 10))
  #View(as.data.frame(res_sorted))
  
  # number of significant tests (note: one test per cluster-marker combination)
  print(table(res_sorted$p_adj <= 0.1))
  
  # runtime (~30 sec on laptop)
  runtime_total <- runtime_preprocessing[["elapsed"]] + runtime_tests[["elapsed"]]
  print(runtime_total)
  
  runtime_diffcyt_DS_limma_null[[s]] <- runtime_total
  
  
  # ---------------------------------------------
  # store results at cluster level (for plotting)
  # ---------------------------------------------
  
  res_clusters <- as.data.frame(rowData(res))
  
  out_clusters_diffcyt_DS_limma_null[[s]] <- res_clusters
  
  
  
  
  ##############################
  # Return results at cell level
  ##############################
  
  # Note: diffcyt methods return results for each cluster-marker combination. To enable
  # performance comparisons between methods at the cell level, we assign the same p-values
  # to all cells within a given cluster-marker combination.
  
  # Note: return cell-level results for marker pS6 only, since the comparative evaluations
  # (ROC curves etc) are based on pS6 in B cells only.
  
  
  # number of cells per sample (including spike-in cells)
  n_cells <- sapply(d_input, nrow)
  
  # identify B cells (these contain the true differential signal; from both 'base' and
  # 'spike' conditions)
  is_B_cell <- unlist(sapply(d_input, function(d) exprs(d)[, "B_cell"]))
  stopifnot(length(is_B_cell) == sum(n_cells))
  
  
  # match cluster-level p-values for marker pS6 to individual cells
  
  stopifnot(nrow(rowData(res)) == nlevels(rowData(d_se)$cluster_id) * length(cols_func), 
            all(levels(rowData(res)$cluster_id) == levels(rowData(d_se)$cluster_id)), 
            all(levels(rowData(res)$cluster_id) %in% rowData(res)$cluster_id))
  
  # select results for pS6
  res_pS6 <- res[rowData(res)$marker == "pS6", ]
  
  # match cells to clusters
  ix_match <- match(rowData(d_se)$cluster_id, rowData(res_pS6)$cluster_id)
  
  p_vals_clusters <- rowData(res_pS6)$p_val
  p_adj_clusters <- rowData(res_pS6)$p_adj
  
  p_vals_cells <- p_vals_clusters[ix_match]
  p_adj_cells <- p_adj_clusters[ix_match]
  
  
  # set up data frame with results (for marker pS6) and true B-cell status at cell level
  
  res_p_vals <- p_vals_cells
  res_p_adj <- p_adj_cells
  
  # replace NAs (due to filtering) to ensure same cells are returned for all methods
  res_p_vals[is.na(res_p_vals)] <- 1
  res_p_adj[is.na(res_p_adj)] <- 1
  
  stopifnot(length(res_p_vals) == length(res_p_adj), 
            length(res_p_vals) == length(is_B_cell))
  
  res <- data.frame(p_val = res_p_vals, 
                    p_adj = res_p_adj, 
                    B_cell = is_B_cell)
  
  # store results
  out_diffcyt_DS_limma_null[[s]] <- res
  
}




#####################
# Save output objects
#####################

save(out_diffcyt_DS_limma_null, runtime_diffcyt_DS_limma_null, 
     file = file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_limma_null.RData"))

save(out_clusters_diffcyt_DS_limma_null, 
     file = file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DS_limma_null.RData"))

save(out_objects_diffcyt_DS_limma_null, 
     file = file.path(DIR_RDATA, "out_objects_BCR_XL_sim_diffcyt_DS_limma_null.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_BCR_XL_sim_diffcyt_DS_limma_null.txt"))
sessionInfo()
sink()



