##########################################################################################
# Script to run methods
# 
# - method: diffcyt-DA-GLMM
# - data set: AML-sim
# 
# - supplementary results: varying random seeds for generating benchmark data
# 
# Lukas Weber, February 2018
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)


DIR_BENCHMARK <- "../../../../../benchmark_data/AML_sim/data/random_seeds"
DIR_RDATA <- "../../../../RData/AML_sim/supp_random_seeds_data"
DIR_SESSION_INFO <- "../../../../session_info/AML_sim/supp_random_seeds_data"




################################
# Loop to run for each threshold
################################

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# condition names
cond_names <- c("CN", "CBF")

# names of random seeds used
seed_names <- c("seed1", "seed2", "seed3")

# contrasts (to compare each of 'CN' and 'CBF' vs. 'healthy')
# note: include random effects for 'patient_IDs' and 'sample_IDs'
contrasts_list <- list(CN = c(0, 1, 0), CBF = c(0, 0, 1))

# lists to store objects and runtime
out_diffcyt_DA_GLMM_supp_random_seeds_data <- runtime_diffcyt_DA_GLMM_supp_random_seeds_data <- 
  out_clusters_diffcyt_DA_GLMM_supp_random_seeds_data <- out_objects_diffcyt_DA_GLMM_supp_random_seeds_data <- 
  vector("list", length(seed_names))
names(out_diffcyt_DA_GLMM_supp_random_seeds_data) <- names(runtime_diffcyt_DA_GLMM_supp_random_seeds_data) <- 
  names(out_clusters_diffcyt_DA_GLMM_supp_random_seeds_data) <- names(out_objects_diffcyt_DA_GLMM_supp_random_seeds_data) <- 
  seed_names




for (s in 1:length(seed_names)) {
  
  out_diffcyt_DA_GLMM_supp_random_seeds_data[[s]] <- runtime_diffcyt_DA_GLMM_supp_random_seeds_data[[s]] <- 
    out_clusters_diffcyt_DA_GLMM_supp_random_seeds_data[[s]] <- out_objects_diffcyt_DA_GLMM_supp_random_seeds_data[[s]] <- 
    vector("list", length(thresholds))
  names(out_diffcyt_DA_GLMM_supp_random_seeds_data[[s]]) <- names(runtime_diffcyt_DA_GLMM_supp_random_seeds_data[[s]]) <- 
    names(out_clusters_diffcyt_DA_GLMM_supp_random_seeds_data[[s]]) <- names(out_objects_diffcyt_DA_GLMM_supp_random_seeds_data[[s]]) <- 
    thresholds
  
  
  for (th in 1:length(thresholds)) {
    
    ###########################
    # Load data, pre-processing
    ###########################
    
    # filenames
    
    files_healthy <- list.files(file.path(DIR_BENCHMARK, seed_names[s], "healthy"), 
                                pattern = "\\.fcs$", full.names = TRUE)
    files_CN <- list.files(file.path(DIR_BENCHMARK, seed_names[s], "CN"), 
                           pattern = paste0("_", thresholds[th], "_randomseed[0-9]+\\.fcs$"), full.names = TRUE)
    files_CBF <- list.files(file.path(DIR_BENCHMARK, seed_names[s], "CBF"), 
                            pattern = paste0("_", thresholds[th], "_randomseed[0-9]+\\.fcs$"), full.names = TRUE)
    
    files_load <- c(files_healthy, files_CN, files_CBF)
    files_load
    
    # load data
    
    d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)
    
    # sample IDs, group IDs, patient IDs
    sample_IDs <- gsub("(_[0-9]+pc$)|(_0\\.[0-9]+pc$)", "", 
                       gsub("^AML_sim_", "", 
                            gsub("_randomseed[0-9]+$", "", 
                                 gsub("\\.fcs$", "", basename(files_load)))))
    sample_IDs
    
    group_IDs <- factor(gsub("_.*$", "", sample_IDs), levels = c("healthy", "CN", "CBF"))
    group_IDs
    
    patient_IDs <- factor(gsub("^.*_", "", sample_IDs))
    patient_IDs
    
    sample_info <- data.frame(group_IDs, patient_IDs, sample_IDs)
    sample_info
    
    # marker information
    
    # indices of all marker columns, lineage markers, and functional markers
    # (16 surface markers / 15 functional markers; see Levine et al. 2015, Supplemental 
    # Information, p. 4)
    cols_markers <- 11:41
    cols_lineage <- c(35, 29, 14, 30, 12, 26, 17, 33, 41, 32, 22, 40, 27, 37, 23, 39)
    cols_func <- setdiff(cols_markers, cols_lineage)
    
    stopifnot(all(sapply(seq_along(d_input), function(i) all(colnames(d_input[[i]]) == colnames(d_input[[1]])))))
    
    marker_names <- colnames(d_input[[1]])
    marker_names <- gsub("\\(.*$", "", marker_names)
    
    is_marker <- is_celltype_marker <- is_state_marker <- rep(FALSE, length(marker_names))
    
    is_marker[cols_markers] <- TRUE
    is_celltype_marker[cols_lineage] <- TRUE
    is_state_marker[cols_func] <- TRUE
    
    marker_info <- data.frame(marker_names, is_marker, is_celltype_marker, is_state_marker)
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
      
      colnames(d_se)[is_celltype_marker]
      colnames(d_se)[is_state_marker]
      
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
      
      # calculate cluster medians by sample
      d_medians <- calcMedians(d_se)
      
      dim(d_medians)
      rowData(d_medians)
      length(assays(d_medians))
      names(assays(d_medians))
      
      # calculate cluster medians across all samples
      d_medians_all <- calcMediansAll(d_se)
      
      dim(d_medians_all)
      length(assays(d_medians_all))
      
    })
    
    
    # ---------------------------------
    # store data objects (for plotting)
    # ---------------------------------
    
    out_objects_diffcyt_DA_GLMM_supp_random_seeds_data[[s]][[th]] <- list(
      d_se = d_se, 
      d_counts = d_counts, 
      d_medians = d_medians, 
      d_medians_all = d_medians_all
    )
    
    
    # -----------------------------------------
    # test for differentially abundant clusters
    # -----------------------------------------
    
    # note: test separately for each condition: CN vs. healthy, CBF vs. healthy
    
    out_diffcyt_DA_GLMM_supp_random_seeds_data[[s]][[th]] <- runtime_diffcyt_DA_GLMM_supp_random_seeds_data[[s]][[th]] <- 
      out_clusters_diffcyt_DA_GLMM_supp_random_seeds_data[[s]][[th]] <- vector("list", length(cond_names))
    names(out_diffcyt_DA_GLMM_supp_random_seeds_data[[s]][[th]]) <- names(runtime_diffcyt_DA_GLMM_supp_random_seeds_data[[s]][[th]]) <- 
      names(out_clusters_diffcyt_DA_GLMM_supp_random_seeds_data[[s]][[th]]) <- cond_names
    
    
    for (j in 1:length(cond_names)) {
      
      runtime_j <- system.time({
        
        # set up model formula
        # note: order of samples has changed
        sample_info_ordered <- as.data.frame(colData(d_counts))
        sample_info_ordered
        # note: include random effects for 'patient_IDs' and 'sample_IDs'
        formula <- createFormula(sample_info_ordered, cols_fixed = 1, cols_random = 2:3)
        formula
        
        # set up contrast matrix
        contrast <- createContrast(contrasts_list[[j]])
        contrast
        
        # run tests
        # note: adjust filtering parameter 'min_samples' (since there are 3 conditions)
        res <- testDA_GLMM(d_counts, formula, contrast, 
                           min_cells = 3, min_samples = nrow(sample_info_ordered) / 3)
        
      })
      
      # show results
      rowData(res)
      
      # sort to show top (most highly significant) clusters first
      res_sorted <- rowData(res)[order(rowData(res)$p_adj), ]
      print(head(res_sorted, 10))
      #View(as.data.frame(res_sorted))
      
      # number of significant tests (note: one test per cluster)
      print(table(res_sorted$p_adj <= 0.1))
      
      # runtime
      runtime_total <- runtime_preprocessing[["elapsed"]] + runtime_j[["elapsed"]]
      print(runtime_total)
      
      runtime_diffcyt_DA_GLMM_supp_random_seeds_data[[s]][[th]][[j]] <- runtime_total
      
      
      # ---------------------------------------------
      # store results at cluster level (for plotting)
      # ---------------------------------------------
      
      res_clusters <- as.data.frame(rowData(res))
      
      out_clusters_diffcyt_DA_GLMM_supp_random_seeds_data[[s]][[th]][[j]] <- res_clusters
      
      
      
      
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
      ix_keep_cnd <- group_IDs %in% c("healthy", cond_names[j])
      
      
      # match cluster-level p-values to individual cells
      
      stopifnot(nrow(rowData(res)) == length(levels(rowData(d_se)$cluster)), 
                all(rowData(res)$cluster == levels(rowData(d_se)$cluster)))
      
      rowData(res)$cluster <- factor(rowData(res)$cluster, levels = levels(rowData(d_se)$cluster))
      
      # match cells to clusters
      ix_match <- match(rowData(d_se)$cluster, rowData(res)$cluster)
      
      p_vals_clusters <- rowData(res)$p_vals
      p_adj_clusters <- rowData(res)$p_adj
      
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
      out_diffcyt_DA_GLMM_supp_random_seeds_data[[s]][[th]][[j]] <- res
      
    }
  }
}




#####################
# Save output objects
#####################

save(out_diffcyt_DA_GLMM_supp_random_seeds_data, runtime_diffcyt_DA_GLMM_supp_random_seeds_data, 
     file = file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_GLMM_supp_random_seeds_data.RData"))

save(out_clusters_diffcyt_DA_GLMM_supp_random_seeds_data, 
     file = file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_GLMM_supp_random_seeds_data.RData"))

save(out_objects_diffcyt_DA_GLMM_supp_random_seeds_data, 
     file = file.path(DIR_RDATA, "out_objects_AML_sim_diffcyt_DA_GLMM_supp_random_seeds_data.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_AML_sim_diffcyt_DA_GLMM_supp_random_seeds_data.txt"))
sessionInfo()
sink()



