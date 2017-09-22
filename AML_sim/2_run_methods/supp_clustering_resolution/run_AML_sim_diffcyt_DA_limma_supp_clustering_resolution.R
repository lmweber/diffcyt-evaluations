##########################################################################################
# Script to run methods
# 
# - method: diffcyt-DA-limma
# - data set: AML-sim
# 
# - supplementary results: clustering resolution
# 
# Lukas Weber, September 2017
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

# clustering resolution: define grid size for FlowSOM (e.g. 30x30 grid)
resolution <- c(3, 5, 10, 20, 30, 40, 50)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")

# contrasts (to compare each of 'CN' and 'CBF' vs. 'healthy')
contrasts_list <- list(CN = c(0, 1, 0), CBF = c(0, 0, 1))

# lists to store objects
out_diffcyt_DA_limma_supp_resolution <- runtime_diffcyt_DA_limma_supp_resolution <- vector("list", length(resolution))
names(out_diffcyt_DA_limma_supp_resolution) <- names(runtime_diffcyt_DA_limma_supp_resolution) <- resolution




for (k in 1:length(resolution)) {
  
  out_diffcyt_DA_limma_supp_resolution[[k]] <- runtime_diffcyt_DA_limma_supp_resolution[[k]] <- vector("list", length(thresholds))
  names(out_diffcyt_DA_limma_supp_resolution[[k]]) <- names(runtime_diffcyt_DA_limma_supp_resolution[[k]]) <- thresholds
  
  
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
    
    # load data
    files_load <- c(files_healthy, files_CN, files_CBF)
    files_load
    
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
    
    # check
    data.frame(sample_IDs, group_IDs, patient_IDs)
    
    # indices of all marker columns, lineage markers, and functional markers
    # (16 surface markers / 15 functional markers; see Levine et al. 2015, Supplemental 
    # Information, p. 4)
    cols_markers <- 11:41
    cols_lineage <- c(35, 29, 14, 30, 12, 26, 17, 33, 41, 32, 22, 40, 27, 37, 23, 39)
    cols_func <- setdiff(cols_markers, cols_lineage)
    
    
    # ------------------------------------
    # choose markers to use for clustering
    # ------------------------------------
    
    cols_clustering <- cols_lineage
    
    
    
    
    ##################
    # diffcyt pipeline
    ##################
    
    # --------------------
    # pre-processing steps
    # --------------------
    
    runtime_preprocessing <- system.time({
      
      # prepare data into required format
      d_se <- prepareData(d_input, sample_IDs, group_IDs, 
                          cols_markers, cols_clustering, cols_func)
      
      colnames(d_se)[cols_clustering]
      colnames(d_se)[cols_func]
      
      # transform data
      d_se <- transformData(d_se, cofactor = 5)
      
      # clustering
      # (runtime: ~60 sec for 30x30 clusters)
      # (note: clustering all samples together)
      seed <- 123
      d_se <- generateClusters(d_se, xdim = resolution[k], ydim = resolution[k], seed = seed)
      
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
      
    })
    
    # following steps not included in runtime since not required for differential testing
    
    # calculate cluster medians by sample
    d_medians <- calcMedians(d_se)
    
    dim(d_medians)
    rowData(d_medians)
    length(assays(d_medians))
    names(assays(d_medians))
    
    # calculate cluster medians across all samples
    d_medians_all <- calcMediansAll(d_se)
    
    dim(d_medians_all)
    
    
    # ----------------------------------------------
    # test for differentially abundant (DA) clusters
    # ----------------------------------------------
    
    # note: test separately for each condition: CN vs. healthy, CBF vs. healthy
    
    out_diffcyt_DA_limma_supp_resolution[[k]][[th]] <- runtime_diffcyt_DA_limma_supp_resolution[[k]][[th]] <- vector("list", length(cond_names))
    names(out_diffcyt_DA_limma_supp_resolution[[k]][[th]]) <- names(runtime_diffcyt_DA_limma_supp_resolution[[k]][[th]]) <- cond_names
    
    
    for (j in 1:length(cond_names)) {
      
      runtime_j <- system.time({
        
        # set up design matrix
        design <- createDesignMatrix(group_IDs)
        design
        
        # set up contrast matrix
        contrast <- createContrast(group_IDs, contrast = contrasts_list[[j]])
        contrast
        
        # run tests
        # - note: include 'patient_IDs' as random effects using limma 'duplicateCorrelation' methodology
        res <- testDA_limma(d_counts, design, contrast, 
                            block_IDs = patient_IDs, plot = FALSE)
        
      })
      
      # show results
      rowData(res)
      
      # sort to show top (most highly significant) clusters first
      res_sorted <- rowData(res)[order(rowData(res)$adj.P.Val), ]
      print(head(res_sorted, 10))
      #View(as.data.frame(res_sorted))
      
      # number of significant DA clusters
      print(table(res_sorted$adj.P.Val <= 0.05))
      
      # runtime (~2 min on laptop)
      runtime_total <- runtime_preprocessing[["elapsed"]] + runtime_j[["elapsed"]]
      print(runtime_total)
      
      runtime_diffcyt_DA_limma_supp_resolution[[k]][[th]][[j]] <- runtime_total
      
      
      
      
      ##############################
      # Return results at cell level
      ##############################
      
      # Note: diffcyt methods return results at cluster level (e.g. 900 small clusters). To
      # enable performance comparisons between methods at the cell level, we assign the
      # cluster-level p-values to all cells within each cluster.
      
      
      # number of cells per sample (including spike-in cells)
      n_cells <- sapply(d_input, nrow)
      
      # spike-in status for each cell
      is_spikein <- unlist(sapply(d_input, function(d) exprs(d)[, "spikein"]))
      stopifnot(length(is_spikein) == sum(n_cells))
      
      # select samples for this condition
      ix_keep_cnd <- group_IDs == cond_names[j]
      
      
      # match cluster-level p-values to individual cells
      
      stopifnot(nrow(rowData(res)) == length(levels(rowData(d_se)$cluster)), 
                all(rowData(res)$cluster == levels(rowData(d_se)$cluster)))
      
      rowData(res)$cluster <- factor(rowData(res)$cluster, levels = levels(rowData(d_se)$cluster))
      
      ix_match <- match(rowData(d_se)$cluster, rowData(res)$cluster)
      
      p_vals_clusters <- rowData(res)$P.Value
      p_adj_clusters <- rowData(res)$adj.P.Val
      
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
      
      # return values for this condition only
      res <- data.frame(p_vals = res_p_vals, 
                        p_adj = res_p_adj, 
                        spikein = is_spikein_cnd)
      
      # store results
      out_diffcyt_DA_limma_supp_resolution[[k]][[th]][[j]] <- res
      
    }
  }
}




#####################
# Save output objects
#####################

save(out_diffcyt_DA_limma_supp_resolution, runtime_diffcyt_DA_limma_supp_resolution, 
     file = file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_limma_supp_clustering_resolution.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_AML_sim_diffcyt_DA_limma_supp_clustering_resolution.txt"))
sessionInfo()
sink()



