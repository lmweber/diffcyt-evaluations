##########################################################################################
# Script to run methods
# 
# - method: diffcyt-DA-limma
# - data set: AML-spike-in
# 
# Lukas Weber, July 2017
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)
library(limma)




################################
# Loop to run for each threshold
################################

# spike-in thresholds (must match filenames)
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")

# lists to store objects
is_spikein_thresholds <- 
  d_se_thresholds <- d_counts_thresholds <- d_medians_all_thresholds <- 
  out_diffcyt_DA_limma <- out_diffcyt_DA_limma_sorted <- 
  vector("list", length(thresholds))

names(is_spikein_thresholds) <- 
  names(d_se_thresholds) <- names(d_counts_thresholds) <- names(d_medians_all_thresholds) <- 
  names(out_diffcyt_DA_limma) <- names(out_diffcyt_DA_limma_sorted) <- 
  thresholds




for (th in 1:length(thresholds)) {
  
  ##############################
  # Load data and pre-processing
  ##############################
  
  # ---------
  # load data
  # ---------
  
  # filenames
  files_healthy <- list.files("../../../../benchmark_data/AML_spike_in/data/healthy", 
                              pattern = "\\.fcs$", full.names = TRUE)
  files_CN <- list.files("../../../../benchmark_data/AML_spike_in/data/CN", 
                         pattern = paste0("_", thresholds[th], "\\.fcs$"), full.names = TRUE)
  files_CBF <- list.files("../../../../benchmark_data/AML_spike_in/data/CBF", 
                          pattern = paste0("_", thresholds[th], "\\.fcs$"), full.names = TRUE)
  
  # load data
  files_load <- c(files_healthy, files_CN, files_CBF)
  files_load
  
  d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)
  
  # sample IDs, group IDs, block IDs
  sample_IDs <- gsub("(_[0-9]+pc$)|(_0\\.[0-9]+pc$)", "", 
                     gsub("^AML_spike_in_", "", 
                          gsub("\\.fcs$", "", basename(files_load))))
  sample_IDs
  
  group_IDs <- gsub("_.*$", "", sample_IDs)
  group_IDs
  
  block_IDs <- gsub("^.*_", "", sample_IDs)
  block_IDs
  
  # set group_IDs reference level (for differential tests)
  group_IDs <- factor(group_IDs, levels = c("healthy", "CN", "CBF"))
  group_IDs
  
  # check all match correctly
  data.frame(sample_IDs, group_IDs, block_IDs)
  
  # indices of all marker columns, lineage markers, and functional markers
  # (16 surface markers / 15 functional markers; see Levine et al. 2015, Supplemental 
  # Information, p. 4)
  cols_markers <- 11:41
  cols_lineage <- c(35, 29, 14, 30, 12, 26, 17, 33, 41, 32, 22, 40, 27, 37, 23, 39)
  cols_func <- setdiff(cols_markers, cols_lineage)
  
  
  # ---------------------
  # store spike-in status
  # ---------------------
  
  # store spike-in status for each cell
  
  is_spikein_thresholds[[th]] <- vector("list", length(sample_IDs))
  names(is_spikein_thresholds[[th]]) <- sample_IDs
  
  for (i in 1:length(sample_IDs)) {
    exprs_i <- exprs(d_input[[i]])
    is_spikein_thresholds[[th]][[i]] <- exprs_i[, "spikein"]
  }
  
  
  
  
  #########################################################################
  # Run initial steps of 'diffcyt' pipeline (prior to differential testing)
  #########################################################################
  
  # prepare data into required format
  # (note: using lineage markers for clustering, and functional markers for DE testing)
  d_se <- prepareData(d_input, sample_IDs, group_IDs, 
                      cols_markers, cols_lineage, cols_func)
  
  # check
  colnames(d_se)[cols_lineage]
  colnames(d_se)[cols_func]
  
  
  # transform data
  d_se <- transformData(d_se, cofactor = 5)
  
  
  # clustering
  # (note: running clustering once on all samples together, including multiple conditions)
  seed <- 123
  runtime_clustering <- system.time(
    d_se <- generateClusters(d_se, xdim = 30, ydim = 30, seed = seed)
  )
  
  runtime_clustering  # ~60 sec (30x30 clusters)
  
  # check
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
  
  
  # calculate cluster medians across all samples (for plotting)
  d_medians_all <- calcMediansAll(d_se)
  
  dim(d_medians_all)
  
  
  # ------------------
  # store data objects
  # ------------------
  
  d_se_thresholds[[th]] <- d_se
  d_counts_thresholds[[th]] <- d_counts
  d_medians_all_thresholds[[th]] <- d_medians_all
  
  
  
  
  ################################################
  # Test for differentially abundant (DA) clusters
  ################################################
  
  # test separately for each condition: CN vs. healthy, CBF vs. healthy
  
  out_diffcyt_DA_limma[[th]] <- out_diffcyt_DA_limma_sorted[[th]] <- vector("list", length(cond_names))
  names(out_diffcyt_DA_limma[[th]]) <- names(out_diffcyt_DA_limma_sorted[[th]]) <- cond_names
  
  
  for (j in 1:length(cond_names)) {
    
    # ------------
    # run DA tests
    # ------------
    
    path <- paste0("../../../plots/AML_spike_in/diffcyt/", thresholds[th], "/diffcyt_DA_limma/", cond_names[j])
    
    # set up contrast
    design <- model.matrix(~ 0 + group_IDs)
    contr_str <- paste0("group_IDs", cond_names[j], " - group_IDshealthy")
    contrast <- makeContrasts(contr_str, levels = design)
    
    
    # test for DA clusters
    # (note: using 'block_IDs' argument for paired tests)
    runtime_DA <- system.time(
      res_DA <- testDA_limma(d_counts, group_IDs, contrast, 
                             block_IDs = block_IDs, path = path)
    )
    
    runtime_DA
    
    
    # -------
    # results
    # -------
    
    # show results
    rowData(res_DA)
    
    # sort to show top (most highly significant) clusters first
    res_DA_sorted <- rowData(res_DA)[order(rowData(res_DA)$adj.P.Val), ]
    
    print(head(res_DA_sorted, 10))
    #View(as.data.frame(res_DA_sorted))
    
    # number of significant DA clusters
    print(table(res_DA_sorted$adj.P.Val < 0.05))
    
    
    # --------------
    # output objects
    # --------------
    
    out_diffcyt_DA_limma[[th]][[j]] <- res_DA
    out_diffcyt_DA_limma_sorted[[th]][[j]] <- res_DA_sorted
  }
}




#####################
# Save output objects
#####################

save.image("../../../RData/AML_spike_in/outputs_AML_spike_in_diffcyt_DA_limma.RData")




#####################
# Session information
#####################

sink("../../../session_info/AML_spike_in/session_info_AML_spike_in_diffcyt_DA_limma.txt")
sessionInfo()
sink()



