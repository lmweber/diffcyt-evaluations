##########################################################################################
# Script to run 'diffcyt' methods for data set 'AML-spike-in'
#
# Lukas Weber, May 2017
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)



########################
# Run for each threshold
########################

# spike-in thresholds (as specified in filenames)
thresholds <- c("1pc", "0.1pc", "0.01pc")



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
  
  # sample IDs, group IDs, and block IDs (for paired tests)
  sample_IDs <- gsub("(_[0-9]+pc$)|(_0\\.[0-9]+pc$)", "", 
                     gsub("^AML_spike_in_", "", 
                          gsub("\\.fcs$", "", basename(files_load))))
  sample_IDs
  
  group_IDs <- gsub("_.*$", "", sample_IDs)
  group_IDs
  
  block_IDs <- gsub("^.*_", "", sample_IDs)
  block_IDs
  
  # check all match correctly
  sample_IDs
  group_IDs
  block_IDs
  
  # set group reference level
  group_IDs <- factor(group_IDs, levels = c("healthy", "CN", "CBF"))
  group_IDs
  
  # indices of all marker columns, lineage markers, and functional markers
  # (16 surface markers / 15 functional markers; see Levine et al. 2015, Supplemental
  # Information, p. 4)
  cols_markers <- 11:41
  cols_lineage <- c(35, 29, 14, 30, 12, 26, 17, 33, 41, 32, 22, 40, 27, 37, 23, 39)
  cols_func <- setdiff(cols_markers, cols_lineage)
  
  # prepare data into required format
  # (note: using lineage markers for clustering, and functional markers for DE testing)
  d_se <- prepareData(d_input, sample_IDs, cols_markers, cols_lineage, cols_func)
  
  # check markers
  colnames(d_se)[cols_lineage]
  colnames(d_se)[cols_func]
  
  
  # --------------
  # transform data
  # --------------
  
  # transform all marker columns
  # (using asinh with cofactor = 5; see Bendall et al. 2011, Supp. Fig. S2)
  d_se <- transformData(d_se, cofactor = 5)
  
  
  # ----------
  # clustering
  # ----------
  
  # generate mini-clusters
  seed <- 123
  runtime_clustering <- system.time(
    d_se <- generateClusters(d_se, xdim = 20, ydim = 20, seed = seed)
  )
  
  # check
  nrow(rowData(d_se))                   # number of cells
  sum(table(rowData(d_se)$cluster))
  length(table(rowData(d_se)$cluster))  # number of clusters
  
  
  # ---------------------------------------------
  # calculate cluster cell counts, medians, ECDFs
  # ---------------------------------------------
  
  # calculate cluster cell counts
  d_counts <- calcCounts(d_se)
  dim(d_counts)
  
  # calculate cluster medians
  d_medians <- calcMedians(d_se)
  dim(d_medians)
  
  # calculate ECDFs
  d_ecdfs <- calcECDFs(d_se)
  dim(d_ecdfs)
  
  # subset marker expresison values
  d_vals <- subsetVals(d_se)
  dim(d_vals)
  
  
  
  
  ###########
  # Contrasts
  ###########
  
  # currently not using contrasts; set up matrix of conditions instead (to subset objects)
  
  cond_names <- c("CNvsH", "CBFvsH")
  
  cond <- matrix(NA, nrow = length(cond_names), ncol = length(group_IDs))
  colnames(cond) <- group_IDs
  rownames(cond) <- cond_names
  
  for (i in 1:length(cond_names)) {
    cond[i, ] <- group_IDs %in% levels(group_IDs)[c(1, i + 1)]
  }
  
  # check
  cond
  
  
  
  
  ################################################
  # Test for differentially abundant (DA) clusters
  ################################################
  
  # test separately for each condition: CN vs. healthy, CBF vs. healthy
  out_DA <- out_DA_sorted <- vector("list", length(cond_names))
  names(out_DA) <- out_DA_sorted <- cond_names
  
  for (i in 1:length(cond_names)) {
    
    # subset objects
    sample_IDs_sub <- sample_IDs[cond[i, ]]
    group_IDs_sub <- group_IDs[cond[i, ]]
    group_IDs_sub <- droplevels(group_IDs_sub)
    block_IDs_sub <- block_IDs[cond[i, ]]
    
    d_counts_sub <- d_counts[, cond[i, ]]
    d_medians_sub <- d_medians[, cond[i, ]]
    d_ecdfs_sub <- d_ecdfs[, cond[i, ]]
    d_vals_sub <- d_vals[, cond[i, ]]
    
    
    # test for differentially abundant (DA) clusters
    runtime_DA <- system.time(
      res_DA <- testDA(d_counts_sub, group_IDs_sub, 
                       paired = TRUE, block_IDs = block_IDs_sub, plot = TRUE, 
                       path = paste0("../../../plots/diffcyt/AML_spike_in/", thresholds[th], "/DA/", cond_names[i]))
    )
    
    # show results
    rowData(res_DA)
    
    # sort to show top (most highly significant) clusters first
    res_DA_sorted <- rowData(res_DA)[order(rowData(res_DA)$adj.P.Val), ]
    
    head(res_DA_sorted, 10)
    #View(res_DA_sorted)
    
    
    # --------
    # analysis
    # --------
    
    # number of significant DA clusters
    print(table(res_DA_sorted$adj.P.Val < 0.05))
    
    # save output objects
    out_DA[[i]] <- res_DA
    out_DA_sorted[[i]] <- res_DA_sorted
  }
  
  
  
  
  #############################################################################
  # Test for differential expression (DE) of functional markers within clusters
  # (method 'diffcyt-med')
  #############################################################################
  
  # test separately for each condition: CN vs. healthy, CBF vs. healthy
  out_DE_med <- out_DE_med_sorted <- vector("list", length(cond_names))
  names(out_DE_med) <- out_DE_med_sorted <- cond_names
  
  for (i in 1:length(cond_names)) {
    
    # subset objects
    sample_IDs_sub <- sample_IDs[cond[i, ]]
    group_IDs_sub <- group_IDs[cond[i, ]]
    group_IDs_sub <- droplevels(group_IDs_sub)
    block_IDs_sub <- block_IDs[cond[i, ]]
    
    d_counts_sub <- d_counts[, cond[i, ]]
    d_medians_sub <- d_medians[, cond[i, ]]
    d_ecdfs_sub <- d_ecdfs[, cond[i, ]]
    d_vals_sub <- d_vals[, cond[i, ]]
    
    # test for differential expression (DE) of functional markers within clusters
    runtime_DE_med <- system.time(
      res_DE_med <- testDE_med(d_counts_sub, d_medians_sub, group_IDs_sub, 
                               paired = TRUE, block_IDs = block_IDs_sub, plot = TRUE, 
                               path = paste0("../../../plots/diffcyt/AML_spike_in/", thresholds[th], "/DE_med/", cond_names[i]))
    )
    
    # show results
    rowData(res_DE_med)
    
    # sort to show top (most highly significant) cluster-marker combinations first
    res_DE_med_sorted <- rowData(res_DE_med)[order(rowData(res_DE_med)$adj.P.Val), ]
    
    head(res_DE_med_sorted, 10)
    #View(res_DE_med_sorted)
    
    
    # --------
    # analysis
    # --------
    
    # number of significant DE cluster-marker combinations
    print(table(res_DE_med_sorted$adj.P.Val < 0.05))
    
    # save output objects
    out_DE_med[[i]] <- res_DE_med
    out_DE_med_sorted[[i]] <- res_DE_med_sorted
  }
  
  
  
  
  #############################################################################
  # Test for differential expression (DE) of functional markers within clusters
  # (method 'diffcyt-FDA')
  #############################################################################
  
  # test separately for each condition: CN vs. healthy, CBF vs. healthy
  out_DE_FDA_unwtd <- out_DE_FDA_unwtd_sorted <- vector("list", length(cond_names))
  out_DE_FDA_wtd <- out_DE_FDA_wtd_sorted <- vector("list", length(cond_names))
  names(out_DE_FDA_unwtd) <- out_DE_FDA_unwtd_sorted <- cond_names
  names(out_DE_FDA_wtd) <- out_DE_FDA_wtd_sorted <- cond_names
  
  for (i in 1:length(cond_names)) {
    
    # subset objects
    sample_IDs_sub <- sample_IDs[cond[i, ]]
    group_IDs_sub <- group_IDs[cond[i, ]]
    group_IDs_sub <- droplevels(group_IDs_sub)
    block_IDs_sub <- block_IDs[cond[i, ]]
    
    d_counts_sub <- d_counts[, cond[i, ]]
    d_medians_sub <- d_medians[, cond[i, ]]
    d_ecdfs_sub <- d_ecdfs[, cond[i, ]]
    d_vals_sub <- d_vals[, cond[i, ]]
    
    
    # ----------
    # unweighted
    # ----------
    
    # set seed (for permutation tests)
    set.seed(123)
    
    # test for differential expression (DE) of functional markers within clusters
    runtime_DE_FDA_unwtd <- system.time(
      res_DE_FDA_unwtd <- testDE_FDA(d_counts_sub, d_medians_sub, d_ecdfs_sub, group_IDs_sub, 
                                     weighted = FALSE, paired = TRUE, block_IDs = block_IDs_sub, 
                                     n_perm = 1000, n_cores = 6)
    )
    
    # show results
    rowData(res_DE_FDA_unwtd)
    
    # sort to show top (most highly significant) cluster-marker combinations first
    res_DE_FDA_unwtd_sorted <- rowData(res_DE_FDA_unwtd)[order(rowData(res_DE_FDA_unwtd)$p_adj), ]
    
    head(res_DE_FDA_unwtd_sorted, 10)
    #View(res_DE_FDA_unwtd_sorted)
    
    
    # --------
    # weighted
    # --------
    
    # set seed (for permutation tests)
    set.seed(123)
    
    # test for differential expression (DE) of functional markers within clusters
    runtime_DE_FDA_wtd <- system.time(
      res_DE_FDA_wtd <- testDE_FDA(d_counts_sub, d_medians_sub, d_ecdfs_sub, group_IDs_sub, 
                                   weighted = TRUE, paired = TRUE, block_IDs = block_IDs_sub, 
                                   n_perm = 1000, n_cores = 6)
    )
    
    # show results
    rowData(res_DE_FDA_wtd)
    
    # sort to show top (most highly significant) cluster-marker combinations first
    res_DE_FDA_wtd_sorted <- rowData(res_DE_FDA_wtd)[order(rowData(res_DE_FDA_wtd)$p_adj), ]
    
    head(res_DE_FDA_wtd_sorted, 10)
    #View(res_DE_FDA_wtd_sorted)
    
    
    # --------
    # analysis
    # --------
    
    # number of significant DE cluster-marker combinations
    
    # unweighted
    print(table(res_DE_FDA_unwtd_sorted$p_adj < 0.05))
    
    # weighted
    print(table(res_DE_FDA_wtd_sorted$p_adj < 0.05))
    
    # save output objects
    out_DE_FDA_unwtd[[i]] <- res_DE_FDA_unwtd
    out_DE_FDA_unwtd_sorted[[i]] <- res_DE_FDA_unwtd_sorted
    out_DE_FDA_wtd[[i]] <- res_DE_FDA_wtd
    out_DE_FDA_wtd_sorted[[i]] <- res_DE_FDA_wtd_sorted
  }
  
  
  
  
  #############################################################################
  # Test for differential expression (DE) of functional markers within clusters
  # (method 'diffcyt-KS')
  #############################################################################
  
  # test separately for each condition: CN vs. healthy, CBF vs. healthy
  out_DE_KS_paired <- out_DE_KS_paired_sorted <- vector("list", length(cond_names))
  out_DE_KS_unpaired <- out_DE_KS_unpaired_sorted <- vector("list", length(cond_names))
  names(out_DE_KS_paired) <- out_DE_KS_paired_sorted <- cond_names
  names(out_DE_KS_unpaired) <- out_DE_KS_unpaired_sorted <- cond_names
  
  for (i in 1:length(cond_names)) {
    
    # subset objects
    sample_IDs_sub <- sample_IDs[cond[i, ]]
    group_IDs_sub <- group_IDs[cond[i, ]]
    group_IDs_sub <- droplevels(group_IDs_sub)
    block_IDs_sub <- block_IDs[cond[i, ]]
    
    d_counts_sub <- d_counts[, cond[i, ]]
    d_medians_sub <- d_medians[, cond[i, ]]
    d_ecdfs_sub <- d_ecdfs[, cond[i, ]]
    d_vals_sub <- d_vals[, cond[i, ]]
    
    
    # ------
    # paired
    # ------
    
    # test for differential expression (DE) of functional markers within clusters
    runtime_DE_KS_paired <- system.time(
      res_DE_KS_paired <- testDE_KS(d_counts_sub, d_medians_sub, d_vals_sub, group_IDs_sub, 
                                    paired = TRUE, block_IDs = block_IDs_sub)
    )
    
    # show results
    rowData(res_DE_KS_paired)
    
    # sort to show top (most highly significant) cluster-marker combinations first
    res_DE_KS_paired_sorted <- rowData(res_DE_KS_paired)[order(rowData(res_DE_KS_paired)$p_adj), ]
    
    head(res_DE_KS_paired_sorted, 10)
    #View(res_DE_KS_paired_sorted)
    
    
    # --------
    # unpaired
    # --------
    
    # test for differential expression (DE) of functional markers within clusters
    runtime_DE_KS_unpaired <- system.time(
      res_DE_KS_unpaired <- testDE_KS(d_counts_sub, d_medians_sub, d_vals_sub, group_IDs_sub, 
                                      n_perm = 1000, n_cores = 6)
    )
    
    # show results
    rowData(res_DE_KS_unpaired)
    
    # sort to show top (most highly significant) cluster-marker combinations first
    res_DE_KS_unpaired_sorted <- rowData(res_DE_KS_unpaired)[order(rowData(res_DE_KS_unpaired)$p_adj), ]
    
    head(res_DE_KS_unpaired_sorted, 10)
    #View(res_DE_KS_unpaired_sorted)
    
    
    # --------
    # analysis
    # --------
    
    # number of significant DE cluster-marker combinations
    
    # paired
    print(table(res_DE_KS_paired_sorted$p_adj < 0.05))
    
    # unpaired
    print(table(res_DE_KS_unpaired_sorted$p_adj < 0.05))
    
    # save output objects
    out_DE_KS_paired[[i]] <- res_DE_KS_paired
    out_DE_KS_paired_sorted[[i]] <- res_DE_KS_paired_sorted
    out_DE_KS_unpaired[[i]] <- res_DE_KS_unpaired
    out_DE_KS_unpaired_sorted[[i]] <- res_DE_KS_unpaired_sorted
  }
  
  
  
  
  #############################################################################
  # Test for differential expression (DE) of functional markers within clusters
  # (method 'diffcyt-LM')
  #############################################################################
  
  # test separately for each condition: CN vs. healthy, CBF vs. healthy
  out_DE_LM <- out_DE_LM_sorted <- vector("list", length(cond_names))
  names(out_DE_LM) <- out_DE_LM_sorted <- cond_names
  
  for (i in 1:length(cond_names)) {
    
    # subset objects
    sample_IDs_sub <- sample_IDs[cond[i, ]]
    group_IDs_sub <- group_IDs[cond[i, ]]
    group_IDs_sub <- droplevels(group_IDs_sub)
    block_IDs_sub <- block_IDs[cond[i, ]]
    
    d_counts_sub <- d_counts[, cond[i, ]]
    d_medians_sub <- d_medians[, cond[i, ]]
    d_ecdfs_sub <- d_ecdfs[, cond[i, ]]
    d_vals_sub <- d_vals[, cond[i, ]]
    
    
    # test for differential expression (DE) of functional markers within clusters
    runtime_DE_LM <- system.time(
      res_DE_LM <- testDE_LM(d_counts_sub, d_medians_sub, d_ecdfs_sub, group_IDs_sub, 
                             paired = TRUE, block_IDs = block_IDs_sub)
    )
    
    # show results
    rowData(res_DE_LM)
    
    # sort to show top (most highly significant) cluster-marker combinations first
    res_DE_LM_sorted <- rowData(res_DE_LM)[order(rowData(res_DE_LM)$p_adj), ]
    
    head(res_DE_LM_sorted, 10)
    #View(res_DE_LM_sorted)
    
    
    # --------
    # analysis
    # --------
    
    # number of significant DE cluster-marker combinations
    print(table(res_DE_LM_sorted$p_adj < 0.05))
    
    # save output objects
    out_DE_LM[[i]] <- res_DE_LM
    out_DE_LM_sorted[[i]] <- res_DE_LM_sorted
  }
}




#####################
# Save output objects
#####################

save.image("../../../RData/outputs_diffcyt_AML_spike_in.RData")




#####################
# Session information
#####################

sink("../../../session_info/session_info_AML_spike_in.txt")
sessionInfo()
sink()



