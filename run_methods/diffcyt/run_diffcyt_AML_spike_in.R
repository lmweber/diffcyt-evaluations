##########################################################################################
# Script to run 'diffcyt' methods for data set 'AML-spike-in'
#
# Lukas Weber, June 2017
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)



########################
# Run for each threshold
########################

# spike-in thresholds (as specified in filenames)
thresholds <- c("1pc", "0.1pc", "0.01pc")


# lists for clustering performance results
clustering_pr <- clustering_re <- clustering_F1 <- vector("list", length(thresholds))
names(clustering_pr) <- names(clustering_re) <- names(clustering_F1) <- thresholds

# lists for DA test results
out_DA <- out_DA_sorted <- vector("list", length(thresholds))
names(out_DA) <- names(out_DA_sorted) <- thresholds



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
  
  # set group reference level
  group_IDs <- factor(group_IDs, levels = c("healthy", "CN", "CBF"))
  group_IDs
  
  # check all match correctly
  data.frame(sample_IDs, group_IDs)
  
  # indices of all marker columns, lineage markers, and functional markers
  # (16 surface markers / 15 functional markers; see Levine et al. 2015, Supplemental 
  # Information, p. 4)
  cols_markers <- 11:41
  cols_lineage <- c(35, 29, 14, 30, 12, 26, 17, 33, 41, 32, 22, 40, 27, 37, 23, 39)
  cols_func <- setdiff(cols_markers, cols_lineage)
  
  
  # ------------------------------------------------------
  # remove spike-in indicator columns and store separately
  # ------------------------------------------------------
  
  # 'AML-spike-in' data set include indicator columns for spike-in cells in conditions
  # 'CN' and 'CBF'; these need to be removed and stored separately
  non_healthy <- which(group_IDs %in% c("CN", "CBF"))
  is_spikein <- vector("list", length(non_healthy))
  names(is_spikein) <- sample_IDs[non_healthy]
  
  for (i in 1:length(non_healthy)) {
    is_spikein[[i]] <- exprs(d_input[[non_healthy[i]]])[, "spikein"]
    exprs_i <- exprs(d_input[[non_healthy[i]]])
    exprs(d_input[[non_healthy[i]]]) <- exprs_i[, -ncol(exprs_i)]
  }
  
  
  # ------------
  # prepare data
  # ------------
  
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
    d_se <- generateClusters(d_se, xdim = 30, ydim = 30, seed = seed)
  )
  
  runtime_clustering  # ~300 sec (40x40 clusters); 170 sec (30x30 clusters); 75 sec (20x20 clusters)
  
  # check
  nrow(rowData(d_se))                   # number of cells
  sum(table(rowData(d_se)$cluster))
  length(table(rowData(d_se)$cluster))  # number of clusters
  
  
  # -----------------------------
  # calculate cluster cell counts
  # -----------------------------
  
  # calculate cluster cell counts
  d_counts <- calcCounts(d_se)
  dim(d_counts)
  
  # calculate cluster medians
  d_medians <- calcMedians(d_se)
  dim(d_medians)
  
  
  
  ########################
  # Contrasts / conditions
  ########################
  
  # currently not using formal contrasts; set up matrix of conditions instead
  
  cond_names <- c("CNvsH", "CBFvsH")
  
  cond <- matrix(NA, nrow = length(cond_names), ncol = length(group_IDs))
  colnames(cond) <- group_IDs
  rownames(cond) <- cond_names
  
  for (i in 1:length(cond_names)) {
    cond[i, ] <- group_IDs %in% levels(group_IDs)[c(1, i + 1)]
  }
  
  # cond <- rbind(cond, nonH = xor(cond[1, ], cond[2, ]))
  
  # check
  cond
  
  
  
  #####################################################
  # Evaluate clustering performance for spiked-in cells
  #####################################################
  
  # number of spiked-in cells per sample
  n_spikein_tbl <- sapply(is_spikein, table)
  n_spikein_tbl
  
  n_healthy <- sapply(d_input, nrow)[group_IDs == "healthy"]
  n_healthy
  
  spikein_rep <- c(rep(0, sum(n_healthy)), 
                    unlist(apply(n_spikein_tbl, 2, function(col) rep(c(0, 1), col))))
  names(spikein_rep) <- NULL
  
  group_IDs_rep <- rep(group_IDs, c(n_healthy, colSums(n_spikein_tbl)))
  
  # store additional information (group IDs, spike-in status) in 'rowData' of 'd_se'
  rowData(d_se) <- cbind(rowData(d_se), 
                         data.frame(spikein = spikein_rep), 
                         data.frame(group_IDs = group_IDs_rep))
  
  
  # loop over conditions
  th_pr <- th_re <- th_F1 <- vector("list", length(cond_names))
  names(th_pr) <- names(th_re) <- names(th_F1) <- rownames(cond)
  
  for (j in 1:length(cond_names)) {
    
    # match cluster label for each sample
    d_split <- split(rowData(d_se), rowData(d_se)$sample)
    
    labels_matched <- sapply(d_split, function(d) unname(which.max(table(d[d$spikein == 1, ]$cluster))))
    labels_matched[group_IDs == "healthy"] <- NA
    labels_matched
    
    
    n_matched_tot <- sapply(labels_matched, function(l) sum(rowData(d_se)$cluster == l))
    n_matched_tot
    
    n_matched_samp <- mapply(function(d, l) {
      sum(d$cluster == l)
    }, d_split, labels_matched)
    n_matched_samp
    
    n_matched_samp_correct <- mapply(function(d, l) {
      sum(d$cluster == l & d$spikein == 1)
    }, d_split, labels_matched)
    n_matched_samp_correct
    
    n_spikein <- c(rep(NA, sum(group_IDs == "healthy")), n_spikein_tbl[2, ])
    n_spikein
    
    
    # calculate precision, recall, F1 score
    pr <- n_matched_samp_correct / n_matched_samp
    re <- n_matched_samp_correct / n_spikein
    F1 <- 2 * (pr * re) / (pr + re)
    
    # store
    th_pr[[j]] <- pr
    th_re[[j]] <- re
    th_F1[[j]] <- F1
  }
  
  # store
  clustering_pr[[th]] <- th_pr
  clustering_re[[th]] <- th_re
  clustering_F1[[th]] <- th_F1
  
  
  
  ################################################
  # Test for differentially abundant (DA) clusters
  ################################################
  
  # test separately for each condition: CN vs. healthy, CBF vs. healthy
  out_DA[[th]] <- out_DA_sorted[[th]] <- vector("list", length(cond_names))
  names(out_DA[[th]]) <- names(out_DA_sorted[[th]]) <- cond_names
  
  for (j in 1:length(cond_names)) {
    
    # subset objects
    sample_IDs_sub <- sample_IDs[cond[j, ]]
    group_IDs_sub <- group_IDs[cond[j, ]]
    group_IDs_sub <- droplevels(group_IDs_sub)
    
    d_counts_sub <- d_counts[, cond[j, ]]
    
    # test for differentially abundant (DA) clusters
    runtime_DA <- system.time(
      res_DA <- testDA(d_counts_sub, group_IDs_sub, plot = TRUE, 
                       path = paste0("../../../plots/diffcyt/AML_spike_in/", thresholds[th], "/DA/", cond_names[j]))
    )
    
    # show results
    rowData(res_DA)
    
    # sort to show top (most highly significant) clusters first
    # (note: using raw p-values since BH adjustment does not work well in this simulated data set)
    res_DA_sorted <- rowData(res_DA)[order(rowData(res_DA)$P.Value), ]
    
    print(head(res_DA_sorted, 10))
    #View(res_DA_sorted)
    
    
    # --------
    # analysis
    # --------
    
    # number of significant DA clusters
    # (note: using raw p-values since BH adjustment does not work well in this simulated data set)
    print(table(res_DA_sorted$P.Value < 0.05))
    
    # save output objects
    out_DA[[th]][[j]] <- res_DA
    out_DA_sorted[[th]][[j]] <- res_DA_sorted
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


