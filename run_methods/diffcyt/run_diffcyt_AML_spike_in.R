##########################################################################################
# Script to run methods
# 
# - method: diffcyt
# - data set: AML-spike-in
# 
# Lukas Weber, July 2017
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(magrittr)
library(limma)
library(scales)




################################
# Loop to run for each threshold
################################

# spike-in thresholds (must match filenames)
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")


# lists to store clustering performance results
clustering_pr <- clustering_re <- clustering_F1 <- vector("list", length(thresholds))
names(clustering_pr) <- names(clustering_re) <- names(clustering_F1) <- thresholds

# lists to store differential abundance (DA) test results
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
  
  
  # ------------------------------------------------------
  # remove spike-in indicator columns and store separately
  # ------------------------------------------------------
  
  # 'AML-spike-in' data set includes indicator columns for spike-in cells in conditions
  # 'CN' and 'CBF'; these need to be removed and stored separately
  is_spikein <- vector("list", length(sample_IDs))
  names(is_spikein) <- sample_IDs
  
  for (i in 1:length(is_spikein)) {
    exprs_i <- exprs(d_input[[i]])
    if (group_IDs[i] == "healthy") {
      is_spikein[[i]] <- rep(0, nrow(exprs_i))
    } else {
      is_spikein[[i]] <- exprs_i[, "spikein"]
      # remove column from expression matrix
      exprs(d_input[[i]]) <- exprs_i[, -match("spikein", colnames(exprs_i))]
    }
  }
  
  
  
  
  ############################################################
  # Run initial steps in 'diffcyt' pipeline (up to clustering)
  ############################################################
  
  # prepare data into required format
  # (note: using lineage markers for clustering, and functional markers for DE testing)
  d_se <- prepareData(d_input, sample_IDs, cols_markers, cols_lineage, cols_func)
  
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
  
  
  
  
  #######################################################
  # Evaluate clustering performance (for spiked-in cells)
  #######################################################
  
  # -----------------
  # construct vectors
  # -----------------
  
  # number of spiked-in cells per sample
  n_spikein_tbl <- lapply(is_spikein, table)
  n_spikein_tbl
  
  vals <- as.numeric(unlist(lapply(n_spikein_tbl, names)))
  n <- unname(unlist(n_spikein_tbl))
  
  # spike-in status for each cell
  spikein_rep <- as.factor(rep(vals, n))
  
  length(spikein_rep)
  
  # group_ID for each cell
  n_tot <- sapply(d_input, nrow)
  group_IDs_rep <- rep(group_IDs, n_tot)
  
  length(group_IDs_rep)
  
  # store spike-in status and group_IDs in rowData of 'd_se'
  rowData(d_se) <- cbind(rowData(d_se), 
                         data.frame(spikein = spikein_rep), 
                         data.frame(group = group_IDs_rep))
  rowData(d_se)
  
  
  # ------------------------------------------------
  # clustering performance for best-matching cluster
  # ------------------------------------------------
  
  # condition names
  cond_names <- c("CN", "CBF")
  
  th_pr <- th_re <- th_F1 <- vector("list", length(cond_names))
  names(th_pr) <- names(th_re) <- names(th_F1) <- cond_names
  
  # loop over conditions
  for (j in 1:length(cond_names)) {
    
    # find best-matching cluster label for each sample
    
    d_split <- split(rowData(d_se), rowData(d_se)$sample)
    
    labels_matched <- sapply(d_split, function(d) unname(which.max(table(d[d$spikein == 1, ]$cluster))))
    
    labels_matched[group_IDs == "healthy"] <- NA
    labels_matched
    
    # number of matching cells in best-matching cluster
    
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
    
    n_spikein <- sapply(d_split, function(d) sum(d$spikein == 1))
    n_spikein
    
    # calculate precision, recall, F1 score
    
    pr <- n_matched_samp_correct / n_matched_samp
    re <- n_matched_samp_correct / n_spikein
    F1 <- 2 * (pr * re) / (pr + re)
    
    th_pr[[j]] <- pr
    th_re[[j]] <- re
    th_F1[[j]] <- F1
  }
  
  # store results
  clustering_pr[[th]] <- th_pr
  clustering_re[[th]] <- th_re
  clustering_F1[[th]] <- th_F1
  
  
  # --------------------------------------------------
  # plot proportion of true spike-in cells per cluster
  # --------------------------------------------------
  
  # plot minimum spanning trees (MST) showing proportion of true spike-in cells per cluster
  
  for (j in 1:length(cond_names)) {
    
    ix_keep <- group_IDs %in% c("healthy", cond_names[j])
    
    group_IDs_sub <- group_IDs[ix_keep]
    group_IDs_sub <- droplevels(group_IDs_sub)
    
    d_se_sub <- d_se[rowData(d_se)$group %in% group_IDs_sub, ]
    
    d_counts_sub <- d_counts[, ix_keep]
    
    mst <- metadata(d_se)$MST
    mst_coords <- as.data.frame(mst$l)
    colnames(mst_coords) <- c("MST_x", "MST_y")
    
    # calculate proportion true spike-in cells per cluster
    rowData(d_se_sub) %>% 
      as.data.frame %>% 
      group_by(cluster) %>% 
      summarize(prop_spikein = mean(as.numeric(as.character(spikein)))) -> 
      d_plot
    
    d_plot <- as.data.frame(d_plot)
    
    # fill in any missing clusters (zero cells)
    if (nrow(d_plot) < nlevels(rowData(d_se)$cluster)) {
      ix_missing <- which(!(levels(rowData(d_se)$cluster) %in% d_plot$cluster))
      d_plot_tmp <- data.frame(factor(ix_missing, levels = levels(rowData(d_se)$cluster)), 0)
      colnames(d_plot_tmp) <- colnames(d_plot)
      rownames(d_plot_tmp) <- ix_missing
      d_plot <- rbind(d_plot, d_plot_tmp)
      # re-order rows
      d_plot <- d_plot[order(d_plot$cluster), ]
      rownames(d_plot) <- d_plot$cluster
    }
    
    n_cells <- rowData(d_counts_sub)$n_cells
    
    if (!(nrow(d_plot) == length(n_cells))) warning("number of clusters does not match")
    
    d_plot <- cbind(d_plot, n_cells, mst_coords)
    
    ggplot(d_plot, aes(x = MST_x, y = MST_y, size = n_cells, color = prop_spikein)) + 
      geom_point(alpha = 0.75) + 
      scale_color_gradient(low = "gray70", high = "orange") + 
      coord_fixed() + 
      ggtitle("MST: Proportion true spike-in cells per cluster") + 
      theme_bw()
    
    path <- paste0("../../../plots/diffcyt/AML_spike_in/", thresholds[th], "/clustering/", cond_names[j])
    
    ggsave(file.path(path, "MST_true_prop_spikein.pdf"), width = 9, height = 9)
  }
  
  
  
  
  ################################################
  # Test for differentially abundant (DA) clusters
  ################################################
  
  # test separately for each condition: CN vs. healthy, CBF vs. healthy
  
  out_DA[[th]] <- out_DA_sorted[[th]] <- vector("list", length(cond_names))
  names(out_DA[[th]]) <- names(out_DA_sorted[[th]]) <- cond_names
  
  
  for (j in 1:length(cond_names)) {
    
    # ------------
    # run DA tests
    # ------------
    
    path <- paste0("../../../plots/diffcyt/AML_spike_in/", thresholds[th], "/DA/", cond_names[j])
    
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
    
    
    # --------
    # analysis
    # --------
    
    # show results
    rowData(res_DA)
    
    # sort to show top (most highly significant) clusters first
    res_DA_sorted <- rowData(res_DA)[order(rowData(res_DA)$adj.P.Val), ]
    
    print(head(res_DA_sorted, 10))
    #View(as.data.frame(res_DA_sorted))
    
    # number of significant DA clusters
    print(table(res_DA_sorted$adj.P.Val < 0.05))
    
    # save output objects
    out_DA[[th]][[j]] <- res_DA
    out_DA_sorted[[th]][[j]] <- res_DA_sorted
    
    
    # ---------
    # MST plots
    # ---------
    
    # identify clusters containing true spike-in cells
    ix_keep <- group_IDs %in% c("healthy", cond_names[j])
    group_IDs_sub <- group_IDs[ix_keep]
    group_IDs_sub <- droplevels(group_IDs_sub)
    d_se_sub <- d_se[rowData(d_se)$group %in% group_IDs_sub, ]
    d_counts_sub <- d_counts[, ix_keep]
    # calculate proportion true spike-in cells per cluster
    rowData(d_se_sub) %>% 
      as.data.frame %>% 
      group_by(cluster) %>% 
      summarize(prop_spikein = mean(as.numeric(as.character(spikein)))) -> 
      d_true
    d_true <- as.data.frame(d_true)
    # fill in any missing clusters (zero cells)
    if (nrow(d_true) < nlevels(rowData(d_se)$cluster)) {
      ix_missing <- which(!(levels(rowData(d_se)$cluster) %in% d_true$cluster))
      d_true_tmp <- data.frame(factor(ix_missing, levels = levels(rowData(d_se)$cluster)), 0)
      colnames(d_true_tmp) <- colnames(d_true)
      rownames(d_true_tmp) <- ix_missing
      d_true <- rbind(d_true, d_true_tmp)
      # re-order rows
      d_true <- d_true[order(d_true$cluster), ]
      rownames(d_true) <- d_true$cluster
    }
    n_cells <- rowData(d_counts_sub)$n_cells
    if (!(nrow(d_true) == length(n_cells))) warning("number of clusters does not match")
    # identify if proportion spike-in > 0.1
    d_true$spikein <- as.numeric(d_true$prop_spikein > 0.1)
    
    
    # MST plot highlighting true spike-in cells on same plot
    
    mst <- metadata(d_se)$MST
    mst_coords <- as.data.frame(mst$l)
    
    cluster <- rowData(d_counts)$cluster
    n_cells <- rowData(d_counts)$n_cells
    
    if (!(nrow(mst_coords) == length(cluster))) {
      stop("minimum spanning tree (MST) does not have correct number of clusters")
    }
    
    d_plot <- data.frame(cluster, n_cells, mst_coords)
    colnames(d_plot) <- c("cluster", "n_cells", "MST_x", "MST_y")
    
    nroot_trans <- function() {
      trans_new("nroot", function(x) x^(1/10), function(x) x^10)
    }
    
    #if (pvalue_type == "adjusted") p_vals_DA <- rowData(res_DA)$adj.P.Val
    #if (pvalue_type == "raw") p_vals_DA <- rowData(res_DA)$P.Value
    p_vals_DA <- rowData(res_DA)$adj.P.Val
    
    names(p_vals_DA) <- rowData(res_DA)$cluster
    
    d_plot <- cbind(d_plot, p_vals = p_vals_DA[as.character(d_plot$cluster)])
    
    min_val <- min(p_vals_DA, na.rm = TRUE)
    max_val <- max(p_vals_DA, na.rm = TRUE) - 0.3  # slightly reduce max value for legend display
                                                   # due to ggplot2 bug (see below)
    
    # add indicator for clusters containing spike-in cells
    d_plot$spikein <- d_true$spikein
    
    ggplot(d_plot, aes(x = MST_x, y = MST_y, size = n_cells, color = p_vals)) + 
      # layer multiple geom_points to outline true spike-in cells
      geom_point(aes(stroke = (2 * spikein) + 0.1), color = "black") + 
      geom_point(color = "white") + 
      geom_point(alpha = 0.75) + 
      scale_color_gradient(low = "red", high = "gray60", trans = nroot_trans(), 
                           breaks = c(min_val, 0.05, max_val), 
                           labels = c(round(min_val), 0.05, round(max_val))) + 
      #guide = guide_colorbar(title.vjust = 0.5)) +  # bug in ggplot: doesn't work
      coord_fixed() + 
      ggtitle("MST: Differential abundance (DA) test results") + 
      theme_bw()
    
    ggsave(file.path(path, "MST_results_DA.pdf"), width = 9, height = 9)
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



