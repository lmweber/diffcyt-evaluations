##########################################################################################
# Analysis and plots
# 
# - method: diffcyt
# - data set: AML-spike-in
# 
# Lukas Weber, July 2017
##########################################################################################


library(flowCore)
library(SummarizedExperiment)
library(ggplot2)
library(reshape2)
library(magrittr)
library(dplyr)
library(scales)
library(Rtsne)
library(ROCR)




####################
# Load saved results
####################

load("../../../RData/outputs_diffcyt_AML_spike_in_DA_limma.RData")




######################
# Loop over thresholds
######################


# lists to store clustering performance results
# note: clustering is performed once on all samples together, including multiple conditions
clustering_pr <- clustering_re <- clustering_F1 <- vector("list", length(thresholds))
names(clustering_pr) <- names(clustering_re) <- names(clustering_F1) <- thresholds


# note: thresholds are defined in previous script
for (th in 1:length(thresholds)) {
  
  
  ######################################################
  # Calculate clustering performance for spiked-in cells
  ######################################################
  
  
  # ----------------
  # get data objects
  # ----------------
  
  d_se <- d_se_thresholds[[th]]
  d_counts <- d_counts_thresholds[[th]]
  d_medians_all <- d_medians_all_thresholds[[th]]
  
  
  # ---------------------------------
  # get spike-in status for each cell
  # ---------------------------------
  
  # number of spiked-in cells per sample
  n_spikein_tbl <- lapply(is_spikein[[th]], table)
  n_spikein_tbl
  
  vals <- as.numeric(unlist(lapply(n_spikein_tbl, names)))
  n <- unname(unlist(n_spikein_tbl))
  
  # spike-in status for each cell
  spikein_rep <- as.factor(rep(vals, n))
  
  length(spikein_rep)
  
  # store spike-in status in rowData of 'd_se'
  rowData(d_se) <- cbind(rowData(d_se), data.frame(spikein = spikein_rep))
  rowData(d_se)
  
  
  # -----------------------------------------------------------
  # calculate clustering performance for best-matching clusters
  # -----------------------------------------------------------
  
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
  
  # total number of spike-in cells
  
  n_spikein <- sapply(d_split, function(d) sum(d$spikein == 1))
  n_spikein
  
  # calculate precision, recall, F1 score
  
  pr <- n_matched_samp_correct / n_matched_samp
  re <- n_matched_samp_correct / n_spikein
  F1 <- 2 * (pr * re) / (pr + re)
  
  # store results
  clustering_pr[[th]] <- pr
  clustering_re[[th]] <- re
  clustering_F1[[th]] <- F1
  
  
  
  
  ###############################
  # Plots: clustering performance
  ###############################
  
  
  # -------------------------------------------------------------
  # Barplots: clustering performance: precision, recall, F1 score
  # -------------------------------------------------------------
  
  # note: thresholds and conditions are defined in previous script
  
  for (j in 1:length(cond_names)) {
    
    d_plot <- data.frame(
      precision = clustering_pr[[th]][group_IDs == cond_names[j]], 
      recall    = clustering_re[[th]][group_IDs == cond_names[j]], 
      F1_score  = clustering_F1[[th]][group_IDs == cond_names[j]]
    )
    
    d_plot$patient <- block_IDs[group_IDs == cond_names[j]]
    
    d_plot <- melt(d_plot, id.vars = "patient")
    
    ggplot(d_plot, aes(x = patient, y = value, group = variable, fill = variable)) + 
      geom_bar(stat = "identity", position = "dodge") + 
      scale_fill_manual(values = c("lightblue", "darkblue", "orange")) + 
      ylim(c(0, 1)) + 
      ggtitle(paste0("Clustering performance: AML-spike-in, ", cond_names[j], ", threshold ", thresholds[th])) + 
      theme_bw() + 
      theme(axis.title.y = element_blank(), 
            legend.title = element_blank())
    
    path <- paste0("../../../plots/diffcyt/AML_spike_in/", thresholds[th], "/clustering/", cond_names[j])
    filename <- file.path(path, paste0("clus_perf_AML_spike_in_", thresholds[th], "_", cond_names[j], ".pdf"))
    
    ggsave(filename, width = 6, height = 5)
  }
  
  
  # -----------------------------------------------------
  # MST plots: proportion true spike-in cells per cluster
  # -----------------------------------------------------
  
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
    filename <- file.path(path, "MST_prop_true_spikein.pdf")
    
    ggsave(filename, width = 9, height = 9)
  }
  
  
  # -------------------------------------------------------
  # t-SNE plots: proportion true spike-in cells per cluster
  # -------------------------------------------------------
  
  for (j in 1:length(cond_names)) {
    
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
      d_prop
    
    d_prop <- as.data.frame(d_prop)
    
    # fill in any missing clusters (zero cells)
    if (nrow(d_prop) < nlevels(rowData(d_se)$cluster)) {
      ix_missing <- which(!(levels(rowData(d_se)$cluster) %in% d_prop$cluster))
      d_prop_tmp <- data.frame(factor(ix_missing, levels = levels(rowData(d_se)$cluster)), 0)
      colnames(d_prop_tmp) <- colnames(d_prop)
      rownames(d_prop_tmp) <- ix_missing
      d_prop <- rbind(d_prop, d_prop_tmp)
      # re-order rows
      d_prop <- d_prop[order(d_prop$cluster), ]
      rownames(d_prop) <- d_prop$cluster
    }
    
    # number of cells
    n_cells <- rowData(d_counts_sub)$n_cells
    
    if (!(nrow(d_prop) == length(n_cells))) warning("number of clusters does not match")
    if (!(nrow(d_medians_all) == nrow(d_prop))) warning("number of clusters does not match")
    
    d_plot <- cbind(d_prop, n_cells)
    
    # run t-SNE
    
    d_tsne <- assay(d_medians_all)[, colData(d_medians_all)$is_clustering_col]
    d_tsne <- as.matrix(d_tsne)
    
    # remove any duplicate rows (required by Rtsne)
    dups <- duplicated(d_tsne)
    d_tsne <- d_tsne[!dups, ]
    
    # also remove duplicated rows from plotting data
    d_plot <- d_plot[!dups, ]
    
    # run Rtsne
    # (note: initial PCA step not required, since we do not have too many dimensions)
    set.seed(123)
    out_tsne <- Rtsne(d_tsne, pca = FALSE, verbose = TRUE)
    
    tsne_coords <- as.data.frame(out_tsne$Y)
    colnames(tsne_coords) <- c("tSNE_1", "tSNE_2")
    
    d_plot <- cbind(d_plot, tsne_coords)
    
    # plot
    ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, size = n_cells, color = prop_spikein)) + 
      geom_point(alpha = 0.75) + 
      scale_size_continuous(range = c(0.25, 4)) + 
      scale_color_gradient(low = "gray70", high = "orange") + 
      coord_fixed() + 
      ggtitle("t-SNE: Proportion true spike-in cells per cluster") + 
      theme_bw()
    
    path <- paste0("../../../plots/diffcyt/AML_spike_in/", thresholds[th], "/clustering/", cond_names[j])
    filename <- file.path(path, "tSNE_prop_true_spikein.pdf")
    
    ggsave(filename, width = 9, height = 9)
  }
  
  
  
  
  #################################################
  # Plots: differential abundance (DA) test results
  #################################################
  
  
  # --------------------------
  # MST plots: DA test results
  # --------------------------
  
  for (j in 1:length(cond_names)) {
    
    # get results object
    
    res_DA <- out_DA[[th]][[j]]
    
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
    
    
    # generate MST plot highlighting true spike-in cells on same plot
    
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
    
    path <- paste0("../../../plots/diffcyt/AML_spike_in/", thresholds[th], "/DA/", cond_names[j])
    filename <- file.path(path, "MST_results_DA.pdf")
    
    ggsave(filename, width = 9, height = 9)
  }
  
  
  # ----------------------
  # t-SNE: DA test results
  # ----------------------
  
  for (j in 1:length(cond_names)) {
    
    # get results object
    
    res_DA <- out_DA[[th]][[j]]
    
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
    
    
    # number of cells
    n_cells <- rowData(d_counts_sub)$n_cells
    
    if (!(nrow(d_true) == length(n_cells))) warning("number of clusters does not match")
    if (!(nrow(d_medians_all) == nrow(d_true))) warning("number of clusters does not match")
    
    d_plot <- cbind(d_true, n_cells)
    
    
    # generate t-SNE plot highlighting true spike-in cells on same plot
    
    # run t-SNE
    
    d_tsne <- assay(d_medians_all)[, colData(d_medians_all)$is_clustering_col]
    d_tsne <- as.matrix(d_tsne)
    
    # remove any duplicate rows (required by Rtsne)
    dups <- duplicated(d_tsne)
    d_tsne <- d_tsne[!dups, ]
    
    # also remove duplicated rows from plotting data
    d_plot <- d_plot[!dups, ]
    
    # run Rtsne
    # (note: initial PCA step not required, since we do not have too many dimensions)
    set.seed(123)
    out_tsne <- Rtsne(d_tsne, pca = FALSE, verbose = TRUE)
    
    tsne_coords <- as.data.frame(out_tsne$Y)
    colnames(tsne_coords) <- c("tSNE_1", "tSNE_2")
    
    d_plot <- cbind(d_plot, tsne_coords)
    
    
    nroot_trans <- function() {
      trans_new("nroot", function(x) x^(1/10), function(x) x^10)
    }
    
    #if (pvalue_type == "adjusted") p_vals_DA <- rowData(res_DA)$adj.P.Val
    #if (pvalue_type == "raw") p_vals_DA <- rowData(res_DA)$P.Value
    p_vals_DA <- rowData(res_DA)$adj.P.Val
    
    names(p_vals_DA) <- rowData(res_DA)$cluster
    
    stopifnot(length(p_vals_DA) == nrow(d_plot))
    
    d_plot <- cbind(d_plot, p_vals = p_vals_DA[as.character(d_plot$cluster)])
    
    min_val <- min(p_vals_DA, na.rm = TRUE)
    max_val <- max(p_vals_DA, na.rm = TRUE) - 0.3  # slightly reduce max value for legend display
                                                   # due to ggplot2 bug (see below)
    
    
    ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, size = n_cells, color = p_vals_DA)) + 
      # layer multiple geom_points to outline true spike-in cells
      # (a bit hacky)
      geom_point(aes(stroke = 1.5 * spikein), shape = 20, color = "black", fill = "white") + 
      geom_point(aes(stroke = 0), shape = 20, color = "white") + 
      geom_point(aes(stroke = 0), shape = 20, alpha = 0.6) + 
      scale_size_continuous(range = c(1, 6)) + 
      scale_color_gradient(low = "red", high = "gray60", trans = nroot_trans(), 
                           breaks = c(min_val, 0.05, max_val), 
                           labels = c(round(min_val), 0.05, round(max_val))) + 
      #guide = guide_colorbar(title.vjust = 0.5)) +  # bug in ggplot: doesn't work
      coord_fixed() + 
      ggtitle("t-SNE: Differential abundance (DA) test results") + 
      theme_bw()
    
    path <- paste0("../../../plots/diffcyt/AML_spike_in/", thresholds[th], "/DA/", cond_names[j])
    filename <- file.path(path, "tSNE_results_DA.pdf")
    
    ggsave(filename, width = 9, height = 9)
  }
  
  
  # ---------------------------
  # ROC curves: DA test results
  # ---------------------------
  
  # note: calculate ROC curves at cell level using cluster-level p-values (i.e. the
  # cluster-level p-value is used for all cells in that cluster)
  
  for (j in 1:length(cond_names)) {
    
    # get results object
    res_DA <- out_DA[[th]][[j]]
    
    # identify clusters containing true spike-in cells
    ix_keep <- group_IDs %in% c("healthy", cond_names[j])
    group_IDs_sub <- group_IDs[ix_keep]
    group_IDs_sub <- droplevels(group_IDs_sub)
    d_se_sub <- d_se[rowData(d_se)$group %in% group_IDs_sub, ]
    
    # get cell-level DA test results
    
    # convert to factor
    rowData(res_DA)$cluster <- factor(rowData(res_DA)$cluster, levels = levels(rowData(d_se_sub)$cluster))
    
    # match cluster-level p-values to individual cells
    
    # use raw p-values (adjusted p-values do not work well for very rare populations)
    p_vals_clusters <- rowData(res_DA)$P.Value
    p_vals_cells <- p_vals_clusters[match(rowData(d_se_sub)$cluster, rowData(res_DA)$cluster)]
    
    rowData(d_se_sub)$p_vals <- p_vals_cells
    
    d_roc <- as.data.frame(rowData(d_se_sub)[, c("cluster", "spikein", "p_vals")])
    
    
    # calculate ROC curve values
    pred <- prediction(1 - d_roc$p_vals, d_roc$spikein)
    perf <- performance(pred, "tpr", "fpr")
    
    FPR <- perf@x.values[[1]]
    TPR <- perf@y.values[[1]]
    x_label <- perf@x.name
    y_label <- perf@y.name
    
    
    # FPR and TPR when using actual cutoff on p-values
    # (note: careful to treat NAs correctly)
    fn_actual_fpr <- function(cutoff) {
      sum(d_roc$p_vals < cutoff & d_roc$spikein == 0, na.rm = TRUE) / 
        (sum(d_roc$spikein == 0) - sum(is.na(d_roc$p_vals) & d_roc$spikein == 0))
    }
    fn_actual_tpr <- function(cutoff) {
      sum(d_roc$p_vals < cutoff & d_roc$spikein == 1, na.rm = TRUE) / 
        (sum(d_roc$spikein == 1) - sum(is.na(d_roc$p_vals) & d_roc$spikein == 1))
    }
    
    cutoffs <- c(0.01, 0.05, 0.1)
    
    actual_fpr <- sapply(cutoffs, fn_actual_fpr)
    actual_tpr <- sapply(cutoffs, fn_actual_tpr)
    
    d_actual <- as.data.frame(cbind(actual_fpr, actual_tpr))
    rownames(d_actual) <- cutoffs
    colnames(d_actual) <- c("FPR", "TPR")
    
    
    # better plot axes (for 5% spike-in threshold only)
    if(th == 1) {
      FPR_max <- 0.2
      TPR_min <- 0.8
    } else {
      FPR_max <- 1
      TPR_min <- 0
    }
    
    # plotting data frame
    d_plot <- data.frame(FPR, TPR)
    
    if (th == 1) {
      # subset
      d_plot <- d_plot[d_plot$FPR <= FPR_max & d_plot$TPR >= TPR_min, ]
      # include min and max values (for plot)
      d_plot <- rbind(d_plot, c(FPR_max, max(d_plot$TPR)), c(min(d_plot$FPR), TPR_min))
    }
    
    # plot
    ggplot(d_plot, aes(x = FPR, y = TPR, lty = "diffcyt-DA-limma")) + 
      geom_line(color = "blue") + 
      geom_vline(xintercept = c(0.01, 0.05, 0.1), color = "red", lty = 2) + 
      geom_point(data = d_actual, shape = 1, size = 3, stroke = 1, col = "orangered") + 
      xlim(0, FPR_max) + 
      ylim(TPR_min, 1) + 
      xlab(x_label) + 
      ylab(y_label) + 
      coord_fixed() + 
      ggtitle("ROC: Differential abundance (DA) test results: AML-spike-in") + 
      theme_bw() + 
      theme(legend.title = element_blank())
    
    path <- paste0("../../../plots/diffcyt/AML_spike_in/", thresholds[th], "/DA/", cond_names[j])
    filename <- file.path(path, "ROC_curves_DA.pdf")
    
    ggsave(filename, width = 9, height = 8)
  }
  
}



