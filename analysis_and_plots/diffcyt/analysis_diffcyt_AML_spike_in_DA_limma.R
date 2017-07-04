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
  # barplots: clustering performance: precision, recall, F1 score
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
}



