##########################################################################################
# Generate plots
# 
# - method: cydar-all-markers
# - data set: AML-spike-in
# 
# Lukas Weber, July 2017
##########################################################################################


library(cydar)
library(ggplot2)
library(ROCR)




####################
# Load saved results
####################

load("../../../RData/AML_spike_in/outputs_AML_spike_in_cydar_all_markers.RData")




######################
# Loop over thresholds
######################


# note: thresholds are defined in previous script
for (th in 1:length(thresholds)) {
  
  #########################################
  # get number of cells and spike-in status
  #########################################
  
  # number of cells per sample (including spike-in cells)
  n_cells <- n_cells_thresholds[[th]]
  
  # spike-in status for each cell
  is_spikein <- is_spikein_thresholds[[th]]
  
  all(sample_IDs == names(is_spikein))
  
  
  
  
  ###########################
  # Plots: cydar test results
  ###########################
  
  
  # ------------------------------
  # ROC curves: cydar test results
  # ------------------------------
  
  # Note: cydar evaluates q-values at the hypersphere level. Since hyperspheres overlap,
  # the q-values are not unique at the cell level. To evaluate performance at the cell
  # level, we first assign a unique q-value to each cell, by selecting the smallest
  # q-value for any hypersphere containing that cell.
  
  
  for (j in 1:length(cond_names)) {
    
    # get objects
    cydar_data <- out_cydar_data[[th]][[j]]
    cydar_pvals <- out_cydar_pvals[[th]][[j]]
    cydar_qvals <- out_cydar_qvals[[th]][[j]]
    
    # choose whether to use raw p-values or q-values
    #vals <- cydar_pvals
    vals <- cydar_qvals
    
    # spike-in status for cells in this condition
    ix_keep_cnd <- group_IDs == cond_names[j]
    is_spikein_cnd <- unlist(is_spikein[ix_keep_cnd])
    
    length(is_spikein_cnd)
    
    # get smallest q-value for each cell, across all hyperspheres
    
    # cell assignments
    cells <- cellAssignments(cydar_data)
    # 'unpack' indices (e.g. '142183 -142188' means all values from 142183 to 142188)
    cells <- unpackIndices(cells)
    
    length(cells)
    length(vals) == length(cells)
    
    # repeat q-values for each hypersphere
    cells_rep <- unlist(cells)
    vals_rep <- rep(vals, sapply(cells, length))
    stopifnot(length(cells_rep) == length(vals_rep))
    # split by cell indices
    vals_rep_split <- split(vals_rep, cells_rep)
    # get minimum q-value for each unique cell
    cells_vals <- sapply(vals_rep_split, function(v) min(v))
    
    # fill in NAs for missing cells
    vals_all <- rep(NA, sum(n_cells))
    names(vals_all) <- 1:length(vals_all)
    vals_all[names(cells_vals)] <- cells_vals
    
    length(vals_all)
    
    # select cells from this condition only
    vals_all_cnd <- vals_all[rep(ix_keep_cnd, n_cells)]
    length(vals_all_cnd)
    
    stopifnot(length(vals_all_cnd) == length(is_spikein_cnd))
    
    # calculate ROC curve values
    
    d_roc <- data.frame(spikein = is_spikein_cnd, vals = vals_all_cnd)
    #View(d_roc)
    
    pred <- prediction(1 - d_roc$vals, d_roc$spikein)
    perf <- performance(pred, "tpr", "fpr")
    
    FPR <- perf@x.values[[1]]
    TPR <- perf@y.values[[1]]
    x_label <- perf@x.name
    y_label <- perf@y.name
    
    # remove point (1, 1) from ROC curve, since this gives misleading diagonal lines
    # ix_ones <- which(FPR == 1 & TPR == 1)
    # FPR <- FPR[-ix_ones]
    # TPR <- TPR[-ix_ones]
    
    # plotting data frame
    d_plot <- data.frame(FPR, TPR)
    
    # plot
    ggplot(d_plot, aes(x = FPR, y = TPR, lty = "cydar-all-markers")) + 
      geom_line(color = "cyan") + 
      geom_vline(xintercept = c(0.01, 0.05, 0.1), color = "red", lty = 2) + 
      xlim(0, 1) + 
      ylim(0, 1) + 
      xlab(x_label) + 
      ylab(y_label) + 
      coord_fixed() + 
      ggtitle(paste0("ROC: Differential test results: AML-spike-in, ", cond_names[j], ", ", thresholds[th])) + 
      theme_bw() + 
      theme(legend.title = element_blank())
    
    path <- paste0("../../../plots/AML_spike_in/cydar/all_markers/", thresholds[th], "/", cond_names[j])
    filename <- file.path(path, paste0("results_cydar_all_markers_ROC_curve_", thresholds[th], "_", cond_names[j], ".pdf"))
    
    ggsave(filename, width = 9, height = 8)
  }
  
}



