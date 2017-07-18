##########################################################################################
# Generate plots
# 
# - method: cydar
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

load("../../../RData/outputs_cydar_AML_spike_in.RData")




######################
# Loop over thresholds
######################


# note: thresholds are defined in previous script
for (th in 1:length(thresholds)) {
  
  #########################################
  # get spike-in status and number of cells
  #########################################
  
  # get spike-in status for each cell
  all(sample_IDs == names(is_spikein_thresholds[[th]]))
  is_spikein <- is_spikein_thresholds[[th]]
  
  # get number of cells per sample (including spike-in cells)
  n_cells <- n_cells_thresholds[[th]]
  
  
  
  
  ###########################
  # Plots: cydar test results
  ###########################
  
  
  # also try a t-SNE plot (at cell level, with highlighting for true spike-in cells)
  
  
  
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
    
    # get spike-in status for cells in this condition
    ix_keep_cnd <- group_IDs == cond_names[j]
    ix_keep_cnd <- rep(ix_keep_cnd, n_cells)
    
    is_spikein_cnd <- unlist(is_spikein)[ix_keep_cnd]
    
    length(is_spikein_cnd)
    
    
    # get smallest q-value for each cell, across all hyperspheres
    cells <- cellAssignments(cydar_data)
    
    # 'unpack' indices (e.g. '142183 -142188' means all values from 142183 to 142188)
    cells <- unpackIndices(cells)
    
    length(cells)
    length(vals)
    stopifnot(length(vals) == length(cells))
    
    cells_rep <- unlist(cells)
    stopifnot(length(cells_rep) == sum(sapply(cells, length)))
    names(cells_rep) <- rep(vals, sapply(cells, length))
    
    cells_rep_split <- split(cells_rep, cells_rep)
    
    cells_vals <- sapply(cells_rep_split, function(s) min(as.numeric(names(s))))
    
    vals_all <- rep(NA, sum(n_cells))
    names(vals_all) <- 1:sum(n_cells)
    vals_all[names(cells_vals)] <- cells_vals
    
    length(vals_all)
    
    # this condition only
    vals_all_cnd <- vals_all[ix_keep_cnd]
    length(vals_all_cnd)
    
    
    # calculate ROC curve values
    stopifnot(length(vals_all_cnd) == length(is_spikein_cnd))
    
    d_roc <- data.frame(spikein = is_spikein_cnd, vals = vals_all_cnd)
    #View(d_roc)
    
    pred <- prediction(1 - d_roc$vals, d_roc$spikein)
    perf <- performance(pred, "tpr", "fpr")
    
    #plot(perf)
    
    FPR <- perf@x.values[[1]]
    TPR <- perf@y.values[[1]]
    x_label <- perf@x.name
    y_label <- perf@y.name
    
    # plotting data frame
    d_plot <- data.frame(FPR, TPR)
    
    
    # plot
    ggplot(d_plot, aes(x = FPR, y = TPR, lty = "cydar")) + 
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
    
    path <- paste0("../../../plots/cydar/AML_spike_in/", thresholds[th], "/", cond_names[j])
    filename <- file.path(path, "ROC_curves_cydar.pdf")
    
    ggsave(filename, width = 9, height = 8)
  }
  
}



