##########################################################################################
# Generate plots
# 
# - method: Citrus-all-markers
# - data set: AML-spike-in
# 
# Lukas Weber, July 2017
##########################################################################################


library(ggplot2)
library(ROCR)




####################
# Load saved results
####################

load("../../../RData/AML_spike_in/outputs_AML_spike_in_Citrus_all_markers.RData")




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
  
  
  
  
  ############################
  # Plots: Citrus test results
  ############################
  
  
  # -------------------------------
  # ROC curves: Citrus test results
  # -------------------------------
  
  for (j in 1:length(cond_names)) {
    
    # get objects
    citrus_cells <- out_Citrus_cells[[th]][[j]]
    
    # spike-in status for cells in this condition
    ix_keep_cnd <- group_IDs == cond_names[j]
    is_spikein_cnd <- unlist(is_spikein[ix_keep_cnd])
    
    length(is_spikein_cnd)
    
    # Citrus results for cells in this condition
    ix_keep_rep <- rep(ix_keep_cnd, n_cells)
    citrus_cells_cnd <- citrus_cells[ix_keep_rep]
    
    length(citrus_cells_cnd)
    
    stopifnot(length(is_spikein_cnd) == length(citrus_cells_cnd))
    
    # calculate ROC curve values
    d_roc <- data.frame(spikein = is_spikein_cnd, detected = citrus_cells_cnd)
    #View(d_roc)
    
    # set NA values to 0 to allow ROC curve to be calculated
    d_roc$detected[is.na(d_roc$detected)] <- 0
    
    pred <- prediction(d_roc$detected, d_roc$spikein)
    perf <- performance(pred, "tpr", "fpr")
    
    FPR <- perf@x.values[[1]]
    TPR <- perf@y.values[[1]]
    x_label <- perf@x.name
    y_label <- perf@y.name
    
    # plotting data frame
    d_plot <- data.frame(FPR, TPR)
    
    # plot
    ggplot(d_plot, aes(x = FPR, y = TPR, lty = "Citrus-all-markers")) + 
      geom_line(color = "magenta") + 
      geom_vline(xintercept = c(0.01, 0.05, 0.1), color = "red", lty = 2) + 
      xlim(0, 1) + 
      ylim(0, 1) + 
      xlab(x_label) + 
      ylab(y_label) + 
      coord_fixed() + 
      ggtitle(paste0("ROC: Differential test results: AML-spike-in, ", cond_names[j], ", ", thresholds[th])) + 
      theme_bw() + 
      theme(legend.title = element_blank())
    
    path <- paste0("../../../plots/AML_spike_in/Citrus/all_markers/", thresholds[th], "/", cond_names[j])
    filename <- file.path(path, paste0("results_Citrus_all_markers_ROC_curve_", thresholds[th], "_", cond_names[j], ".pdf"))
    
    ggsave(filename, width = 9, height = 8)
  }
  
}



