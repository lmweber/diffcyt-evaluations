##########################################################################################
# Analysis and plots
# 
# - method: CellCnn
# - data set: AML-spike-in
# 
# Lukas Weber, July 2017
##########################################################################################


library(ggplot2)
library(ROCR)




####################
# Load saved results
####################

load("../../../RData/outputs_CellCnn_AML_spike_in_all_markers.RData")




######################
# Loop over thresholds
######################


# note: thresholds are defined in previous script
for (th in 1:length(thresholds)) {
  
  ###################################
  # get filenames and spike-in status
  ###################################
  
  # -------------
  # get filenames
  # -------------
  
  files_load <- files_load_thresholds[[th]]
  
  
  # ---------------------------------
  # get spike-in status for each cell
  # ---------------------------------
  
  all(sample_IDs == names(is_spikein[[th]]))
  
  spikein <- is_spikein[[th]]
  
  
  
  
  #############################
  # Plots: CellCnn test results
  #############################
  
  
  # also try a t-SNE plot (at cell level, with highlighting for true spike-in cells)
  
  
  
  # --------------------------------
  # ROC curves: CellCnn test results
  # --------------------------------
  
  # note: calculate ROC curves at cell level using cluster-level p-values (i.e. the
  # cluster-level p-value is used for all cells in that cluster)
  
  for (j in 1:length(cond_names)) {
    
    # select samples for this condition
    ix_keep_cnd <- group_IDs == cond_names[j]
    
    path_out <- paste0("../../../CellCnn_files/out_CellCnn/AML_spike_in/", thresholds[th], "/", cond_names[j], "/selected_cells")
    
    
    # skip if no files exist (CellCnn did not run correctly)
    if (length(list.files(path_out)) == 0) next
    
    
    # filenames for this condition
    files_cnd <- paste0(path_out, "/", gsub("\\.fcs$", "", basename(files_load[ix_keep_cnd])), "_transf_selected_cells.csv")
    
    # spike-in status for cells in this condition
    spikein_cnd <- spikein[ix_keep_cnd]
    spikein_cnd <- unlist(spikein_cnd)
    
    # significant cells in top filter for this condition
    
    # assume only one filter for now (to do: allow multiple filters)
    filter_continuous_cnd <- lapply(files_cnd, function(f) {
      read.csv(f)[, 1]  ## column index 1 contains continuous values for the top filter
    })
    
    filter_continuous_cnd <- unlist(filter_continuous_cnd)
    
    length(spikein_cnd) == length(filter_continuous_cnd)
    
    # calculate ROC curve values
    d_roc <- data.frame(spikein = spikein_cnd, filter = filter_continuous_cnd)
    #View(d_roc)
    
    pred <- prediction(d_roc$filter, d_roc$spikein)
    perf <- performance(pred, "tpr", "fpr")
    
    # to do: replace last row value=c(1.0, 1.0) with value=c(FPR=1.0, TPR=TPR_max)
    # this will remove the straight diagonal line, which is misleading
    
    
    FPR <- perf@x.values[[1]]
    TPR <- perf@y.values[[1]]
    x_label <- perf@x.name
    y_label <- perf@y.name
    
    
    # note: CellCnn: no p-values, so cannot calculate FPR and TPR using actual p-value cutoffs
    
    
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
    ggplot(d_plot, aes(x = FPR, y = TPR, lty = "CellCnn-all-markers")) + 
      geom_line(color = "forestgreen") + 
      geom_vline(xintercept = c(0.01, 0.05, 0.1), color = "red", lty = 2) + 
      xlim(0, FPR_max) + 
      ylim(TPR_min, 1) + 
      xlab(x_label) + 
      ylab(y_label) + 
      coord_fixed() + 
      ggtitle(paste0("ROC: Differential test results: AML-spike-in, ", cond_names[j], ", ", thresholds[th])) + 
      theme_bw() + 
      theme(legend.title = element_blank())
    
    path <- paste0("../../../plots/CellCnn/AML_spike_in/", thresholds[th], "/", cond_names[j])
    filename <- file.path(path, "ROC_curves_CellCnn.pdf")
    
    ggsave(filename, width = 9, height = 8)
  }
  
}



