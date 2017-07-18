##########################################################################################
# Generate plots
# 
# - method: CellCnn-lineage-markers
# - data set: AML-spike-in
# 
# Lukas Weber, July 2017
##########################################################################################


# also add t-SNE plots (at cell level, with highlighting for true spike-in cells)


library(ggplot2)
library(ROCR)




####################
# Load saved results
####################

load("../../../RData/AML_spike_in/outputs_AML_spike_in_CellCnn_lineage_markers.RData")




######################
# Loop over thresholds
######################


# note: thresholds are defined in previous script
for (th in 1:length(thresholds)) {
  
  ###################################
  # get filenames and spike-in status
  ###################################
  
  # get filenames
  files_load <- files_load_thresholds[[th]]
  
  # get spike-in status for each cell
  is_spikein <- is_spikein_thresholds[[th]]
  
  all(sample_IDs == names(is_spikein))
  
  
  
  
  #############################
  # Plots: CellCnn test results
  #############################
  
  
  # --------------------------------
  # ROC curves: CellCnn test results
  # --------------------------------
  
  # note: calculate ROC curves at cell level using cluster-level p-values (i.e. the
  # cluster-level p-value is used for all cells in that cluster)
  
  for (j in 1:length(cond_names)) {
    
    # select samples for this condition
    ix_keep_cnd <- group_IDs == cond_names[j]
    
    path_out <- paste0("../../../CellCnn_files/AML_spike_in/lineage_markers/out_CellCnn/", thresholds[th], "/", cond_names[j], "/selected_cells")
    
    # skip if no files exist (CellCnn did not run correctly)
    if (length(list.files(path_out)) == 0) next
    
    # filenames for this condition
    files_cnd <- paste0(path_out, "/", gsub("\\.fcs$", "", basename(files_load[ix_keep_cnd])), "_transf_selected_cells.csv")
    
    # spike-in status for cells in this condition
    is_spikein_cnd <- unlist(is_spikein[ix_keep_cnd])
    
    # get significant cells in top filter for this condition
    
    filter_continuous_cnd <- vector("list", length(files_cnd))
    
    for (f in 1:length(files_cnd)) {
      d <- try(read.csv(files_cnd[f]), silent = TRUE)
      
      if (!(class(d) == "try-error")) {
        # note: if there are multiple filters for one sample, combine them using 'rowSums'
        # (filters are stored in odd-numbered columns of the .csv file)
        ix <- seq(1, ncol(d), by = 2)
        filt_sum <- rowSums(d[, ix, drop = FALSE])
        filter_continuous_cnd[[f]] <- filt_sum
        
      } else {
        # if .csv file is empty, fill with zeros instead
        filter_continuous_cnd[[f]] <- rep(0, length(is_spikein[ix_keep_cnd][[f]]))
      }
    }
    
    filter_continuous_cnd <- unlist(filter_continuous_cnd)
    length(is_spikein_cnd) == length(filter_continuous_cnd)
    
    # skip if all empty
    if (all(filter_continuous_cnd == 0)) next
    
    # calculate ROC curve values
    
    d_roc <- data.frame(spikein = is_spikein_cnd, filter = filter_continuous_cnd)
    #View(d_roc)
    
    pred <- prediction(d_roc$filter, d_roc$spikein)
    perf <- performance(pred, "tpr", "fpr")
    
    FPR <- perf@x.values[[1]]
    TPR <- perf@y.values[[1]]
    x_label <- perf@x.name
    y_label <- perf@y.name
    
    # remove point (1, 1) from ROC curve, since this gives misleading diagonal lines
    ix_ones <- which(FPR == 1 & TPR == 1)
    FPR <- FPR[-ix_ones]
    TPR <- TPR[-ix_ones]
    
    # plotting data frame
    d_plot <- data.frame(FPR, TPR)
    
    # plot
    ggplot(d_plot, aes(x = FPR, y = TPR, lty = "CellCnn-lineage-markers")) + 
      geom_line(color = "forestgreen") + 
      geom_vline(xintercept = c(0.01, 0.05, 0.1), color = "red", lty = 2) + 
      xlim(0, 1) + 
      ylim(0, 1) + 
      xlab(x_label) + 
      ylab(y_label) + 
      coord_fixed() + 
      ggtitle(paste0("ROC: Differential test results: AML-spike-in, ", cond_names[j], ", ", thresholds[th])) + 
      theme_bw() + 
      theme(legend.title = element_blank())
    
    path <- paste0("../../../plots/AML_spike_in/CellCnn/lineage_markers/", thresholds[th], "/", cond_names[j])
    filename <- file.path(path, paste0("results_CellCnn_lineage_markers_ROC_curve_", thresholds[th], "_", cond_names[j], ".pdf"))
    
    ggsave(filename, width = 9, height = 8)
  }
  
}



