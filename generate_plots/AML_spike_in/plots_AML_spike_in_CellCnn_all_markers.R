##########################################################################################
# Generate plots
# 
# - method: CellCnn-all-markers
# - data set: AML-spike-in
# 
# Lukas Weber, July 2017
##########################################################################################


# also add t-SNE plots (at cell level, with highlighting for true spike-in cells)


library(iCOBRA)
library(ggplot2)


# load saved results
load("../../../RData/AML_spike_in/outputs_AML_spike_in_CellCnn_all_markers.RData")




######################
# Loop over thresholds
######################

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")




for (th in 1:length(thresholds)) {
  
  #############################
  # Plots: CellCnn test results
  #############################
  
  
  # --------------------------------
  # ROC curves: CellCnn test results
  # --------------------------------
  
  # generate ROC curves using iCOBRA package
  
  for (j in 1:length(cond_names)) {
    
    # create 'COBRAData' object
    data <- out_CellCnn_all_markers[[th]][[j]]
    cobradata <- COBRAData(score = data.frame(CellCnn_all_markers = data[, "scores"]), 
                           truth = data.frame(spikein = data[, "spikein"]))
    
    # calculate performance scores
    cobraperf <- calculate_performance(cobradata, aspects = "roc", binary_truth = "spikein")
    
    # prepare object for plotting
    cobraplot <- prepare_data_for_plot(cobraperf)
    
    # create plot
    p <- plot_roc(cobraplot, linewidth = 0.75)
    
    p + 
      geom_vline(xintercept = c(0.01, 0.05, 0.1), lty = 2) + 
      coord_fixed() + 
      xlab("False positive rate") + 
      ylab("True positive rate") + 
      ggtitle(paste0("ROC: Differential test results: AML-spike-in, ", cond_names[j], ", ", thresholds[th])) + 
      theme_bw() + 
      theme(strip.text.x = element_blank())
    
    path <- paste0("../../../plots/AML_spike_in/CellCnn/all_markers/", thresholds[th], "/", cond_names[j])
    filename <- file.path(path, paste0("results_CellCnn_all_markers_ROC_curve_", thresholds[th], "_", cond_names[j], ".pdf"))
    
    ggsave(filename, width = 9, height = 8)
  }
  
}



