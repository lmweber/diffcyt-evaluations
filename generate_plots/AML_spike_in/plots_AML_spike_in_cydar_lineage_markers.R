##########################################################################################
# Generate plots
# 
# - method: cydar-lineage-markers
# - data set: AML-spike-in
# 
# Lukas Weber, July 2017
##########################################################################################


library(iCOBRA)
library(ggplot2)


# load saved results
load("../../../RData/AML_spike_in/outputs_AML_spike_in_cydar_lineage_markers.RData")




######################
# Loop over thresholds
######################

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")




for (th in 1:length(thresholds)) {
  
  ###########################
  # Plots: cydar test results
  ###########################
  
  
  # ------------------------------
  # ROC curves: cydar test results
  # ------------------------------
  
  # generate ROC curves using iCOBRA package
  
  for (j in 1:length(cond_names)) {
    
    # create 'COBRAData' object
    data <- out_cydar_lineage_markers[[th]][[j]]
    cobradata <- COBRAData(pval = data.frame(cydar_lineage_markers = data[, "p_vals"]), 
                           padj = data.frame(cydar_lineage_markers = data[, "q_vals"]), 
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
    
    path <- paste0("../../../plots/AML_spike_in/cydar/lineage_markers/", thresholds[th], "/", cond_names[j])
    filename <- file.path(path, paste0("results_cydar_lineage_markers_ROC_curve_", thresholds[th], "_", cond_names[j], ".pdf"))
    
    ggsave(filename, width = 9, height = 8)
  }
  
}



