##########################################################################################
# Generate plots
# 
# - method: all methods; lineage markers
# - data set: AML-spike-in
# 
# Lukas Weber, July 2017
##########################################################################################


library(iCOBRA)
library(ggplot2)


# load saved results
load("../../../RData/AML_spike_in/outputs_AML_spike_in_diffcyt_DA_limma_fixed_lineage_markers.RData")
load("../../../RData/AML_spike_in/outputs_AML_spike_in_diffcyt_DA_limma_random_lineage_markers.RData")
load("../../../RData/AML_spike_in/outputs_AML_spike_in_diffcyt_DA_edgeR_lineage_markers.RData")
load("../../../RData/AML_spike_in/outputs_AML_spike_in_diffcyt_DA_GLMM_lineage_markers.RData")
load("../../../RData/AML_spike_in/outputs_AML_spike_in_CellCnn_lineage_markers.RData")
load("../../../RData/AML_spike_in/outputs_AML_spike_in_cydar_lineage_markers.RData")
load("../../../RData/AML_spike_in/outputs_AML_spike_in_Citrus_lineage_markers.RData")




######################
# Loop over thresholds
######################

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")




for (th in 1:length(thresholds)) {
  
  #####################################
  # Plots: test results for each method
  #####################################
  
  
  for (j in 1:length(cond_names)) {
    
    # -------------------------------------
    # Pre-processing steps for iCOBRA plots
    # -------------------------------------
    
    # create 'COBRAData' object
    data_diffcyt_DA_limma_fixed <- out_diffcyt_DA_limma_fixed_lineage_markers[[th]][[j]]
    data_diffcyt_DA_limma_random <- out_diffcyt_DA_limma_random_lineage_markers[[th]][[j]]
    data_diffcyt_DA_edgeR <- out_diffcyt_DA_edgeR_lineage_markers[[th]][[j]]
    data_diffcyt_DA_GLMM <- out_diffcyt_DA_GLMM_lineage_markers[[th]][[j]]
    data_CellCnn <- out_CellCnn_lineage_markers[[th]][[j]]
    data_cydar <- out_cydar_lineage_markers[[th]][[j]]
    data_Citrus <- out_Citrus_lineage_markers[[th]][[j]]
    
    stopifnot(all.equal(data_diffcyt_DA_limma_fixed$spikein, data_diffcyt_DA_limma_random$spikein))
    stopifnot(all.equal(data_diffcyt_DA_limma_fixed$spikein, data_diffcyt_DA_edgeR$spikein))
    stopifnot(all.equal(data_diffcyt_DA_limma_fixed$spikein, data_diffcyt_DA_GLMM$spikein))
    stopifnot(all.equal(data_diffcyt_DA_limma_fixed$spikein, data_CellCnn$spikein))
    stopifnot(all.equal(data_diffcyt_DA_limma_fixed$spikein, data_cydar$spikein))
    stopifnot(all.equal(data_diffcyt_DA_limma_fixed$spikein, data_Citrus$spikein))
    
    # note: provide all available values here:
    # - 'padj' is required for threshold points on TPR-FDR curves
    # - depending on availability, plotting functions use 'score', then 'pval', then 'padj'
    cobradata <- COBRAData(pval = data.frame(diffcyt_DA_limma_fixed = data_diffcyt_DA_limma_fixed[, "p_vals"], 
                                             diffcyt_DA_limma_random = data_diffcyt_DA_limma_random[, "p_vals"], 
                                             diffcyt_DA_edgeR = data_diffcyt_DA_edgeR[, "p_vals"], 
                                             diffcyt_DA_GLMM = data_diffcyt_DA_GLMM[, "p_vals"], 
                                             cydar = data_cydar[, "p_vals"]), 
                           padj = data.frame(diffcyt_DA_limma_fixed = data_diffcyt_DA_limma_fixed[, "p_adj"], 
                                             diffcyt_DA_limma_random = data_diffcyt_DA_limma_random[, "p_adj"], 
                                             diffcyt_DA_edgeR = data_diffcyt_DA_edgeR[, "p_adj"], 
                                             diffcyt_DA_GLMM = data_diffcyt_DA_GLMM[, "p_adj"], 
                                             cydar = data_cydar[, "q_vals"]), 
                           score = data.frame(CellCnn = data_CellCnn[, "scores"], 
                                              Citrus = data_Citrus[, "scores"]), 
                           truth = data.frame(spikein = data_diffcyt_DA_limma_fixed[, "spikein"]))
    
    # calculate performance scores
    # (note: warning messages when 'padj' not available)
    cobraperf <- calculate_performance(cobradata, 
                                       binary_truth = "spikein", 
                                       aspects = c("fdrtpr", "fdrtprcurve", "roc"))
    
    # prepare plotting object
    cobraplot <- prepare_data_for_plot(cobraperf)
    
    
    # ----------
    # ROC curves
    # ----------
    
    # create plot
    p <- plot_roc(cobraplot, linewidth = 0.75)
    
    p + 
      coord_fixed() + 
      xlab("False positive rate") + 
      ylab("True positive rate") + 
      ggtitle(paste0("ROC curves: AML-spike-in, lineage markers, ", cond_names[j], ", ", thresholds[th])) + 
      theme_bw() + 
      theme(strip.text.x = element_blank())
    
    path <- paste0("../../../plots/AML_spike_in/all_methods/lineage_markers/", thresholds[th], "/", cond_names[j])
    filename <- file.path(path, paste0("results_all_methods_lineage_markers_ROC_curves_", thresholds[th], "_", cond_names[j], ".pdf"))
    
    ggsave(filename, width = 9, height = 8)
    
    
    # --------------
    # TPR-FDR curves
    # --------------
    
    # create plot
    p <- plot_fdrtprcurve(cobraplot, linewidth = 0.75, pointsize = 4)
    
    p + 
      coord_fixed() + 
      xlab("False discovery rate") + 
      ylab("True positive rate") + 
      ggtitle(paste0("TPR-FDR curves: AML-spike-in, lineage markers, ", cond_names[j], ", ", thresholds[th])) + 
      theme_bw() + 
      theme(strip.text.x = element_blank())
    
    path <- paste0("../../../plots/AML_spike_in/all_methods/lineage_markers/", thresholds[th], "/", cond_names[j])
    filename <- file.path(path, paste0("results_all_methods_lineage_markers_TPR_FDR_curves_", thresholds[th], "_", cond_names[j], ".pdf"))
    
    ggsave(filename, width = 9, height = 8)
    
  }
}



