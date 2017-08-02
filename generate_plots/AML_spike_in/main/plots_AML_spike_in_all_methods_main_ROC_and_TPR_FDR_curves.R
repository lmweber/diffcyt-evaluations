##########################################################################################
# Generate plots
# 
# - plot type: ROC and TPR-FDR curves
# - methods: all methods, main results
# - data set: AML-spike-in
# 
# Lukas Weber, July 2017
##########################################################################################


library(iCOBRA)
library(ggplot2)


# load saved results

DIR_RDATA <- "../../../../RData/AML_spike_in/main"

load(file.path(DIR_RDATA, "outputs_AML_spike_in_CellCnn_main.RData"))
load(file.path(DIR_RDATA, "outputs_AML_spike_in_Citrus_main.RData"))
load(file.path(DIR_RDATA, "outputs_AML_spike_in_cydar_main.RData"))
load(file.path(DIR_RDATA, "outputs_AML_spike_in_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA, "outputs_AML_spike_in_diffcyt_DA_GLMM_main.RData"))
load(file.path(DIR_RDATA, "outputs_AML_spike_in_diffcyt_DA_limma_main.RData"))

DIR_PLOTS <- "../../../../plots/AML_spike_in/all_methods/main"




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
    data_CellCnn <- out_CellCnn_main[[th]][[j]]
    data_Citrus <- out_Citrus_main[[th]][[j]]
    data_cydar <- out_cydar_main[[th]][[j]]
    data_diffcyt_DA_edgeR <- out_diffcyt_DA_edgeR_main[[th]][[j]]
    data_diffcyt_DA_GLMM <- out_diffcyt_DA_GLMM_main[[th]][[j]]
    data_diffcyt_DA_limma <- out_diffcyt_DA_limma_main[[th]][[j]]
    
    # check
    data_all <- list(data_CellCnn, data_Citrus, data_cydar, data_diffcyt_DA_edgeR, 
                     data_diffcyt_DA_GLMM, data_diffcyt_DA_limma)
    
    stopifnot(all(sapply(data_all, function(d) all(d$spikein == data_diffcyt_DA_limma$spikein))))
    
    # note: provide all available values
    # - 'padj' is required for threshold points on TPR-FDR curves
    # - depending on availability, plotting functions use 'score', then 'pval', then 'padj'
    cobradata <- COBRAData(pval = data.frame(cydar = data_cydar[, "p_vals"], 
                                             diffcyt_DA_edgeR = data_diffcyt_DA_edgeR[, "p_vals"], 
                                             diffcyt_DA_GLMM = data_diffcyt_DA_GLMM[, "p_vals"], 
                                             diffcyt_DA_limma = data_diffcyt_DA_limma[, "p_vals"]), 
                           padj = data.frame(cydar = data_cydar[, "q_vals"], 
                                             diffcyt_DA_edgeR = data_diffcyt_DA_edgeR[, "p_adj"], 
                                             diffcyt_DA_GLMM = data_diffcyt_DA_GLMM[, "p_adj"], 
                                             diffcyt_DA_limma = data_diffcyt_DA_limma[, "p_adj"]), 
                           score = data.frame(CellCnn = data_CellCnn[, "scores"], 
                                              Citrus = data_Citrus[, "scores"]), 
                           truth = data.frame(spikein = data_diffcyt_DA_limma[, "spikein"]))
    
    # calculate performance scores
    # (note: can ignore warning messages when 'padj' not available)
    cobraperf <- calculate_performance(cobradata, 
                                       binary_truth = "spikein", 
                                       aspects = c("fdrtpr", "fdrtprcurve", "roc"))
    
    # color scheme
    # modifed default "Set1" to use different yellow (#FFD92F) from colorbrewer2.org
    colors <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFD92F', '#A65628', '#F781BF')
    # re-order for legend
    #colors <- colors[c(3:7, 1:2)]
    
    # prepare plotting object
    cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colors)
    
    
    # ----------
    # ROC curves
    # ----------
    
    # create plot
    p <- plot_roc(cobraplot, linewidth = 0.75)
    
    p + 
      coord_fixed() + 
      xlab("False positive rate") + 
      ylab("True positive rate") + 
      ggtitle(paste0("ROC curves: AML-spike-in, main results, ", cond_names[j], ", ", thresholds[th])) + 
      theme_bw() + 
      theme(strip.text.x = element_blank())
    
    path <- paste0(file.path(DIR_PLOTS, thresholds[th], cond_names[j]))
    filename <- file.path(path, paste0("results_all_methods_main_ROC_curves_", thresholds[th], "_", cond_names[j], ".pdf"))
    
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
      ggtitle(paste0("TPR-FDR curves: AML-spike-in, main results, ", cond_names[j], ", ", thresholds[th])) + 
      theme_bw() + 
      theme(strip.text.x = element_blank())
    
    path <- paste0(file.path(DIR_PLOTS, thresholds[th], cond_names[j]))
    filename <- file.path(path, paste0("results_all_methods_main_TPR_FDR_curves_", thresholds[th], "_", cond_names[j], ".pdf"))
    
    ggsave(filename, width = 9, height = 8)
    
  }
}



