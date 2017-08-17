##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: ROC and TPR-FDR curves
# 
# - supplementary: clustering resolution
# 
# Lukas Weber, August 2017
##########################################################################################


library(iCOBRA)
library(ggplot2)


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/supp_resolution"

load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_edgeR_supp_resolution.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/supp_resolution"




######################################################################
# Loop over clustering resolution (k), thresholds (th), conditions (j)
######################################################################

# clustering resolution
resolution <- c(3, 5, 10, 20, 30, 40)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")


for (k in 1:length(resolution)) {
  
  for (th in 1:length(thresholds)) {
    
    for (j in 1:length(cond_names)) {
      
      # -------------------------------------
      # Pre-processing steps for iCOBRA plots
      # -------------------------------------
      
      # create 'COBRAData' object
      data <- list(diffcyt_DA_edgeR = out_diffcyt_DA_edgeR_supp_resolution[[k]][[th]][[j]])
      
      # check
      stopifnot(all(sapply(data, function(d) all(d$spikein == data[[1]]$spikein))))
      
      # note: provide all available values
      # - 'padj' is required for threshold points on TPR-FDR curves
      # - depending on availability, plotting functions use 'score', then 'pval', then 'padj'
      cobradata <- COBRAData(pval = data.frame(diffcyt_DA_edgeR = data[["diffcyt_DA_edgeR"]][, "p_vals"]), 
                             padj = data.frame(diffcyt_DA_edgeR = data[["diffcyt_DA_edgeR"]][, "p_adj"]), 
                             #score = data.frame(), 
                             truth = data.frame(spikein = data[["diffcyt_DA_edgeR"]][, "spikein"]))
      
      # calculate performance scores
      # (note: can ignore warning messages when 'padj' not available)
      cobraperf <- calculate_performance(cobradata, 
                                         binary_truth = "spikein", 
                                         aspects = c("fdrtpr", "fdrtprcurve", "roc"))
      
      # color scheme
      # modifed default "Set1" to use different yellow (#FFD92F) from colorbrewer2.org
      colors <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFD92F', '#A65628', '#F781BF')
      
      colors <- colors[1:length(data)]
      names(colors) <- names(data)
      
      # x axis labels
      x_min <- 0
      x_max <- 1
      x_labels <- seq(x_min, x_max, by = 0.1)
      
      # prepare plotting object
      cobraplot <- prepare_data_for_plot(cobraperf, 
                                         colorscheme = colors, 
                                         conditionalfill = FALSE)
      
      
      # ----------
      # ROC curves
      # ----------
      
      # create plot
      p <- plot_roc(cobraplot, linewidth = 0.75)
      
      p + 
        coord_fixed() + 
        xlab("False positive rate") + 
        ylab("True positive rate") + 
        ggtitle(paste0("ROC curves: AML-sim, ", resolution[k]^2, " clusters, ", cond_names[j], ", ", thresholds[th])) + 
        theme_bw() + 
        theme(strip.text.x = element_blank())
      
      path <- DIR_PLOTS
      filename <- file.path(path, paste0("results_supp_resolution_ROC_", resolution[k]^2, "_clusters_", thresholds[th], "_", cond_names[j], ".pdf"))
      
      ggsave(filename, width = 9, height = 8)
      
      
      # --------------
      # TPR-FDR curves
      # --------------
      
      # create plot
      p <- plot_fdrtprcurve(cobraplot, linewidth = 0.75, pointsize = 4)
      
      p + 
        scale_shape_manual(values = c(22, 21, 23)) + 
        scale_x_continuous(breaks = x_labels, labels = x_labels) + 
        coord_fixed() + 
        xlab("False discovery rate") + 
        ylab("True positive rate") + 
        ggtitle(paste0("TPR-FDR curves: AML-sim, ", resolution[k]^2, " clusters, ", cond_names[j], ", ", thresholds[th])) + 
        theme_bw() + 
        theme(strip.text.x = element_blank()) + 
        guides(color = guide_legend(override.aes = list(shape = NA)), shape = FALSE)
      
      path <- DIR_PLOTS
      filename <- file.path(path, paste0("results_supp_resolution_TPR_FDR_", resolution[k]^2, "_clusters_", thresholds[th], "_", cond_names[j], ".pdf"))
      
      ggsave(filename, width = 9, height = 8)
      
    }
  }
}




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



