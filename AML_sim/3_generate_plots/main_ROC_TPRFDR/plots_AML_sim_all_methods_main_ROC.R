##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: ROC curves
# - method: all methods
# 
# - main results
# 
# Lukas Weber, October 2017
##########################################################################################


# note: all methods except CellCnn use lineage markers only; CellCnn uses all markers


library(iCOBRA)
library(ggplot2)
library(cowplot)


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/main"

load(file.path(DIR_RDATA, "outputs_AML_sim_CellCnn_main.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_Citrus_main.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_cydar_main.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_GLMM_main.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_limma_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/main_ROC_TPRFDR"




################
# Generate plots
################

# loop over thresholds (th) and conditions (j)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")


# store plots in list
plots_ROC <- vector("list", length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  for (j in 1:length(cond_names)) {
    
    # -------------------------------------
    # Pre-processing steps for iCOBRA plots
    # -------------------------------------
    
    # index to store plots sequentially in list
    ix <- (th * length(cond_names)) - (length(cond_names) - j)
    
    # create 'COBRAData' object
    data <- list(CellCnn = out_CellCnn_main[[th]][[j]], 
                 Citrus = out_Citrus_main[[th]][[j]], 
                 cydar = out_cydar_main[[th]][[j]], 
                 diffcyt_DA_edgeR = out_diffcyt_DA_edgeR_main[[th]][[j]], 
                 diffcyt_DA_GLMM = out_diffcyt_DA_GLMM_main[[th]][[j]], 
                 diffcyt_DA_limma = out_diffcyt_DA_limma_main[[th]][[j]])
    
    # check
    stopifnot(all(sapply(data, function(d) all(d$spikein == data[[1]]$spikein))))
    
    # note: provide all available values
    # - 'padj' is required for threshold points on TPR-FDR curves
    # - depending on availability, plotting functions use 'score', then 'pval', then 'padj'
    cobradata <- COBRAData(pval = data.frame(cydar = data[["cydar"]][, "p_vals"], 
                                             diffcyt_DA_edgeR = data[["diffcyt_DA_edgeR"]][, "p_vals"], 
                                             diffcyt_DA_GLMM = data[["diffcyt_DA_GLMM"]][, "p_vals"], 
                                             diffcyt_DA_limma = data[["diffcyt_DA_limma"]][, "p_vals"]), 
                           padj = data.frame(cydar = data[["cydar"]][, "q_vals"], 
                                             diffcyt_DA_edgeR = data[["diffcyt_DA_edgeR"]][, "p_adj"], 
                                             diffcyt_DA_GLMM = data[["diffcyt_DA_GLMM"]][, "p_adj"], 
                                             diffcyt_DA_limma = data[["diffcyt_DA_limma"]][, "p_adj"]), 
                           score = data.frame(CellCnn = data[["CellCnn"]][, "scores"], 
                                              Citrus = data[["Citrus"]][, "scores"]), 
                           truth = data.frame(spikein = data[["diffcyt_DA_limma"]][, "spikein"]))
    
    # calculate performance scores
    # (note: can ignore warning messages when 'padj' not available)
    cobraperf <- calculate_performance(cobradata, 
                                       binary_truth = "spikein", 
                                       aspects = "roc")
    
    # color scheme
    colors <- c("mediumorchid3", "gold", "salmon", "darkblue", "deepskyblue2", "darkslategray2")
    
    # alternative: greens
    #colors <- c("mediumorchid3", "gold", "salmon", "darkgreen", "limegreen", "olivedrab1")
    # alternative: modifed default "Set1" to use different yellow (#FFD92F) from colorbrewer2.org
    #colors <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFD92F', '#A65628', '#F781BF')
    
    colors <- colors[1:length(data)]
    names(colors) <- names(data)
    
    # prepare plotting object
    cobraplot <- prepare_data_for_plot(cobraperf, 
                                       colorscheme = colors, 
                                       conditionalfill = FALSE)
    
    # re-order legend
    cobraplot <- reorder_levels(cobraplot, levels = names(data))
    
    
    
    # ----------
    # ROC curves
    # ----------
    
    # create plot
    p <- plot_roc(cobraplot, linewidth = 0.75)
    
    # with short title for multi-panel plot
    p <- p + 
      coord_fixed() + 
      xlab("False positive rate") + 
      ylab("True positive rate") + 
      ggtitle(paste0(cond_names[j], ", threshold ", gsub("pc$", "\\%", thresholds[th]))) + 
      theme_bw() + 
      theme(strip.text.x = element_blank()) + 
      guides(color = guide_legend("method"))
    
    plots_ROC[[ix]] <- p
    
    # save individual panel plot
    p <- p + 
      ggtitle(paste0("AML-sim, main results, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th]), ": ROC curves"))
    
    fn <- file.path(DIR_PLOTS, "panels", 
                    paste0("results_AML_sim_all_methods_main_ROC_curves_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 7.5, height = 6)
    
  }
}




########################
# Save multi-panel plots
########################

# re-order plots to fill each condition by row
ord <- c(2 * (1:4) - 1, 2 * (1:4))
plots_ROC <- plots_ROC[ord]

# modify plot elements
plots_ROC <- lapply(plots_ROC, function(p) {
  p + theme(legend.position = "none", 
            axis.title.x = element_blank(), axis.title.y = element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

# format into grid
grid_ROC <- do.call(plot_grid, append(plots_ROC, list(nrow = 2, ncol = 4, 
                                                      align = "hv", axis = "bl", 
                                                      scale = 0.975)))

# add combined axis titles
xaxis_ROC <- ggdraw() + draw_label("False positive rate", size = 14)
yaxis_ROC <- ggdraw() + draw_label("True positive rate", size = 14, angle = 90)

grid_ROC <- plot_grid(grid_ROC, xaxis_ROC, ncol = 1, rel_heights = c(12, 1))
grid_ROC <- plot_grid(yaxis_ROC, grid_ROC, nrow = 1, rel_widths = c(1, 30))

# add combined title
title_ROC <- ggdraw() + draw_label("AML-sim, main results: ROC curves", fontface = "bold")
grid_ROC <- plot_grid(title_ROC, grid_ROC, ncol = 1, rel_heights = c(1, 13))

# add combined legend
legend_ROC <- get_legend(plots_ROC[[1]] + theme(legend.position = "right", 
                                                legend.title = element_text(size = 12, face = "bold"), 
                                                legend.text = element_text(size = 12)))
grid_ROC <- plot_grid(grid_ROC, legend_ROC, nrow = 1, rel_widths = c(4.8, 1))

# save plots
fn_ROC <- file.path(DIR_PLOTS, "results_AML_sim_all_methods_main_ROC_curves.pdf")
ggsave(fn_ROC, grid_ROC, width = 11, height = 5.4)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



