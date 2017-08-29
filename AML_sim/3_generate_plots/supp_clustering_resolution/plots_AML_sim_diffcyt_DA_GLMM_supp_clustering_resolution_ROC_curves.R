##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: ROC curves
# 
# - supplementary: clustering resolution; diffcyt-DA-GLMM
# 
# Lukas Weber, August 2017
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/supp_clustering_resolution"

load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_GLMM_supp_clustering_resolution.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/supp_clustering_resolution"




################
# Generate plots
################

# loop over thresholds (th) and conditions (j)

# clustering resolution
resolution <- c(3, 5, 10, 20, 30, 40, 50)

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
    data <- list(k_9    = out_diffcyt_DA_GLMM_supp_resolution[[1]][[th]][[j]], 
                 k_25   = out_diffcyt_DA_GLMM_supp_resolution[[2]][[th]][[j]], 
                 k_100  = out_diffcyt_DA_GLMM_supp_resolution[[3]][[th]][[j]], 
                 k_400  = out_diffcyt_DA_GLMM_supp_resolution[[4]][[th]][[j]], 
                 k_900  = out_diffcyt_DA_GLMM_supp_resolution[[5]][[th]][[j]], 
                 k_1600 = out_diffcyt_DA_GLMM_supp_resolution[[6]][[th]][[j]], 
                 k_2500 = out_diffcyt_DA_GLMM_supp_resolution[[7]][[th]][[j]])
    
    # check
    stopifnot(all(sapply(data, function(d) all(d$spikein == data[[1]]$spikein))))
    
    # note: provide all available values
    # - 'padj' is required for threshold points on TPR-FDR curves
    # - depending on availability, plotting functions use 'score', then 'pval', then 'padj'
    cobradata <- COBRAData(pval = data.frame(k_9    = data[["k_9"]][, "p_vals"], 
                                             k_25   = data[["k_25"]][, "p_vals"], 
                                             k_100  = data[["k_100"]][, "p_vals"], 
                                             k_400  = data[["k_400"]][, "p_vals"], 
                                             k_900  = data[["k_900"]][, "p_vals"], 
                                             k_1600 = data[["k_1600"]][, "p_vals"], 
                                             k_2500 = data[["k_2500"]][, "p_vals"]), 
                           padj = data.frame(k_9    = data[["k_9"]][, "p_adj"], 
                                             k_25   = data[["k_25"]][, "p_adj"], 
                                             k_100  = data[["k_100"]][, "p_adj"], 
                                             k_400  = data[["k_400"]][, "p_adj"], 
                                             k_900  = data[["k_900"]][, "p_adj"], 
                                             k_1600 = data[["k_1600"]][, "p_adj"], 
                                             k_2500 = data[["k_2500"]][, "p_adj"]), 
                           truth = data.frame(spikein = data[["k_9"]][, "spikein"]))
    
    # calculate performance scores
    # (note: can ignore warning messages when 'padj' not available)
    cobraperf <- calculate_performance(cobradata, 
                                       binary_truth = "spikein", 
                                       aspects = c("roc"))
    
    # color scheme
    colors <- colorRampPalette(c("#c7e9c0", "darkgreen"))(7)
    
    colors <- colors[1:length(data)]
    names(colors) <- names(data)
    
    # linetypes
    #linetypes <- 1:length(data)
    #names(linetypes) <- names(data)
    
    # axis ranges
    x_range <- c(0, 0.25)
    y_range <- c(0.75, 1)
    
    # prepare plotting object
    cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colors)
    
    # re-order legend
    cobraplot <- iCOBRA:::reorder_levels(cobraplot, levels = names(data))
    
    
    
    # ----------
    # ROC curves
    # ----------
    
    # create plot
    p <- plot_roc(cobraplot, linewidth = 0.75)
    
    # with short title for multi-panel plot
    p <- p + 
      #scale_linetype_manual(values = linetypes) + 
      coord_fixed(xlim = x_range, ylim = y_range) + 
      xlab("False positive rate") + 
      ylab("True positive rate") + 
      ggtitle(paste0(cond_names[j], ", threshold ", gsub("pc$", "\\%", thresholds[th]))) + 
      theme_bw() + 
      theme(strip.text.x = element_blank())
    
    plots_ROC[[ix]] <- p
    
    # save individual panel plot
    p <- p + 
      ggtitle(paste0("diffcyt-DA-GLMM: AML-sim, clustering resolution, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th]), ": ROC curves"))
    
    fn <- file.path(DIR_PLOTS, "panels", 
                    paste0("results_diffcyt_DA_GLMM_supp_clustering_resolution_ROC_curves_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 7.5, height = 6)
    
  }
}




########################
# Save multi-panel plots
########################

# ----------
# ROC curves
# ----------

# remove duplicated annotation
plots_ROC <- lapply(plots_ROC, function(p) {
  p + theme(legend.position = "none", 
            axis.title.x = element_blank(), 
            axis.title.y = element_blank())
})

# format into grid
grid_ROC <- do.call(plot_grid, append(plots_ROC, list(labels = "AUTO", nrow = 4, ncol = 2, 
                                                      scale = 0.95, label_y = 0.975)))

# add combined axis titles
xaxis_ROC <- ggdraw() + draw_label("False positive rate", size = 12)
yaxis_ROC <- ggdraw() + draw_label("True positive rate", size = 12, angle = 90)

grid_ROC <- plot_grid(grid_ROC, xaxis_ROC, ncol = 1, rel_heights = c(50, 1))
grid_ROC <- plot_grid(yaxis_ROC, grid_ROC, nrow = 1, rel_widths = c(1, 30))

# add combined legend
legend_ROC <- get_legend(plots_ROC[[1]] + theme(legend.position = "right"))
grid_ROC <- plot_grid(grid_ROC, legend_ROC, nrow = 1, rel_widths = c(5, 1))

# add combined title
title_ROC <- ggdraw() + draw_label("diffcyt-DA-GLMM: AML-sim, clustering resolution: ROC curves", fontface = "bold")
grid_ROC <- plot_grid(title_ROC, grid_ROC, ncol = 1, rel_heights = c(1, 32))

# save plots
fn_ROC <- file.path(DIR_PLOTS, "results_diffcyt_DA_GLMM_supp_clustering_resolution_ROC_curves.pdf")
ggsave(fn_ROC, grid_ROC, width = 10, height = 14.14)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



