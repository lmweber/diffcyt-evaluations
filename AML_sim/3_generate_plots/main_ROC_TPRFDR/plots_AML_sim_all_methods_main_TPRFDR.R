##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: TPR-FDR curves
# - method: diffcyt methods
# 
# - main results
# 
# Lukas Weber, October 2017
##########################################################################################


# note: showing 'diffcyt' methods only


library(iCOBRA)
library(ggplot2)
library(cowplot)


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/main"

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
plots_TPRFDR <- plots_TPRFDR_zoom <- vector("list", length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  for (j in 1:length(cond_names)) {
    
    # -------------------------------------
    # Pre-processing steps for iCOBRA plots
    # -------------------------------------
    
    # index to store plots sequentially in list
    ix <- (th * length(cond_names)) - (length(cond_names) - j)
    
    # create 'COBRAData' object
    data <- list(diffcyt_DA_edgeR = out_diffcyt_DA_edgeR_main[[th]][[j]], 
                 diffcyt_DA_GLMM = out_diffcyt_DA_GLMM_main[[th]][[j]], 
                 diffcyt_DA_limma = out_diffcyt_DA_limma_main[[th]][[j]])
    
    # check
    stopifnot(all(sapply(data, function(d) all(d$spikein == data[[1]]$spikein))))
    
    # note: provide all available values
    # - 'padj' is required for threshold points on TPR-FDR curves
    # - depending on availability, plotting functions use 'score', then 'pval', then 'padj'
    cobradata <- COBRAData(pval = data.frame(diffcyt_DA_edgeR = data[["diffcyt_DA_edgeR"]][, "p_vals"], 
                                             diffcyt_DA_GLMM = data[["diffcyt_DA_GLMM"]][, "p_vals"], 
                                             diffcyt_DA_limma = data[["diffcyt_DA_limma"]][, "p_vals"]), 
                           padj = data.frame(diffcyt_DA_edgeR = data[["diffcyt_DA_edgeR"]][, "p_adj"], 
                                             diffcyt_DA_GLMM = data[["diffcyt_DA_GLMM"]][, "p_adj"], 
                                             diffcyt_DA_limma = data[["diffcyt_DA_limma"]][, "p_adj"]), 
                           truth = data.frame(spikein = data[["diffcyt_DA_limma"]][, "spikein"]))
    
    # calculate performance scores
    # (note: can ignore warning messages when 'padj' not available)
    cobraperf <- calculate_performance(cobradata, 
                                       binary_truth = "spikein", 
                                       aspects = c("fdrtpr", "fdrtprcurve"))
    
    # color scheme
    colors <- c("darkblue", "deepskyblue2", "darkslategray2")
    
    colors <- colors[1:length(data)]
    names(colors) <- names(data)
    
    # prepare plotting object
    cobraplot <- prepare_data_for_plot(cobraperf, 
                                       colorscheme = colors, 
                                       conditionalfill = FALSE)
    
    # re-order legend
    cobraplot <- reorder_levels(cobraplot, levels = names(data))
    
    
    
    # --------------------
    # TPR-FDR curves: full
    # --------------------
    
    # axis limits
    x_max <- 1
    y_min <- 0
    
    # create plot
    p <- plot_fdrtprcurve(cobraplot, linewidth = 0.75, pointsize = 4)
    
    # with short title for multi-panel plot
    p <- p + 
      scale_shape_manual(values = c(22, 21, 23)) + 
      coord_fixed() + 
      xlim(c(0, x_max)) + 
      ylim(c(y_min, 1)) + 
      xlab("False discovery rate") + 
      ylab("True positive rate") + 
      ggtitle(paste0(cond_names[j], ", threshold ", gsub("pc$", "\\%", thresholds[th]))) + 
      theme_bw() + 
      theme(strip.text.x = element_blank()) + 
      guides(color = guide_legend("method", override.aes = list(shape = NA)), shape = FALSE)
    
    plots_TPRFDR[[ix]] <- p
    
    # save individual panel plot
    p <- p + 
      ggtitle(paste0("AML-sim, main results, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th]), ": TPR vs. FDR"))
    
    fn <- file.path(DIR_PLOTS, "panels", 
                    paste0("results_AML_sim_diffcyt_methods_main_TPRFDR_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 6.5, height = 5.25)
    
    
    
    # ----------------------
    # TPR-FDR curves: zoomed
    # ----------------------
    
    # axis limits
    x_max_zoom <- 0.2
    y_min_zoom <- 0.8
    
    # create plot
    p <- plot_fdrtprcurve(cobraplot, linewidth = 0.75, pointsize = 4)
    
    # with short title for multi-panel plot
    p <- p + 
      scale_shape_manual(values = c(22, 21, 23)) + 
      coord_fixed() + 
      xlim(c(0, x_max_zoom)) + 
      ylim(c(y_min_zoom, 1)) + 
      xlab("False discovery rate") + 
      ylab("True positive rate") + 
      ggtitle(paste0(cond_names[j], ", threshold ", gsub("pc$", "\\%", thresholds[th]))) + 
      theme_bw() + 
      theme(strip.text.x = element_blank()) + 
      guides(color = guide_legend("method", override.aes = list(shape = NA)), shape = FALSE)
    
    plots_TPRFDR_zoom[[ix]] <- p
    
    # save individual panel plot
    p <- p + 
      ggtitle(paste0("AML-sim, main results, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th]), ": TPR vs. FDR"))
    
    fn <- file.path(DIR_PLOTS, "panels", 
                    paste0("results_AML_sim_diffcyt_methods_main_TPRFDR_", thresholds[th], "_", cond_names[j], "_zoom.pdf"))
    ggsave(fn, width = 6.5, height = 5.25)
    
  }
}




########################
# Save multi-panel plots
########################

# --------------------
# TPR-FDR curves: full
# --------------------

# re-order plots to fill each condition by row
ord <- c(2 * (1:4) - 1, 2 * (1:4))
plots_TPRFDR <- plots_TPRFDR[ord]

# modify plot elements
plots_TPRFDR <- lapply(plots_TPRFDR, function(p) {
  p + theme(legend.position = "none", 
            axis.title.x = element_blank(), axis.title.y = element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

# format into grid
grid_TPRFDR <- do.call(plot_grid, append(plots_TPRFDR, list(nrow = 2, ncol = 4, 
                                                            align = "hv", axis = "bl", 
                                                            scale = 0.975)))

# add combined axis titles
xaxis_TPRFDR <- ggdraw() + draw_label("False discovery rate", size = 14)
yaxis_TPRFDR <- ggdraw() + draw_label("True positive rate", size = 14, angle = 90)

grid_TPRFDR <- plot_grid(grid_TPRFDR, xaxis_TPRFDR, ncol = 1, rel_heights = c(12, 1))
grid_TPRFDR <- plot_grid(yaxis_TPRFDR, grid_TPRFDR, nrow = 1, rel_widths = c(1, 30))

# add combined title
title_TPRFDR <- ggdraw() + draw_label("AML-sim, main results: TPR vs. FDR", fontface = "bold")
grid_TPRFDR <- plot_grid(title_TPRFDR, grid_TPRFDR, ncol = 1, rel_heights = c(1, 13))

# add combined legend
legend_TPRFDR <- get_legend(plots_TPRFDR[[1]] + theme(legend.position = "right", 
                                                legend.title = element_text(size = 12, face = "bold"), 
                                                legend.text = element_text(size = 12)))
grid_TPRFDR <- plot_grid(grid_TPRFDR, legend_TPRFDR, nrow = 1, rel_widths = c(4.8, 1))

# save plots
fn_TPRFDR <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_methods_main_TPRFDR.pdf")
ggsave(fn_TPRFDR, grid_TPRFDR, width = 11, height = 5.4)



# ----------------------
# TPR-FDR curves: zoomed
# ----------------------

# re-order plots to fill each condition by row
ord <- c(2 * (1:4) - 1, 2 * (1:4))
plots_TPRFDR_zoom <- plots_TPRFDR_zoom[ord]

# modify plot elements
plots_TPRFDR_zoom <- lapply(plots_TPRFDR_zoom, function(p) {
  p + theme(legend.position = "none", 
            axis.title.x = element_blank(), axis.title.y = element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

# format into grid
grid_TPRFDR_zoom <- do.call(plot_grid, append(plots_TPRFDR_zoom, list(nrow = 2, ncol = 4, 
                                                                      align = "hv", axis = "bl", 
                                                                      scale = 0.975)))

# add combined axis titles
xaxis_TPRFDR_zoom <- ggdraw() + draw_label("False discovery rate", size = 14)
yaxis_TPRFDR_zoom <- ggdraw() + draw_label("True positive rate", size = 14, angle = 90)

grid_TPRFDR_zoom <- plot_grid(grid_TPRFDR_zoom, xaxis_TPRFDR_zoom, ncol = 1, rel_heights = c(12, 1))
grid_TPRFDR_zoom <- plot_grid(yaxis_TPRFDR_zoom, grid_TPRFDR_zoom, nrow = 1, rel_widths = c(1, 30))

# add combined title
title_TPRFDR_zoom <- ggdraw() + draw_label(paste0("AML-sim, main results: TPR vs. FDR (FDR < ", x_max_zoom, ", TPR > ", y_min_zoom, ")"), fontface = "bold")
grid_TPRFDR_zoom <- plot_grid(title_TPRFDR_zoom, grid_TPRFDR_zoom, ncol = 1, rel_heights = c(1, 13))

# add combined legend
legend_TPRFDR_zoom <- get_legend(plots_TPRFDR_zoom[[1]] + theme(legend.position = "right", 
                                                                legend.title = element_text(size = 12, face = "bold"), 
                                                                legend.text = element_text(size = 12)))
grid_TPRFDR_zoom <- plot_grid(grid_TPRFDR_zoom, legend_TPRFDR_zoom, nrow = 1, rel_widths = c(4.8, 1))

# save plots
fn_TPRFDR_zoom <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_methods_main_TPRFDR_zoom.pdf")
ggsave(fn_TPRFDR_zoom, grid_TPRFDR_zoom, width = 11, height = 5.4)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



