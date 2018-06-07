##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: performance metrics
# - method: diffcyt-DA-voom
# 
# - supplementary results: using random effects instead of fixed effects for patient IDs
# 
# Lukas Weber, June 2018
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2


# load saved results
DIR_RDATA_MAIN <- "../../../../RData/AML_sim/main"
DIR_RDATA_SUPP_RANDOM_EFFECTS_VOOM <- "../../../../RData/AML_sim/supp_random_effects_voom"

load(file.path(DIR_RDATA_MAIN, "outputs_AML_sim_diffcyt_DA_voom_main.RData"))
load(file.path(DIR_RDATA_SUPP_RANDOM_EFFECTS_VOOM, "outputs_AML_sim_diffcyt_DA_voom_supp_random_effects.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/supp_random_effects_voom"




################
# Generate plots
################

# loop over thresholds (th) and conditions (j)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# condition names
cond_names <- c("CN", "CBF")


# store plots in list
plots_all <- plots_ROC <- plots_TPRFDR <- title_objs <- legend_objs <- vector("list", length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  for (j in 1:length(cond_names)) {
    
    
    # index to store plots in list
    ix <- (j * length(thresholds)) - (length(thresholds) - th)
    
    
    # -------------------------------------
    # Pre-processing steps for iCOBRA plots
    # -------------------------------------
    
    # create 'COBRAData' object
    data <- list(main = out_diffcyt_DA_voom_main[[th]][[j]], 
                 random_effects = out_diffcyt_DA_voom_supp_random_effects[[th]][[j]])
    
    # check
    stopifnot(all(sapply(data, function(d) all(d$spikein == data[[1]]$spikein))))
    
    # note: provide all available values
    # 'padj' is required for threshold points on TPR-FDR curves
    # depending on availability, plotting functions use 'score', then 'pval', then 'padj'
    cobradata <- COBRAData(pval = data.frame(main = data[["main"]][, "p_val"], 
                                             random_effects = data[["random_effects"]][, "p_val"]), 
                           padj = data.frame(main = data[["main"]][, "p_adj"], 
                                             random_effects = data[["random_effects"]][, "p_adj"]), 
                           truth = data.frame(spikein = data[["main"]][, "spikein"]))
    
    # calculate performance scores
    # (note: can ignore warning messages when 'padj' not available)
    cobraperf <- calculate_performance(cobradata, 
                                       binary_truth = "spikein", 
                                       aspects = c("roc", "fdrtpr", "fdrtprcurve", "tpr", "fpr"))
    
    # color scheme
    colors <- rep("dodgerblue", 2)
    
    colors <- colors[1:length(data)]
    names(colors) <- names(data)
    
    # linetypes
    linetypes <- c("solid", "dashed")
    names(linetypes) <- names(data)
    
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
    p_ROC <- 
      plot_roc(cobraplot, linewidth = 0.75) + 
      scale_linetype_manual(values = linetypes, guide = FALSE) + 
      coord_fixed() + 
      xlab("False positive rate") + 
      ylab("True positive rate") + 
      ggtitle(paste0("AML-sim: diffcyt-DA-voom, random effects, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th])), subtitle = "ROC curve") + 
      theme_bw() + 
      theme(strip.text.x = element_blank()) + 
      guides(color = guide_legend("method"), linetype = guide_legend("method"))
    
    plots_ROC[[ix]] <- p_ROC
    
    # save individual panel plot
    fn <- file.path(DIR_PLOTS, "panels", paste0("results_AML_sim_diffcyt_DA_voom_supp_random_effects_ROC_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 4.75, height = 3.5)
    
    
    
    # --------------
    # TPR-FDR curves
    # --------------
    
    # create plot
    p_TPRFDR <- 
      plot_fdrtprcurve(cobraplot, linewidth = 0.75, pointsize = 4) + 
      scale_shape_manual(values = c(15, 19, 17), labels = c(0.01, 0.05, 0.1)) + 
      scale_linetype_manual(values = linetypes, guide = FALSE) + 
      coord_fixed() + 
      xlab("False discovery rate") + 
      ylab("True positive rate") + 
      scale_x_continuous(breaks = seq(0, 1, by = 0.2)) + 
      ggtitle(paste0("AML-sim: diffcyt-DA-voom, random effects, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th])), subtitle = "TPR vs. FDR") + 
      theme_bw() + 
      theme(strip.text.x = element_blank()) + 
      guides(shape = guide_legend("FDR threshold", override.aes = list(size = 4), order = 1), 
             color = guide_legend("method", order = 2, override.aes = list(shape = NA)), 
             linetype = guide_legend("method", order = 2))
    
    plots_TPRFDR[[ix]] <- p_TPRFDR
    
    # save individual panel plot
    fn <- file.path(DIR_PLOTS, "panels", paste0("results_AML_sim_diffcyt_DA_voom_supp_random_effects_TPRFDR_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 4.75, height = 3.5)
    
  }
}




#######################################
# Combine plots: ROC and TPR-FDR curves
#######################################

# ----------
# ROC curves
# ----------

plots_ROC_keep <- plots_ROC

# modify plot elements
plots_ROC <- lapply(plots_ROC, function(p) {
  p + 
    labs(title = gsub("^.*effects, ", "", p$labels$title)) + 
    theme(legend.position = "none")
})

# format into grid
plots_ROC <- plot_grid(plotlist = plots_ROC, nrow = 2, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_ROC <- ggdraw() + draw_label(gsub("effects,.*$", "effects", plots_ROC_keep[[1]]$labels$title), fontface = "bold")
grid_ROC <- plot_grid(title_ROC, plots_ROC, ncol = 1, rel_heights = c(1, 20))

# add combined legend
legend_ROC <- get_legend(plots_ROC_keep[[1]] + theme(legend.position = "right", 
                                                     legend.title = element_text(size = 12, face = "bold"), 
                                                     legend.text = element_text(size = 12)))
grid_ROC <- plot_grid(grid_ROC, legend_ROC, nrow = 1, rel_widths = c(3.5, 1))

# save plots
fn_ROC <- file.path(DIR_PLOTS, paste0("results_AML_sim_diffcyt_DA_voom_supp_random_effects_performance_ROC.pdf"))
ggsave(fn_ROC, grid_ROC, width = 8, height = 4.9)


# --------------
# TPR-FDR curves
# --------------

plots_TPRFDR_keep <- plots_TPRFDR

# modify plot elements
plots_TPRFDR <- lapply(plots_TPRFDR, function(p) {
  p + 
    labs(title = gsub("^.*effects, ", "", p$labels$title)) + 
    theme(legend.position = "none")
})

# format into grid
plots_TPRFDR <- plot_grid(plotlist = plots_TPRFDR, nrow = 2, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_TPRFDR <- ggdraw() + draw_label(gsub("effects,.*$", "effects", plots_TPRFDR_keep[[1]]$labels$title), fontface = "bold")
grid_TPRFDR <- plot_grid(title_TPRFDR, plots_TPRFDR, ncol = 1, rel_heights = c(1, 20))

# add combined legend
legend_TPRFDR <- get_legend(plots_TPRFDR_keep[[1]] + theme(legend.position = "right", 
                                                           legend.title = element_text(size = 12, face = "bold"), 
                                                           legend.text = element_text(size = 12)))
grid_TPRFDR <- plot_grid(grid_TPRFDR, legend_TPRFDR, nrow = 1, rel_widths = c(3.5, 1))

# save plots
fn_TPRFDR <- file.path(DIR_PLOTS, paste0("results_AML_sim_diffcyt_DA_voom_supp_random_effects_performance_TPRFDR.pdf"))
ggsave(fn_TPRFDR, grid_TPRFDR, width = 8, height = 4.9)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



