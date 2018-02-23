##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: performance metrics
# - method: diffcyt-DA-GLMM
# 
# - supplementary results: 'less distinct' data sets
# 
# Lukas Weber, February 2018
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2


# load saved results
DIR_RDATA_MAIN <- "../../../../RData/AML_sim/main"
DIR_RDATA_SUPP_LESS_DISTINCT <- "../../../../RData/AML_sim/supp_less_distinct"

load(file.path(DIR_RDATA_MAIN, "outputs_AML_sim_diffcyt_DA_GLMM_main.RData"))
load(file.path(DIR_RDATA_SUPP_LESS_DISTINCT, "outputs_AML_sim_diffcyt_DA_GLMM_supp_less_distinct.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/supp_less_distinct"




################
# Generate plots
################

# loop over thresholds (th) and conditions (j)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# condition names
cond_names <- c("CN", "CBF")


# store plots in list
plots_all <- vector("list", length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  for (j in 1:length(cond_names)) {
    
    
    # index to store plots in list
    ix <- (j * length(thresholds)) - (length(thresholds) - th)
    
    
    # -------------------------------------
    # Pre-processing steps for iCOBRA plots
    # -------------------------------------
    
    # create 'COBRAData' object
    data <- list(main = out_diffcyt_DA_GLMM_main[[th]][[j]], 
                 less_50pc = out_diffcyt_DA_GLMM_supp_less_distinct[["less_50pc"]][[th]][[j]], 
                 less_75pc = out_diffcyt_DA_GLMM_supp_less_distinct[["less_75pc"]][[th]][[j]])
    
    # check
    stopifnot(all(sapply(data, function(d) all(d$spikein == data[[1]]$spikein))))
    
    # note: provide all available values
    # 'padj' is required for threshold points on TPR-FDR curves
    # depending on availability, plotting functions use 'score', then 'pval', then 'padj'
    cobradata <- COBRAData(pval = data.frame(main = data[["main"]][, "p_vals"], 
                                             less_50pc = data[["less_50pc"]][, "p_vals"], 
                                             less_75pc = data[["less_75pc"]][, "p_vals"]), 
                           padj = data.frame(main = data[["main"]][, "p_adj"], 
                                             less_50pc = data[["less_50pc"]][, "p_adj"], 
                                             less_75pc = data[["less_75pc"]][, "p_adj"]), 
                           truth = data.frame(spikein = data[["main"]][, "spikein"]))
    
    # calculate performance scores
    # (note: can ignore warning messages when 'padj' not available)
    cobraperf <- calculate_performance(cobradata, 
                                       binary_truth = "spikein", 
                                       aspects = "roc")
    
    # color scheme
    colors <- c("darkslategray3")
    
    colors <- colors[1:length(data)]
    names(colors) <- names(data)
    
    # linetypes
    linetypes <- c("solid", "dashed", "dotted")
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
      ggtitle(paste0("AML-sim: 'less distinct' benchmark data, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th])), subtitle = "ROC curve") + 
      theme_bw() + 
      theme(strip.text.x = element_blank()) + 
      guides(color = guide_legend("distinctness"), linetype = guide_legend("distinctness"))
    
    
    # modify plot to allow multiple legends
    d_plot <- p_ROC$data
    distinctness_names <- c("main", "less_50pc", "less_75pc")
    d_plot$distinctness_names <- factor(d_plot$method, levels = distinctness_names)
    d_plot$method_names <- as.factor("diffcyt-DA-GLMM")
    
    p <- 
      ggplot(d_plot, aes(x = FPR, y = TPR, linetype = distinctness_names, color = method_names)) + 
      geom_line(lwd = 0.75) + 
      scale_color_manual(values = unname(colors[1])) + 
      scale_linetype_manual(values = linetypes) + 
      coord_fixed() + 
      xlim(c(0, 1)) + 
      ylim(c(0, 1)) + 
      xlab("False positive rate") + 
      ylab("True positive rate") + 
      ggtitle(paste0("AML-sim: 'less distinct' benchmark data, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th])), subtitle = "ROC curve") + 
      theme_bw() + 
      guides(linetype = guide_legend("distinctness", order = 1), 
             color = guide_legend("method", order = 2))
    
    plots_all[[ix]] <- p
    
    # save individual panel plot
    fn <- file.path(DIR_PLOTS, "panels", 
                    paste0("results_AML_sim_diffcyt_DA_GLMM_supp_less_distinct_ROC_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 4.75, height = 3.5)
    
  }
}




########################
# Save multi-panel plots
########################

# modify plot elements
plots_multi <- lapply(plots_all, function(p) {
  p + 
    labs(title = gsub("^.*data, ", "", p$labels$title)) + 
    theme(legend.position = "none")
})

# format into grid
plots_multi <- plot_grid(plotlist = plots_multi, nrow = 2, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_multi <- ggdraw() + draw_label(gsub(",.*$", "", plots_all[[1]]$labels$title), fontface = "bold")
grid_multi <- plot_grid(title_multi, plots_multi, ncol = 1, rel_heights = c(1, 20))

# add combined legend
legend_multi <- get_legend(plots_all[[1]] + theme(legend.position = "right", 
                                                  legend.title = element_text(size = 12, face = "bold"), 
                                                  legend.text = element_text(size = 12)))
grid_multi <- plot_grid(grid_multi, legend_multi, nrow = 1, rel_widths = c(3.5, 1))

# save plots
fn_multi <- file.path(DIR_PLOTS, paste0("results_AML_sim_diffcyt_DA_GLMM_supp_less_distinct_ROC.pdf"))
ggsave(fn_multi, grid_multi, width = 8, height = 4.9)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



