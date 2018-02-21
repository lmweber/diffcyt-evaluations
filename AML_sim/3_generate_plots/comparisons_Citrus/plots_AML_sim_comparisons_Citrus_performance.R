##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: performance metrics
# - method: Citrus
# 
# - main results
# 
# Lukas Weber, February 2018
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2


# load saved results
DIR_RDATA_MAIN <- "../../../../RData/AML_sim/main"
DIR_RDATA_CITRUS <- "../../../../RData/AML_sim/comparisons_Citrus"

load(file.path(DIR_RDATA_MAIN, "outputs_AML_sim_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA_MAIN, "outputs_AML_sim_diffcyt_DA_limma_main.RData"))
load(file.path(DIR_RDATA_MAIN, "outputs_AML_sim_diffcyt_DA_GLMM_main.RData"))
load(file.path(DIR_RDATA_CITRUS, "outputs_AML_sim_Citrus_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/comparisons_Citrus"

# path where automatically generated Citrus plots are saved
DIR_CITRUS_FILES <- "../../../../Citrus_files/AML_sim/main"




##########################
# Copy Citrus output plots
##########################

# copy automatically generated Citrus plots (for one threshold and condition only)

plots_Citrus <- c(file.path(DIR_CITRUS_FILES, "5pc", "CN", "citrusOutput", "defaultCondition", 
                            "markerPlotsAll.pdf"), 
                  file.path(DIR_CITRUS_FILES, "5pc", "CN", "citrusOutput", "defaultCondition", "glmnet_results", 
                            c("clusters-cv_min.pdf", "featurePlots_cv.min.pdf", "features-cv_min.pdf")))

cmds <- paste("cp", plots_Citrus, DIR_PLOTS)

for (i in 1:length(cmds)) {
  system(cmds[i])
}




################
# Generate plots
################

# loop over thresholds (th) and conditions (j)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# condition names
cond_names <- c("CN", "CBF")


# store plots in list
plots_all <- title_objs <- legend_objs <- vector("list", length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  for (j in 1:length(cond_names)) {
    
    
    # index to store plots in list
    ix <- (j * length(thresholds)) - (length(thresholds) - th)
    
    
    # -------------------------------------
    # Pre-processing steps for iCOBRA plots
    # -------------------------------------
    
    # create 'COBRAData' object
    data <- list(diffcyt_DA_edgeR = out_diffcyt_DA_edgeR_main[[th]][[j]], 
                 diffcyt_DA_limma = out_diffcyt_DA_limma_main[[th]][[j]], 
                 diffcyt_DA_GLMM = out_diffcyt_DA_GLMM_main[[th]][[j]], 
                 Citrus = out_Citrus_main[[th]][[j]])
    
    # check
    stopifnot(all(sapply(data, function(d) all(d$spikein == data[[1]]$spikein))))
    
    # note: provide all available values
    # 'padj' is required for threshold points on TPR-FDR curves
    # depending on availability, plotting functions use 'score', then 'pval', then 'padj'
    cobradata <- COBRAData(pval = data.frame(diffcyt_DA_edgeR = data[["diffcyt_DA_edgeR"]][, "p_vals"], 
                                             diffcyt_DA_limma = data[["diffcyt_DA_limma"]][, "p_vals"], 
                                             diffcyt_DA_GLMM = data[["diffcyt_DA_GLMM"]][, "p_vals"]), 
                           padj = data.frame(diffcyt_DA_edgeR = data[["diffcyt_DA_edgeR"]][, "p_adj"], 
                                             diffcyt_DA_limma = data[["diffcyt_DA_limma"]][, "p_adj"], 
                                             diffcyt_DA_GLMM = data[["diffcyt_DA_GLMM"]][, "p_adj"]), 
                           score = data.frame(Citrus = data[["Citrus"]][, "scores"]), 
                           truth = data.frame(spikein = data[["diffcyt_DA_edgeR"]][, "spikein"]))
    
    # calculate performance scores
    # (note: can ignore warning messages when 'padj' not available)
    cobraperf <- calculate_performance(cobradata, 
                                       binary_truth = "spikein", 
                                       aspects = c("roc", "fdrtpr", "fdrtprcurve", "tpr", "fpr"))
    
    # color scheme
    colors <- c("darkblue", "dodgerblue", "darkslategray3", "gold2")
    
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
    p_ROC <- 
      plot_roc(cobraplot, linewidth = 0.75) + 
      coord_fixed() + 
      xlab("False positive rate") + 
      ylab("True positive rate") + 
      ggtitle(paste0("AML-sim: comparisons, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th])), subtitle = "ROC curve") + 
      theme_bw() + 
      theme(strip.text.x = element_blank()) + 
      guides(color = guide_legend("method"))
    
    # save individual panel plot
    fn <- file.path(DIR_PLOTS, "panels", paste0("results_AML_sim_comparisons_Citrus_main_ROC_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 4.75, height = 3.5)
    
    
    
    # --------------
    # TPR-FDR curves
    # --------------
    
    # create plot
    p_TPRFDR <- 
      plot_fdrtprcurve(cobraplot, linewidth = 0.75, pointsize = 4) + 
      scale_shape_manual(values = c(15, 19, 17), labels = c(0.01, 0.05, 0.1)) + 
      coord_fixed() + 
      xlab("False discovery rate") + 
      ylab("True positive rate") + 
      scale_x_continuous(breaks = seq(0, 1, by = 0.2)) + 
      ggtitle(paste0("AML-sim: comparisons, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th])), subtitle = "TPR vs. FDR") + 
      theme_bw() + 
      theme(strip.text.x = element_blank()) + 
      guides(shape = guide_legend("FDR threshold", override.aes = list(size = 4), order = 1), 
             color = guide_legend("method", order = 2))
    
    # save individual panel plot
    fn <- file.path(DIR_PLOTS, "panels", paste0("results_AML_sim_comparisons_Citrus_main_TPRFDR_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 4.75, height = 3.5)
    
    
    
    
    ###################
    # Multi-panel plots
    ###################
    
    # note: ROC and TPR-FDR plots only (for Citrus)
    
    
    plots_list <- list(p_ROC, p_TPRFDR)
    
    # modify plot elements
    plots_list <- lapply(plots_list, function(p) {
      p + 
        labs(title = p$labels$subtitle, subtitle = gsub("^.*comparisons, ", "", p$labels$title)) + 
        theme(legend.position = "none")
    })
    
    plots_multi <- plot_grid(plotlist = plots_list, 
                             nrow = 1, ncol = 2, align = "hv", axis = "bl")
    
    # store in list
    plots_all[[ix]] <- plots_multi
    
    # add combined title
    title_single <- gsub("comparisons,.*$", "performance comparisons", p_ROC$labels$title)
    plots_title <- ggdraw() + draw_label(title_single)
    plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 7))
    
    # add combined legend
    legend_single <- get_legend(plots_list[[2]] + theme(legend.position = "right"))
    plots_multi <- plot_grid(plots_multi, legend_single, nrow = 1, rel_widths = c(3, 1))
    
    # store title and legend objects
    title_objs[[ix]] <- title_single
    legend_objs[[ix]] <- legend_single
    
    # save multi-panel plot
    fn <- file.path(DIR_PLOTS, paste0("results_AML_sim_comparisons_Citrus_main_performance_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 6.5, height = 2.625)
    
  }
}




###############################################
# Combine plots: one set of plots per condition
###############################################

ix_cnd <- list(CN = 1:3, CBF = 4:6)


for (j in 1:length(cond_names)) {
  
  # combine plots
  plots_j <- plot_grid(plotlist = plots_all[ix_cnd[[j]]], 
                       nrow = 3, ncol = 1, align = "hv", axis = "bl")
  
  # add combined title
  title_j <- paste0(title_objs[[1]], ", ", cond_names[j], " vs. healthy")
  plots_title <- ggdraw() + draw_label(title_j)
  plots_j <- plot_grid(plots_title, plots_j, ncol = 1, rel_heights = c(1, 16))
  
  # add combined legend
  legend_single <- legend_objs[[1]]
  plots_j <- plot_grid(plots_j, legend_single, nrow = 1, rel_widths = c(3, 1))
  
  # save combined plot
  fn <- file.path(DIR_PLOTS, paste0("results_AML_sim_comparisons_Citrus_main_performance_", cond_names[j], ".pdf"))
  ggsave(fn, width = 6.5, height = 7.25)

}




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



