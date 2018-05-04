##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: performance metrics
# - method: diffcyt methods
# 
# - supplementary results: using FlowSOM meta-clustering
# 
# Lukas Weber, May 2018
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/supp_metaclustering"

load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_edgeR_supp_metaclustering.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_voom_supp_metaclustering.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_GLMM_supp_metaclustering.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/supp_metaclustering"




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
    data <- list(diffcyt_DA_edgeR = out_diffcyt_DA_edgeR_supp_metaclustering[[th]][[j]], 
                 diffcyt_DA_voom = out_diffcyt_DA_voom_supp_metaclustering[[th]][[j]], 
                 diffcyt_DA_GLMM = out_diffcyt_DA_GLMM_supp_metaclustering[[th]][[j]])
    
    # check
    stopifnot(all(sapply(data, function(d) all(d$spikein == data[[1]]$spikein))))
    
    # note: provide all available values
    # 'padj' is required for threshold points on TPR-FDR curves
    # depending on availability, plotting functions use 'score', then 'pval', then 'padj'
    cobradata <- COBRAData(pval = data.frame(diffcyt_DA_edgeR = data[["diffcyt_DA_edgeR"]][, "p_val"], 
                                             diffcyt_DA_voom = data[["diffcyt_DA_voom"]][, "p_val"], 
                                             diffcyt_DA_GLMM = data[["diffcyt_DA_GLMM"]][, "p_val"]), 
                           padj = data.frame(diffcyt_DA_edgeR = data[["diffcyt_DA_edgeR"]][, "p_adj"], 
                                             diffcyt_DA_voom = data[["diffcyt_DA_voom"]][, "p_adj"], 
                                             diffcyt_DA_GLMM = data[["diffcyt_DA_GLMM"]][, "p_adj"]), 
                           truth = data.frame(spikein = data[["diffcyt_DA_edgeR"]][, "spikein"]))
    
    # calculate performance scores
    # (note: can ignore warning messages when 'padj' not available)
    cobraperf <- calculate_performance(cobradata, 
                                       binary_truth = "spikein", 
                                       aspects = c("roc", "fdrtpr", "fdrtprcurve", "tpr", "fpr"))
    
    # color scheme
    colors <- c("darkblue", "dodgerblue", "darkslategray3")
    
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
      ggtitle(paste0("AML-sim: meta-clustering, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th])), subtitle = "ROC curve") + 
      theme_bw() + 
      theme(strip.text.x = element_blank()) + 
      guides(color = guide_legend("method"))
    
    # save individual panel plot
    fn <- file.path(DIR_PLOTS, "panels", paste0("results_AML_sim_diffcyt_supp_metaclustering_ROC_", thresholds[th], "_", cond_names[j], ".pdf"))
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
      ggtitle(paste0("AML-sim: meta-clustering, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th])), subtitle = "TPR vs. FDR") + 
      theme_bw() + 
      theme(strip.text.x = element_blank()) + 
      guides(shape = guide_legend("FDR threshold", override.aes = list(size = 4), order = 1), 
             color = guide_legend("method", order = 2))
    
    # save individual panel plot
    fn <- file.path(DIR_PLOTS, "panels", paste0("results_AML_sim_diffcyt_supp_metaclustering_TPRFDR_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 4.75, height = 3.5)
    
    
    
    # ---------
    # TPR plots
    # ---------
    
    # create plot
    p_TPR <- 
      plot_tpr(cobraplot, pointsize = 4) + 
      scale_shape_manual(values = c(15, 19, 17), labels = c(0.01, 0.05, 0.1)) + 
      #coord_fixed() + 
      xlab("True positive rate") + 
      ggtitle(paste0("AML-sim: meta-clustering, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th])), subtitle = "TPR") + 
      theme_bw() + 
      theme(strip.text.x = element_blank(), 
            axis.text.y = element_blank()) + 
      guides(shape = guide_legend("FDR threshold", override.aes = list(size = 4), order = 1), 
             color = guide_legend("method", override.aes = list(shape = 19, size = 4), order = 2))
    
    # save individual panel plot
    fn <- file.path(DIR_PLOTS, "panels", paste0("results_AML_sim_diffcyt_supp_metaclustering_TPR_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 4.75, height = 3.5)
    
    
    
    # ---------
    # FPR plots
    # ---------
    
    # create plot
    p_FPR <- 
      plot_fpr(cobraplot, pointsize = 4) + 
      scale_shape_manual(values = c(15, 19, 17), labels = c(0.01, 0.05, 0.1)) + 
      #coord_fixed() + 
      xlab("False positive rate") + 
      ggtitle(paste0("AML-sim: meta-clustering, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th])), subtitle = "FPR") + 
      theme_bw() + 
      theme(strip.text.x = element_blank(), 
            axis.text.y = element_blank()) + 
      guides(shape = guide_legend("FDR threshold", override.aes = list(size = 4), order = 1), 
             color = guide_legend("method", override.aes = list(shape = 19, size = 4), order = 2))
    
    # save individual panel plot
    fn <- file.path(DIR_PLOTS, "panels", paste0("results_AML_sim_diffcyt_supp_metaclustering_FPR_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 4.75, height = 3.5)
    
    
    
    
    ###################
    # Multi-panel plots
    ###################
    
    plots_list <- list(p_ROC, p_TPRFDR, p_TPR, p_FPR)
    
    # modify plot elements
    plots_list <- lapply(plots_list, function(p) {
      p + 
        labs(title = p$labels$subtitle, subtitle = gsub("^.*clustering, ", "", p$labels$title)) + 
        theme(legend.position = "none")
    })
    
    plots_multi <- plot_grid(plotlist = plots_list, 
                             nrow = 1, ncol = 4, align = "hv", axis = "bl")
    
    # store in list
    plots_all[[ix]] <- plots_multi
    
    # add combined title
    title_single <- gsub(",.*$", "", p_ROC$labels$title)
    plots_title <- ggdraw() + draw_label(title_single)
    plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 7))
    
    # add combined legend
    legend_single <- get_legend(plots_list[[2]] + theme(legend.position = "right"))
    plots_multi <- plot_grid(plots_multi, legend_single, nrow = 1, rel_widths = c(6, 1))
    
    # store title and legend objects
    title_objs[[ix]] <- title_single
    legend_objs[[ix]] <- legend_single
    
    # save multi-panel plot
    fn <- file.path(DIR_PLOTS, paste0("results_AML_sim_diffcyt_supp_metaclustering_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 10, height = 2.625)
    
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
  plots_j <- plot_grid(plots_j, legend_single, nrow = 1, rel_widths = c(6, 1))
  
  # save combined plot
  fn <- file.path(DIR_PLOTS, paste0("results_AML_sim_diffcyt_supp_metaclustering_", cond_names[j], ".pdf"))
  ggsave(fn, width = 10, height = 7.25)

}




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



