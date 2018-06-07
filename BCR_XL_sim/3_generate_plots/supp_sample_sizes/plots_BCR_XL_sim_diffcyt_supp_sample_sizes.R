##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: performance metrics
# - method: diffcyt methods
# 
# - supplementary results: smaller sample sizes
# 
# Lukas Weber, June 2018
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL_sim/supp_sample_sizes"

load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_limma_supp_sample_sizes.RData"))
load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_LMM_supp_sample_sizes.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/supp_sample_sizes"




################
# Generate plots
################

# -------------------------------------
# Pre-processing steps for iCOBRA plots
# -------------------------------------

# note: create separate objects for each sample size


# create 'COBRAData' objects

data_2vs2 <- list(diffcyt_DS_limma = out_diffcyt_DS_limma_supp_sample_sizes[[1]], 
                  diffcyt_DS_LMM = out_diffcyt_DS_LMM_supp_sample_sizes[[1]])

data_4vs4 <- list(diffcyt_DS_limma = out_diffcyt_DS_limma_supp_sample_sizes[[2]], 
                  diffcyt_DS_LMM = out_diffcyt_DS_LMM_supp_sample_sizes[[2]])

cobradata_2vs2 <- COBRAData(pval = data.frame(diffcyt_DS_limma = data_2vs2[["diffcyt_DS_limma"]][, "p_val"], 
                                              diffcyt_DS_LMM = data_2vs2[["diffcyt_DS_LMM"]][, "p_val"]), 
                            padj = data.frame(diffcyt_DS_limma = data_2vs2[["diffcyt_DS_limma"]][, "p_adj"], 
                                              diffcyt_DS_LMM = data_2vs2[["diffcyt_DS_LMM"]][, "p_adj"]), 
                            truth = data.frame(B_cell = data_2vs2[["diffcyt_DS_limma"]][, "B_cell"]))

cobradata_4vs4 <- COBRAData(pval = data.frame(diffcyt_DS_limma = data_4vs4[["diffcyt_DS_limma"]][, "p_val"], 
                                              diffcyt_DS_LMM = data_4vs4[["diffcyt_DS_LMM"]][, "p_val"]), 
                            padj = data.frame(diffcyt_DS_limma = data_4vs4[["diffcyt_DS_limma"]][, "p_adj"], 
                                              diffcyt_DS_LMM = data_4vs4[["diffcyt_DS_LMM"]][, "p_adj"]), 
                            truth = data.frame(B_cell = data_4vs4[["diffcyt_DS_limma"]][, "B_cell"]))

# calculate performance scores

cobraperf_2vs2 <- calculate_performance(cobradata_2vs2, 
                                        binary_truth = "B_cell", 
                                        aspects = c("roc", "fdrtpr", "fdrtprcurve", "tpr", "fpr"))

cobraperf_4vs4 <- calculate_performance(cobradata_4vs4, 
                                        binary_truth = "B_cell", 
                                        aspects = c("roc", "fdrtpr", "fdrtprcurve", "tpr", "fpr"))

# color scheme

colors <- c("firebrick1", "darkviolet")

# prepare plotting objects

cobraplot_2vs2 <- prepare_data_for_plot(cobraperf_2vs2, 
                                        colorscheme = colors, 
                                        conditionalfill = FALSE)

cobraplot_4vs4 <- prepare_data_for_plot(cobraperf_4vs4, 
                                        colorscheme = colors, 
                                        conditionalfill = FALSE)



# -----------------------------------
# Generate plots for each sample size
# -----------------------------------

cobraplot_list <- list(size_2vs2 = cobraplot_2vs2, 
                       size_4vs4 = cobraplot_4vs4)


for (i in 1:length(cobraplot_list)) {
  
  
  # select 'cobraplot' object from correct sample size
  cobraplot <- cobraplot_list[[i]]
  
  
  # ----------
  # ROC curves
  # ----------
  
  # create plot
  p_ROC <- 
    plot_roc(cobraplot, linewidth = 0.75) + 
    coord_fixed() + 
    xlab("False positive rate") + 
    ylab("True positive rate") + 
    ggtitle(paste0("BCR-XL-sim: smaller sample sizes: ", 
                   gsub("vs", " vs. ", gsub("^size_", "", names(cobraplot_list)[i]))), 
            subtitle = "ROC curve") + 
    theme_bw() + 
    theme(strip.text.x = element_blank()) + 
    guides(color = guide_legend("method"))
  
  # save plot
  fn <- file.path(DIR_PLOTS, "panels", 
                  paste0("results_BCR_XL_sim_diffcyt_supp_sample_sizes_ROC_", names(cobraplot_list)[i], ".pdf"))
  ggsave(fn, width = 4.25, height = 3)
  
  
  
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
    ggtitle(paste0("BCR-XL-sim: smaller sample sizes: ", 
                   gsub("vs", " vs. ", gsub("^size_", "", names(cobraplot_list)[i]))), 
            subtitle = "TPR vs. FDR") + 
    theme_bw() + 
    theme(strip.text.x = element_blank()) + 
    guides(shape = guide_legend("FDR threshold", override.aes = list(size = 4), order = 1), 
           color = guide_legend("method", order = 2))
  
  # save plot
  fn <- file.path(DIR_PLOTS, "panels", 
                  paste0("results_BCR_XL_sim_diffcyt_supp_sample_sizes_TPRFDR_", names(cobraplot_list)[i], ".pdf"))
  ggsave(fn, width = 4.25, height = 3)
  
  
  
  # ---------
  # TPR plots
  # ---------
  
  # create plot
  p_TPR <- 
    plot_tpr(cobraplot, pointsize = 4) + 
    scale_shape_manual(values = c(15, 19, 17), labels = c(0.01, 0.05, 0.1)) + 
    #coord_fixed() + 
    xlab("True positive rate") + 
    ggtitle(paste0("BCR-XL-sim: smaller sample sizes: ", 
                   gsub("vs", " vs. ", gsub("^size_", "", names(cobraplot_list)[i]))), 
            subtitle = "TPR") + 
    theme_bw() + 
    theme(strip.text.x = element_blank(), 
          axis.text.y = element_blank()) + 
    guides(shape = guide_legend("FDR threshold", override.aes = list(size = 4), order = 1), 
           color = guide_legend("method", override.aes = list(shape = 19, size = 4), order = 2))
  
  # save plot
  fn <- file.path(DIR_PLOTS, "panels", 
                  paste0("results_BCR_XL_sim_diffcyt_supp_sample_sizes_TPR_", names(cobraplot_list)[i], ".pdf"))
  ggsave(fn, width = 4.5, height = 3.5)
  
  
  
  # ---------
  # FPR plots
  # ---------
  
  # create plot
  p_FPR <- 
    plot_fpr(cobraplot, pointsize = 4) + 
    scale_shape_manual(values = c(15, 19, 17), labels = c(0.01, 0.05, 0.1)) + 
    #coord_fixed() + 
    xlab("False positive rate") + 
    ggtitle(paste0("BCR-XL-sim: smaller sample sizes: ", 
                   gsub("vs", " vs. ", gsub("^size_", "", names(cobraplot_list)[i]))), 
            subtitle = "FPR") + 
    theme_bw() + 
    theme(strip.text.x = element_blank(), 
          axis.text.y = element_blank()) + 
    guides(shape = guide_legend("FDR threshold", override.aes = list(size = 4), order = 1), 
           color = guide_legend("method", override.aes = list(shape = 19, size = 4), order = 2))
  
  # save plot
  fn <- file.path(DIR_PLOTS, "panels", 
                  paste0("results_BCR_XL_sim_diffcyt_supp_sample_sizes_FPR_", names(cobraplot_list)[i], ".pdf"))
  ggsave(fn, width = 4.5, height = 3.5)
  
  
  
  
  ###################
  # Multi-panel plots
  ###################
  
  # ----------
  # All panels
  # ----------
  
  plots_list <- list(p_ROC, p_TPRFDR, p_TPR, p_FPR)
  
  # modify plot elements
  plots_list <- lapply(plots_list, function(p) {
    p + 
      labs(title = p$labels$subtitle, subtitle = element_blank()) + 
      theme(legend.position = "none")
  })
  
  plots_multi <- plot_grid(plotlist = plots_list, 
                           nrow = 1, ncol = 4, align = "hv", axis = "bl")
  
  # add combined title
  title_single <- p_ROC$labels$title
  plots_title <- ggdraw() + draw_label(title_single)
  plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 7))
  
  # add combined legend
  legend_single <- get_legend(plots_list[[2]] + theme(legend.position = "right"))
  plots_multi <- plot_grid(plots_multi, legend_single, nrow = 1, rel_widths = c(6, 1))
  
  # save multi-panel plot
  fn <- file.path(DIR_PLOTS, 
                  paste0("results_BCR_XL_sim_diffcyt_supp_sample_sizes_", 
                         gsub("^size_", "", names(cobraplot_list)[i]), ".pdf"))
  ggsave(fn, width = 10, height = 2.625)
  
  
  # ----------------------
  # ROC and TPR-FDR curves
  # ----------------------
  
  plots_list <- list(p_ROC, p_TPRFDR)
  
  # modify plot elements
  plots_list <- lapply(plots_list, function(p) {
    p + 
      labs(title = p$labels$subtitle, subtitle = element_blank()) + 
      theme(legend.position = "none")
  })
  
  plots_multi <- plot_grid(plotlist = plots_list, 
                           nrow = 1, ncol = 2, align = "hv", axis = "bl")
  
  # add combined title
  title_single <- p_ROC$labels$title
  plots_title <- ggdraw() + draw_label(title_single)
  plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 7))
  
  # add combined legend
  legend_single <- get_legend(plots_list[[2]] + theme(legend.position = "right"))
  plots_multi <- plot_grid(plots_multi, legend_single, nrow = 1, rel_widths = c(3.25, 1))
  
  # save multi-panel plot
  fn <- file.path(DIR_PLOTS, 
                  paste0("results_BCR_XL_sim_diffcyt_supp_sample_sizes_", 
                         gsub("^size_", "", names(cobraplot_list)[i]), "_2_panels.pdf"))
  ggsave(fn, width = 6, height = 2.625)
  
}



###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



