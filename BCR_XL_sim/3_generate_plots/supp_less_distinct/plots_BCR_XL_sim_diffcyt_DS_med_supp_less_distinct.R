##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: performance metrics
# - method: diffcyt-DS-med
# 
# - supplementary results: 'less distinct' data sets
# 
# Lukas Weber, November 2017
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2


# load saved results
DIR_RDATA_SUPP_MAIN <- "../../../../RData/BCR_XL_sim/main"
DIR_RDATA_SUPP_LESS_DISTINCT <- "../../../../RData/BCR_XL_sim/supp_less_distinct"

load(file.path(DIR_RDATA_SUPP_MAIN, "outputs_BCR_XL_sim_diffcyt_DS_med_main.RData"))
load(file.path(DIR_RDATA_SUPP_LESS_DISTINCT, "outputs_BCR_XL_sim_diffcyt_DS_med_supp_less_distinct.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/supp_less_distinct"




#########################################
# Generate plots: all performance metrics
#########################################

# -------------------------------------
# Pre-processing steps for iCOBRA plots
# -------------------------------------

# note: create separate objects for each level of 'distinctness'


# create 'COBRAData' objects
data_less_50pc <- list(diffcyt_DS_med = out_diffcyt_DS_med_supp_less_distinct[["less_50pc"]])
data_less_75pc <- list(diffcyt_DS_med = out_diffcyt_DS_med_supp_less_distinct[["less_75pc"]])

cobradata_less_50pc <- COBRAData(pval = data.frame(diffcyt_DS_med = data_less_50pc[["diffcyt_DS_med"]][, "p_vals"]), 
                                 padj = data.frame(diffcyt_DS_med = data_less_50pc[["diffcyt_DS_med"]][, "p_adj"]), 
                                 truth = data.frame(B_cell = data_less_50pc[["diffcyt_DS_med"]][, "B_cell"]))

cobradata_less_75pc <- COBRAData(pval = data.frame(diffcyt_DS_med = data_less_75pc[["diffcyt_DS_med"]][, "p_vals"]), 
                                 padj = data.frame(diffcyt_DS_med = data_less_75pc[["diffcyt_DS_med"]][, "p_adj"]), 
                                 truth = data.frame(B_cell = data_less_75pc[["diffcyt_DS_med"]][, "B_cell"]))

# calculate performance scores
cobraperf_less_50pc <- calculate_performance(cobradata_less_50pc, 
                                             binary_truth = "B_cell", 
                                             aspects = c("roc", "fdrtpr", "fdrtprcurve", "tpr", "fpr"))

cobraperf_less_75pc <- calculate_performance(cobradata_less_75pc, 
                                             binary_truth = "B_cell", 
                                             aspects = c("roc", "fdrtpr", "fdrtprcurve", "tpr", "fpr"))

# color scheme
colors <- c("darkblue")

# prepare plotting objects
cobraplot_less_50pc <- prepare_data_for_plot(cobraperf_less_50pc, 
                                             colorscheme = colors, 
                                             conditionalfill = FALSE)

cobraplot_less_75pc <- prepare_data_for_plot(cobraperf_less_75pc, 
                                             colorscheme = colors, 
                                             conditionalfill = FALSE)



# -----------------------------------------------
# Generate plots for each level of 'distinctness'
# -----------------------------------------------

cobraplot_list <- list(less_50pc = cobraplot_less_50pc, 
                       less_75pc = cobraplot_less_75pc)


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
    ggtitle(paste0("BCR-XL-sim: 'less distinct' data: reduced ", 
                   gsub("pc$", "\\%", gsub("^less_", "", names(cobraplot_list)[i]))), 
            subtitle = "ROC curve") + 
    theme_bw() + 
    theme(strip.text.x = element_blank()) + 
    guides(color = guide_legend("method"))
  
  # save plot
  fn <- file.path(DIR_PLOTS, "panels", 
                  paste0("results_BCR_XL_sim_diffcyt_DS_med_supp_less_distinct_ROC_", names(cobraplot_list)[i], ".pdf"))
  ggsave(fn, width = 4.75, height = 3.5)
  
  
  
  # --------------
  # TPR-FDR curves
  # --------------
  
  # create plot
  p_TPRFDR <- 
    plot_fdrtprcurve(cobraplot, linewidth = 0.75, pointsize = 4) + 
    #scale_shape_manual(values = c(22, 21, 24), labels = c(0.01, 0.05, 0.1)) + 
    scale_shape_manual(values = c(15, 19, 17), labels = c(0.01, 0.05, 0.1)) + 
    coord_fixed() + 
    xlab("False discovery rate") + 
    ylab("True positive rate") + 
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) + 
    ggtitle(paste0("BCR-XL-sim: 'less distinct' data: reduced ", 
                   gsub("pc$", "\\%", gsub("^less_", "", names(cobraplot_list)[i]))), 
            subtitle = "TPR vs. FDR") + 
    theme_bw() + 
    theme(strip.text.x = element_blank()) + 
    guides(shape = guide_legend("FDR threshold", override.aes = list(size = 4), order = 1), 
           color = guide_legend("method", order = 2))
  
  # save plot
  fn <- file.path(DIR_PLOTS, "panels", 
                  paste0("results_BCR_XL_sim_diffcyt_DS_med_supp_less_distinct_TPRFDR_", names(cobraplot_list)[i], ".pdf"))
  ggsave(fn, width = 4.75, height = 3.5)
  
  
  
  # ---------
  # TPR plots
  # ---------
  
  # create plot
  p_TPR <- 
    plot_tpr(cobraplot, pointsize = 4) + 
    #scale_shape_manual(values = c(22, 21, 24), labels = c(0.01, 0.05, 0.1)) + 
    scale_shape_manual(values = c(15, 19, 17), labels = c(0.01, 0.05, 0.1)) + 
    coord_fixed() + 
    xlab("True positive rate") + 
    ggtitle(paste0("BCR-XL-sim: 'less distinct' data: reduced ", 
                   gsub("pc$", "\\%", gsub("^less_", "", names(cobraplot_list)[i]))), 
            subtitle = "TPR") + 
    theme_bw() + 
    theme(strip.text.x = element_blank(), 
          axis.text.y = element_blank()) + 
    guides(shape = guide_legend("FDR threshold", override.aes = list(size = 4), order = 1), 
           color = guide_legend("method", override.aes = list(shape = 19, size = 4), order = 2))
  
  # save plot
  fn <- file.path(DIR_PLOTS, "panels", 
                  paste0("results_BCR_XL_sim_diffcyt_DS_med_supp_less_distinct_TPR_", names(cobraplot_list)[i], ".pdf"))
  ggsave(fn, width = 4.5, height = 3.5)
  
  
  
  # ---------
  # FPR plots
  # ---------
  
  # create plot
  p_FPR <- 
    plot_fpr(cobraplot, pointsize = 4) + 
    #scale_shape_manual(values = c(22, 21, 24), labels = c(0.01, 0.05, 0.1)) + 
    scale_shape_manual(values = c(15, 19, 17), labels = c(0.01, 0.05, 0.1)) + 
    coord_fixed() + 
    xlab("False positive rate") + 
    ggtitle(paste0("BCR-XL-sim: 'less distinct' data: reduced ", 
                   gsub("pc$", "\\%", gsub("^less_", "", names(cobraplot_list)[i]))), 
            subtitle = "FPR") + 
    theme_bw() + 
    theme(strip.text.x = element_blank(), 
          axis.text.y = element_blank()) + 
    guides(shape = guide_legend("FDR threshold", override.aes = list(size = 4), order = 1), 
           color = guide_legend("method", override.aes = list(shape = 19, size = 4), order = 2))
  
  # save plot
  fn <- file.path(DIR_PLOTS, "panels", 
                  paste0("results_BCR_XL_sim_diffcyt_DS_med_supp_less_distinct_FPR_", names(cobraplot_list)[i], ".pdf"))
  ggsave(fn, width = 4.5, height = 3.5)
  
  
  
  
  ##################
  # Multi-panel plot
  ##################
  
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
                  paste0("results_BCR_XL_sim_diffcyt_DS_med_supp_less_distinct_", 
                         gsub("^less_", "", names(cobraplot_list)[i]), ".pdf"))
  ggsave(fn, width = 10, height = 2.75)
  
}




#############################################
# Generate plots: combined plot of ROC curves
#############################################

# -------------------------------------
# Pre-processing steps for iCOBRA plots
# -------------------------------------

# create 'COBRAData' object
data <- list(main = out_diffcyt_DS_med_main, 
             less_50pc = out_diffcyt_DS_med_supp_less_distinct[[1]], 
             less_75pc = out_diffcyt_DS_med_supp_less_distinct[[2]])

# check
stopifnot(all(sapply(data, function(d) all(d$B_cell == data[[1]]$B_cell))))

# note: provide all available values
# 'padj' is required for threshold points on TPR-FDR curves
# depending on availability, plotting functions use 'score', then 'pval', then 'padj'
cobradata <- COBRAData(pval = data.frame(main = data[["main"]][, "p_vals"], 
                                         less_50pc = data[["less_50pc"]][, "p_vals"], 
                                         less_75pc = data[["less_75pc"]][, "p_vals"]), 
                       padj = data.frame(main = data[["main"]][, "p_adj"], 
                                         less_50pc = data[["less_50pc"]][, "p_adj"], 
                                         less_75pc = data[["less_75pc"]][, "p_adj"]), 
                       truth = data.frame(B_cell = data[["main"]][, "B_cell"]))

# calculate performance scores
# (note: can ignore warning messages when 'padj' not available)
cobraperf <- calculate_performance(cobradata, 
                                   binary_truth = "B_cell", 
                                   aspects = "roc")

# color scheme
colors <- c("darkblue", "darkblue", "darkblue")

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



# -------------------------
# ROC curves: combined plot
# -------------------------

# create plot
p_ROC <- 
  plot_roc(cobraplot, linewidth = 0.75) + 
  scale_linetype_manual(values = linetypes, guide = FALSE) + 
  coord_fixed() + 
  xlab("False positive rate") + 
  ylab("True positive rate") + 
  ggtitle("BCR-XL-sim: 'less distinct' benchmark data", subtitle = "ROC curves") + 
  theme_bw() + 
  theme(strip.text.x = element_blank()) + 
  guides(color = guide_legend("distinctness"), linetype = guide_legend("distinctness"))

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_med_supp_less_distinct_ROC_combined.pdf")
ggsave(fn, width = 4.75, height = 3.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



