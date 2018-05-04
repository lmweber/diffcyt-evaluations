##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: performance metrics
# - method: diffcyt-DS-limma
# 
# - supplementary results: 'less distinct' data sets
# 
# Lukas Weber, May 2018
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2


# load saved results
DIR_RDATA_MAIN <- "../../../../RData/BCR_XL_sim/main"
DIR_RDATA_SUPP_LESS_DISTINCT <- "../../../../RData/BCR_XL_sim/supp_less_distinct"

load(file.path(DIR_RDATA_MAIN, "outputs_BCR_XL_sim_diffcyt_DS_limma_main.RData"))
load(file.path(DIR_RDATA_SUPP_LESS_DISTINCT, "outputs_BCR_XL_sim_diffcyt_DS_limma_supp_less_distinct.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/supp_less_distinct"




################
# Generate plots
################

# -------------------------------------
# Pre-processing steps for iCOBRA plots
# -------------------------------------

# note: create separate objects for each level of 'distinctness'


# create 'COBRAData' objects

data_main      <- out_diffcyt_DS_limma_main
data_less_50pc <- out_diffcyt_DS_limma_supp_less_distinct[["less_50pc"]]
data_less_75pc <- out_diffcyt_DS_limma_supp_less_distinct[["less_75pc"]]

cobradata_main      <- COBRAData(pval = data.frame(main = data_main[, "p_val"]), 
                                 padj = data.frame(main = data_main[, "p_adj"]), 
                                 truth = data.frame(B_cell = data_main[, "B_cell"]))

cobradata_less_50pc <- COBRAData(pval = data.frame(less_50pc = data_less_50pc[, "p_val"]), 
                                 padj = data.frame(less_50pc = data_less_50pc[, "p_adj"]), 
                                 truth = data.frame(B_cell = data_less_50pc[, "B_cell"]))

cobradata_less_75pc <- COBRAData(pval = data.frame(less_75pc = data_less_75pc[, "p_val"]), 
                                 padj = data.frame(less_75pc = data_less_75pc[, "p_adj"]), 
                                 truth = data.frame(B_cell = data_less_75pc[, "B_cell"]))

# calculate performance scores

cobraperf_main <- calculate_performance(cobradata_main, binary_truth = "B_cell", aspects = "roc")
cobraperf_less_50pc <- calculate_performance(cobradata_less_50pc, binary_truth = "B_cell", aspects = "roc")
cobraperf_less_75pc <- calculate_performance(cobradata_less_75pc, binary_truth = "B_cell", aspects = "roc")

# color scheme

colors <- "firebrick1"

# prepare plotting objects

cobraplot_main      <- prepare_data_for_plot(cobraperf_main, colorscheme = colors[1])
cobraplot_less_50pc <- prepare_data_for_plot(cobraperf_less_50pc, colorscheme = colors[1])
cobraplot_less_75pc <- prepare_data_for_plot(cobraperf_less_75pc, colorscheme = colors[1])

# linetypes

linetypes <- c("solid", "dashed", "dotted")



# ----------
# ROC curves
# ----------

cobraplot_list <- list(
  main = cobraplot_main, 
  less_50pc = cobraplot_less_50pc, 
  less_75pc = cobraplot_less_75pc
)

p_ROC_list <- vector("list", length(cobraplot_list))

distinctness_names <- c("main", "less_50pc", "less_75pc")

for (i in 1:length(cobraplot_list)) {
  
  # create plot
  p_ROC_list[[i]] <- 
    plot_roc(cobraplot_list[[i]], linewidth = 0.75) + 
    scale_linetype_manual(values = linetypes[i], guide = FALSE) + 
    coord_fixed() + 
    xlab("False positive rate") + 
    ylab("True positive rate") + 
    ggtitle("BCR-XL-sim: 'less distinct' benchmark data", subtitle = "ROC curves") + 
    theme_bw() + 
    theme(strip.text.x = element_blank()) + 
    guides(color = guide_legend("distinctness"), linetype = guide_legend("distinctness"))
  
  # save plot
  fn <- file.path(DIR_PLOTS, "panels", paste0("results_BCR_XL_sim_diffcyt_DS_limma_supp_less_distinct_ROC_", distinctness_names[i], ".pdf"))
  ggsave(fn, width = 4.75, height = 3.5)
  
}




###############
# Combined plot
###############

# combined plot showing all levels of distinctness

d_combined <- rbind(
  p_ROC_list[[1]]$data, 
  p_ROC_list[[2]]$data, 
  p_ROC_list[[3]]$data
)

d_combined$distinctness <- factor(d_combined$method, levels = distinctness_names)
d_combined$method_names <- as.factor("diffcyt-DS-limma")

# create plot
p_combined <- 
  ggplot(d_combined, aes(x = FPR, y = TPR, linetype = distinctness, color = method_names)) + 
  geom_line(lwd = 0.75) + 
  scale_color_manual(values = colors) + 
  scale_linetype_manual(values = linetypes) + 
  coord_fixed() + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  xlab("False positive rate") + 
  ylab("True positive rate") + 
  ggtitle("BCR-XL-sim: 'less distinct' benchmark data", subtitle = "ROC curves") + 
  theme_bw() + 
  guides(linetype = guide_legend("distinctness", order = 1), 
         color = guide_legend("method", order = 2))

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_limma_supp_less_distinct_ROC.pdf")
ggsave(fn, width = 4.75, height = 3.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



