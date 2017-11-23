##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: performance metrics
# - method: diffcyt-DS-med
# 
# - supplementary results: varying random seeds for generating benchmark data
# 
# Lukas Weber, November 2017
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2


# load saved results
DIR_RDATA_SUPP_MAIN <- "../../../../RData/BCR_XL_sim/main"
DIR_RDATA_SUPP_RANDOM_SEEDS_DATA <- "../../../../RData/BCR_XL_sim/supp_random_seeds_data"

load(file.path(DIR_RDATA_SUPP_MAIN, "outputs_BCR_XL_sim_diffcyt_DS_med_main.RData"))
load(file.path(DIR_RDATA_SUPP_RANDOM_SEEDS_DATA, "outputs_BCR_XL_sim_diffcyt_DS_med_supp_random_seeds_data.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/supp_random_seeds_data"




################
# Generate plots
################

# -------------------------------------
# Pre-processing steps for iCOBRA plots
# -------------------------------------

# note: need to create separate objects for each seed, since each seed contains a
# different number of cells, and cells are not in the same order


# create 'COBRAData' objects

data_main   <- out_diffcyt_DS_med_main
data_seed_1 <- out_diffcyt_DS_med_supp_random_seeds_data[[1]]
data_seed_2 <- out_diffcyt_DS_med_supp_random_seeds_data[[2]]
data_seed_3 <- out_diffcyt_DS_med_supp_random_seeds_data[[3]]

cobradata_main   <- COBRAData(pval = data.frame(main = data_main[, "p_vals"]), 
                              padj = data.frame(main = data_main[, "p_adj"]), 
                              truth = data.frame(B_cell = data_main[, "B_cell"]))

cobradata_seed_1 <- COBRAData(pval = data.frame(seed_1 = data_seed_1[, "p_vals"]), 
                              padj = data.frame(seed_1 = data_seed_1[, "p_adj"]), 
                              truth = data.frame(B_cell = data_seed_1[, "B_cell"]))

cobradata_seed_2 <- COBRAData(pval = data.frame(seed_2 = data_seed_2[, "p_vals"]), 
                              padj = data.frame(seed_2 = data_seed_2[, "p_adj"]), 
                              truth = data.frame(B_cell = data_seed_2[, "B_cell"]))

cobradata_seed_3 <- COBRAData(pval = data.frame(seed_3 = data_seed_3[, "p_vals"]), 
                              padj = data.frame(seed_3 = data_seed_3[, "p_adj"]), 
                              truth = data.frame(B_cell = data_seed_3[, "B_cell"]))

# calculate performance scores
cobraperf_main   <- calculate_performance(cobradata_main, binary_truth = "B_cell", aspects = "roc")
cobraperf_seed_1 <- calculate_performance(cobradata_seed_1, binary_truth = "B_cell", aspects = "roc")
cobraperf_seed_2 <- calculate_performance(cobradata_seed_2, binary_truth = "B_cell", aspects = "roc")
cobraperf_seed_3 <- calculate_performance(cobradata_seed_3, binary_truth = "B_cell", aspects = "roc")

# color scheme
colors <- "darkblue"

# prepare plotting objects
cobraplot_main   <- prepare_data_for_plot(cobraperf_main, colorscheme = colors)
cobraplot_seed_1 <- prepare_data_for_plot(cobraperf_seed_1, colorscheme = colors)
cobraplot_seed_2 <- prepare_data_for_plot(cobraperf_seed_2, colorscheme = colors)
cobraplot_seed_3 <- prepare_data_for_plot(cobraperf_seed_3, colorscheme = colors)

# linetypes
linetypes <- c("solid", "dashed", "dashed", "dashed")



# ----------
# ROC curves
# ----------

cobraplot_list <- list(main = cobraplot_main, 
                       seed_1 = cobraplot_seed_1, 
                       seed_2 = cobraplot_seed_2, 
                       seed_3 = cobraplot_seed_3)

p_ROC_list <- vector("list", length(cobraplot_list))

seed_names <- c("main", "seed_1", "seed_2", "seed_3")


for (i in 1:length(cobraplot_list)) {
  
  # create plot
  p_ROC_list[[i]] <- 
    plot_roc(cobraplot_list[[i]], linewidth = 0.75) + 
    scale_linetype_manual(values = linetypes[i], guide = FALSE) + 
    coord_fixed() + 
    xlab("False positive rate") + 
    ylab("True positive rate") + 
    ggtitle("BCR-XL-sim: random seeds for data generation", subtitle = "ROC curve") + 
    theme_bw() + 
    theme(strip.text.x = element_blank()) + 
    guides(color = guide_legend("random seed"), linetype = guide_legend("random seed"))
  
  # save plot
  fn <- file.path(DIR_PLOTS, "panels", paste0("results_BCR_XL_sim_diffcyt_DS_med_supp_random_seeds_clustering_ROC_", seed_names[i], ".pdf"))
  ggsave(fn, width = 4.75, height = 3.5)
  
}



###############
# Combined plot
###############

# combined plot showing all seeds

d_combined <- rbind(p_ROC_list[[1]]$data, 
                    p_ROC_list[[2]]$data, 
                    p_ROC_list[[3]]$data, 
                    p_ROC_list[[4]]$data)

d_combined$method <- factor(d_combined$method, levels = seed_names)

# create plot
p_combined <- 
  ggplot(d_combined, aes(x = FPR, y = TPR, linetype = method)) + 
  geom_line(color = colors, lwd = 0.75) + 
  scale_linetype_manual(values = linetypes) + 
  coord_fixed() + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  xlab("False positive rate") + 
  ylab("True positive rate") + 
  ggtitle("BCR-XL-sim: random seeds for data generation", subtitle = "ROC curves") + 
  theme_bw() + 
  guides(color = guide_legend("random seed"), linetype = guide_legend("random seed"))

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_med_supp_random_seeds_data_ROC.pdf")
ggsave(fn, width = 4.75, height = 3.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



