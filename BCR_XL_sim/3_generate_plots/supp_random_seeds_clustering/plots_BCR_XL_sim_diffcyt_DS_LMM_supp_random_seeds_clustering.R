##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: performance metrics
# - method: diffcyt-DS-LMM
# 
# - supplementary results: varying random seeds for clustering
# 
# Lukas Weber, May 2018
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2


# load saved results
DIR_RDATA_MAIN <- "../../../../RData/BCR_XL_sim/main"
DIR_RDATA_SUPP_RANDOM_SEEDS_CLUSTERING <- "../../../../RData/BCR_XL_sim/supp_random_seeds_clustering"

load(file.path(DIR_RDATA_MAIN, "outputs_BCR_XL_sim_diffcyt_DS_LMM_main.RData"))
load(file.path(DIR_RDATA_SUPP_RANDOM_SEEDS_CLUSTERING, "outputs_BCR_XL_sim_diffcyt_DS_LMM_supp_random_seeds_clustering.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/supp_random_seeds_clustering"




################
# Generate plots
################

# -------------------------------------
# Pre-processing steps for iCOBRA plots
# -------------------------------------

# create 'COBRAData' object
data <- list(main = out_diffcyt_DS_LMM_main, 
             seed_1 = out_diffcyt_DS_LMM_supp_random_seeds_clustering[[1]], 
             seed_2 = out_diffcyt_DS_LMM_supp_random_seeds_clustering[[2]], 
             seed_3 = out_diffcyt_DS_LMM_supp_random_seeds_clustering[[3]])

# check
stopifnot(all(sapply(data, function(d) all(d$B_cell == data[[1]]$B_cell))))

# note: provide all available values
# 'padj' is required for threshold points on TPR-FDR curves
# depending on availability, plotting functions use 'score', then 'pval', then 'padj'
cobradata <- COBRAData(pval = data.frame(main = data[["main"]][, "p_val"], 
                                         seed_1 = data[["seed_1"]][, "p_val"], 
                                         seed_2 = data[["seed_2"]][, "p_val"], 
                                         seed_3 = data[["seed_3"]][, "p_val"]), 
                       padj = data.frame(main = data[["main"]][, "p_adj"], 
                                         seed_1 = data[["seed_1"]][, "p_adj"], 
                                         seed_2 = data[["seed_2"]][, "p_adj"], 
                                         seed_3 = data[["seed_3"]][, "p_adj"]), 
                       truth = data.frame(B_cell = data[["main"]][, "B_cell"]))

# calculate performance scores
# (note: can ignore warning messages when 'padj' not available)
cobraperf <- calculate_performance(cobradata, 
                                   binary_truth = "B_cell", 
                                   aspects = "roc")

# color scheme
colors <- rep("darkviolet", length(data))

colors <- colors[1:length(data)]
names(colors) <- names(data)

# linetypes
linetypes <- c("solid", "dashed", "dashed", "dashed")
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
  ggtitle("BCR-XL-sim: random seeds for clustering", subtitle = "ROC curves") + 
  theme_bw() + 
  theme(strip.text.x = element_blank()) + 
  guides(color = guide_legend("random seed"), linetype = guide_legend("random seed"))


# modify plot to allow multiple legends
d_plot <- p_ROC$data
seed_names <- c("main", "seed_1", "seed_2", "seed_3")
d_plot$seed_names <- factor(d_plot$method, levels = seed_names)
d_plot$method_names <- as.factor("diffcyt-DS-LMM")

ggplot(d_plot, aes(x = FPR, y = TPR, linetype = seed_names, color = method_names)) + 
  geom_line(lwd = 0.75) + 
  scale_color_manual(values = unname(colors[1])) + 
  scale_linetype_manual(values = linetypes) + 
  coord_fixed() + 
  xlim(c(0, 1)) + 
  ylim(c(0, 1)) + 
  xlab("False positive rate") + 
  ylab("True positive rate") + 
  ggtitle("BCR-XL-sim: random seeds for clustering", subtitle = "ROC curves") + 
  theme_bw() + 
  guides(linetype = guide_legend("random seed", order = 1), 
         color = guide_legend("method", order = 2))


# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_LMM_supp_random_seeds_clustering_ROC.pdf")
ggsave(fn, width = 4.75, height = 3.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



