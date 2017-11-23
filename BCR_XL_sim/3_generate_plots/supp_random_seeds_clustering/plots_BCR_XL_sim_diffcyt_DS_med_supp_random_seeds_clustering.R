##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: performance metrics
# - method: diffcyt-DS-med
# 
# - supplementary results: varying random seeds for clustering
# 
# Lukas Weber, November 2017
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2


# load saved results
DIR_RDATA_SUPP_MAIN <- "../../../../RData/BCR_XL_sim/main"
DIR_RDATA_SUPP_RANDOM_SEEDS_CLUSTERING <- "../../../../RData/BCR_XL_sim/supp_random_seeds_clustering"

load(file.path(DIR_RDATA_SUPP_MAIN, "outputs_BCR_XL_sim_diffcyt_DS_med_main.RData"))
load(file.path(DIR_RDATA_SUPP_RANDOM_SEEDS_CLUSTERING, "outputs_BCR_XL_sim_diffcyt_DS_med_supp_random_seeds_clustering.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/supp_random_seeds_clustering"




################
# Generate plots
################

# -------------------------------------
# Pre-processing steps for iCOBRA plots
# -------------------------------------

# create 'COBRAData' object
data <- list(seed_main = out_diffcyt_DS_med_main, 
             seed_1 = out_diffcyt_DS_med_supp_random_seeds_clustering[[1]], 
             seed_2 = out_diffcyt_DS_med_supp_random_seeds_clustering[[2]], 
             seed_3 = out_diffcyt_DS_med_supp_random_seeds_clustering[[3]])

# check
stopifnot(all(sapply(data, function(d) all(d$B_cell == data[[1]]$B_cell))))

# note: provide all available values
# 'padj' is required for threshold points on TPR-FDR curves
# depending on availability, plotting functions use 'score', then 'pval', then 'padj'
cobradata <- COBRAData(pval = data.frame(seed_main = data[["seed_main"]][, "p_vals"], 
                                         seed_1 = data[["seed_1"]][, "p_vals"], 
                                         seed_2 = data[["seed_2"]][, "p_vals"], 
                                         seed_3 = data[["seed_3"]][, "p_vals"]), 
                       padj = data.frame(seed_main = data[["seed_main"]][, "p_adj"], 
                                         seed_1 = data[["seed_1"]][, "p_adj"], 
                                         seed_2 = data[["seed_2"]][, "p_adj"], 
                                         seed_3 = data[["seed_3"]][, "p_adj"]), 
                       truth = data.frame(B_cell = data[["seed_main"]][, "B_cell"]))

# calculate performance scores
# (note: can ignore warning messages when 'padj' not available)
cobraperf <- calculate_performance(cobradata, 
                                   binary_truth = "B_cell", 
                                   aspects = "roc")

# color scheme
colors <- c("darkblue", "darkblue", "darkblue", "darkblue")

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

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_med_supp_random_seeds_clustering_ROC.pdf")
ggsave(fn, width = 4.75, height = 3.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



