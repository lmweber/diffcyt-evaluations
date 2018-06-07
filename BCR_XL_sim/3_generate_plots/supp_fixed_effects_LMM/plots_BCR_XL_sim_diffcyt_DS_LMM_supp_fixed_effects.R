##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: performance metrics
# - method: diffcyt-DS-LMM
# 
# - supplementary results: using fixed effects instead of random effects for patient IDs
# 
# Lukas Weber, June 2018
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2


# load saved results
DIR_RDATA_MAIN <- "../../../../RData/BCR_XL_sim/main"
DIR_RDATA_SUPP_FIXED_EFFECTS <- "../../../../RData/BCR_XL_sim/supp_fixed_effects_LMM"

load(file.path(DIR_RDATA_MAIN, "outputs_BCR_XL_sim_diffcyt_DS_LMM_main.RData"))
load(file.path(DIR_RDATA_SUPP_FIXED_EFFECTS, "outputs_BCR_XL_sim_diffcyt_DS_LMM_supp_fixed_effects.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/supp_fixed_effects_LMM"




################
# Generate plots
################

# -------------------------------------
# Pre-processing steps for iCOBRA plots
# -------------------------------------

# create 'COBRAData' object
data <- list(main = out_diffcyt_DS_LMM_main, 
             fixed_effects = out_diffcyt_DS_LMM_supp_fixed_effects)

# check
stopifnot(all(sapply(data, function(d) all(d$B_cell == data[[1]]$B_cell))))

# note: provide all available values
# 'padj' is required for threshold points on TPR-FDR curves
# depending on availability, plotting functions use 'score', then 'pval', then 'padj'
cobradata <- COBRAData(pval = data.frame(main = data[["main"]][, "p_val"], 
                                         fixed_effects = data[["fixed_effects"]][, "p_val"]), 
                       padj = data.frame(main = data[["main"]][, "p_adj"], 
                                         fixed_effects = data[["fixed_effects"]][, "p_adj"]), 
                       truth = data.frame(B_cell = data[["main"]][, "B_cell"]))

# calculate performance scores
# (note: can ignore warning messages when 'padj' not available)
cobraperf <- calculate_performance(cobradata, 
                                   binary_truth = "B_cell", 
                                   aspects = c("roc", "fdrtpr", "fdrtprcurve", "tpr", "fpr"))

# color scheme
colors <- rep("darkviolet", 2)

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
  ggtitle("BCR-XL-sim: diffcyt-DS-LMM", subtitle = "ROC curve") + 
  theme_bw() + 
  theme(strip.text.x = element_blank()) + 
  guides(color = guide_legend("method"), linetype = guide_legend("method"))

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_LMM_supp_fixed_effects_ROC.pdf")
ggsave(fn, width = 4.25, height = 3)



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
  ggtitle("BCR-XL-sim: diffcyt-DS-LMM", subtitle = "TPR vs. FDR") + 
  theme_bw() + 
  theme(strip.text.x = element_blank()) + 
  guides(shape = guide_legend("FDR threshold", override.aes = list(size = 4), order = 1), 
         color = guide_legend("method", order = 2, override.aes = list(shape = NA)), 
         linetype = guide_legend("method", order = 2))

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_LMM_supp_fixed_effects_TPRFDR.pdf")
ggsave(fn, width = 4.25, height = 3)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



