##########################################################################################
# Generate plots
# 
# - plot type: ROC and TPR-FDR curves
# - methods: all methods
# - data set: BCR-XL-sim
# 
# Lukas Weber, August 2017
##########################################################################################


library(iCOBRA)
library(ggplot2)


# load saved results
DIR_RDATA <- "../../../RData/BCR_XL_sim"

load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DE_med.RData"))


# path to save plots
DIR_PLOTS <- "../../../plots/BCR_XL_sim/all_methods"

DIR_TIMESTAMP <- "../../../plots/BCR_XL_sim"




#####################################
# Plots: test results for each method
#####################################

## to do: add other markers (not just 'pS6')
## to do: add other simulations (not just 'sim_full')


# -------------------------------------
# Pre-processing steps for iCOBRA plots
# -------------------------------------

# create 'COBRAData' object
data <- list(diffcyt_DE_med = out_diffcyt_DE_med[["pS6(Yb172)Dd"]])

# check
stopifnot(all(sapply(data, function(d) all(d$spikein == data[[1]]$spikein))))

# note: provide all available values
# - 'padj' is required for threshold points on TPR-FDR curves
# - depending on availability, plotting functions use 'score', then 'pval', then 'padj'
cobradata <- COBRAData(pval = data.frame(diffcyt_DE_med = data[["diffcyt_DE_med"]][, "p_vals"]), 
                       padj = data.frame(diffcyt_DE_med = data[["diffcyt_DE_med"]][, "p_adj"]), 
                       #score = data.frame(), 
                       truth = data.frame(spikein = data[["diffcyt_DE_med"]][, "spikein"]))

# calculate performance scores
# (note: can ignore warning messages when 'padj' not available)
cobraperf <- calculate_performance(cobradata, 
                                   binary_truth = "spikein", 
                                   aspects = c("fdrtpr", "fdrtprcurve", "roc"))

# color scheme
# modifed default "Set1" to use different yellow (#FFD92F) from colorbrewer2.org
colors <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFD92F', '#A65628', '#F781BF')

colors <- colors[1:length(data)]
names(colors) <- names(data)

# x axis labels
x_min <- 0
x_max <- 1
x_labels <- seq(x_min, x_max, by = 0.1)

# prepare plotting object
cobraplot <- prepare_data_for_plot(cobraperf, 
                                   colorscheme = colors, 
                                   conditionalfill = FALSE)


# ----------
# ROC curves
# ----------

# create plot
p <- plot_roc(cobraplot, linewidth = 0.75)

p + 
  coord_fixed() + 
  xlab("False positive rate") + 
  ylab("True positive rate") + 
  ggtitle("ROC curves: BCR-XL-sim, sim_full") + 
  theme_bw() + 
  theme(strip.text.x = element_blank())

path <- DIR_PLOTS
filename <- file.path(path, "results_BCR_XL_sim_all_methods_ROC_curves.pdf")

ggsave(filename, width = 9, height = 8)


# --------------
# TPR-FDR curves
# --------------

# create plot
p <- plot_fdrtprcurve(cobraplot, linewidth = 0.75, pointsize = 4)

p + 
  scale_shape_manual(values = c(22, 21, 23)) + 
  scale_x_continuous(breaks = x_labels, labels = x_labels) + 
  coord_fixed() + 
  xlab("False discovery rate") + 
  ylab("True positive rate") + 
  ggtitle("TPR-FDR curves: BCR-XL-sim, sim_full") + 
  theme_bw() + 
  theme(strip.text.x = element_blank()) + 
  guides(color = guide_legend(override.aes = list(shape = NA)), shape = FALSE)

path <- DIR_PLOTS
filename <- file.path(path, "results_BCR_XL_sim_all_methods_TPR_FDR_curves.pdf")

ggsave(filename, width = 9, height = 8)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_TIMESTAMP, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



