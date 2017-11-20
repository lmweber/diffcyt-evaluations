##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: ROC curves
# - method: all methods
# 
# - main results
# 
# Lukas Weber, November 2017
##########################################################################################


library(iCOBRA)
library(ggplot2)


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL_sim/main"

load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_med_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/main_ROC_TPRFDR"




################
# Generate plots
################

# -------------------------------------
# Pre-processing steps for iCOBRA plots
# -------------------------------------

# create 'COBRAData' object
data <- list(diffcyt_DS_med = out_diffcyt_DS_med_main)

# check
stopifnot(all(sapply(data, function(d) all(d$B_cell == data[[1]]$B_cell))))

# note: provide all available values
# 'padj' is required for threshold points on TPR-FDR curves
# depending on availability, plotting functions use 'score', then 'pval', then 'padj'
cobradata <- COBRAData(pval = data.frame(diffcyt_DS_med = data[["diffcyt_DS_med"]][, "p_vals"]), 
                       padj = data.frame(diffcyt_DS_med = data[["diffcyt_DS_med"]][, "p_adj"]), 
                       truth = data.frame(B_cell = data[["diffcyt_DS_med"]][, "B_cell"]))

# calculate performance scores
# (note: can ignore warning messages when 'padj' not available)
cobraperf <- calculate_performance(cobradata, 
                                   binary_truth = "B_cell", 
                                   aspects = "roc")

# color scheme
#colors <- c("mediumorchid3", "gold", "salmon", "darkblue", "deepskyblue2", "darkslategray2")
colors <- c("darkblue")

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
plot_roc(cobraplot, linewidth = 0.75) + 
  coord_fixed() + 
  xlab("False positive rate") + 
  ylab("True positive rate") + 
  ggtitle(paste0("BCR-XL-sim, main results: ROC curves")) + 
  theme_bw() + 
  theme(strip.text.x = element_blank()) + 
  guides(color = guide_legend("method"))

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_all_methods_main_ROC_curves.pdf")
ggsave(fn, width = 6, height = 5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



