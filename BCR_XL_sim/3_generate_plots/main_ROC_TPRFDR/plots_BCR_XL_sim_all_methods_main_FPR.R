##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: FPR plots
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
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/main_ROC_TPR_FDR"




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
                                   aspects = "fpr")

# color scheme
#colors <- c("darkblue", "deepskyblue2", "darkslategray2")
colors <- c("darkblue")

colors <- colors[1:length(data)]
names(colors) <- names(data)

# prepare plotting object
cobraplot <- prepare_data_for_plot(cobraperf, 
                                   colorscheme = colors, 
                                   conditionalfill = FALSE)

# re-order legend
cobraplot <- reorder_levels(cobraplot, levels = names(data))



# ---------
# FPR plots
# ---------

# create plot
plot_fpr(cobraplot, pointsize = 5) + 
  scale_shape_manual(values = c(15, 19, 17), labels = c(0.01, 0.05, 0.1)) + 
  xlab("False positive rate") + 
  ggtitle(paste0("BCR-XL-sim, main results: FPR")) + 
  theme_bw() + 
  theme(strip.text.x = element_blank()) + 
  guides(shape = guide_legend("FDR threshold", override.aes = list(size = 4), order = 1), 
         color = guide_legend("method", override.aes = list(shape = 19, size = 4), order = 2))

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_all_methods_main_FPR.pdf")
ggsave(fn, width = 6.5, height = 4.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



