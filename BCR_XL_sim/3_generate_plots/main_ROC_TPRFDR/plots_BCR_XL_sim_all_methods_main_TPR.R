##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: TPR plots
# - method: diffcyt methods
# 
# - main results
# 
# Lukas Weber, November 2017
##########################################################################################


# note: showing 'diffcyt' methods only


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
stopifnot(all(sapply(data, function(d) all(d$spikein == data[[1]]$spikein))))

# note: provide all available values
# 'padj' is required for threshold points on TPR-FDR curves
# depending on availability, plotting functions use 'score', then 'pval', then 'padj'
cobradata <- COBRAData(pval = data.frame(diffcyt_DS_med = data[["diffcyt_DS_med"]][, "p_vals"]), 
                       padj = data.frame(diffcyt_DS_med = data[["diffcyt_DS_med"]][, "p_adj"]), 
                       truth = data.frame(spikein = data[["diffcyt_DS_med"]][, "spikein"]))

# calculate performance scores
# (note: can ignore warning messages when 'padj' not available)
cobraperf <- calculate_performance(cobradata, 
                                   binary_truth = "spikein", 
                                   aspects = "tpr")

# color scheme
#colors <- c("darkblue", "deepskyblue2", "darkslategray2")
colors <- "darkblue"

colors <- colors[1:length(data)]
names(colors) <- names(data)

# prepare plotting object
cobraplot <- prepare_data_for_plot(cobraperf, 
                                   colorscheme = colors, 
                                   conditionalfill = FALSE)

# re-order legend
cobraplot <- reorder_levels(cobraplot, levels = names(data))



# ---------
# TPR plots
# ---------

# create plot
p <- 
  plot_tpr(cobraplot, pointsize = 5) + 
  scale_shape_manual(values = c(15, 19, 17), labels = c(0.01, 0.05, 0.1)) + 
  coord_fixed() + 
  xlab("True positive rate") + 
  ggtitle(paste0("BCR-XL-sim, main results: TPR")) + 
  theme_bw() + 
  theme(strip.text.x = element_blank()) + 
  guides(shape = guide_legend("FDR threshold", override.aes = list(size = 4), order = 1), 
         color = guide_legend("method", override.aes = list(shape = 19, size = 4), order = 2))

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_methods_main_TPR.pdf")
ggsave(fn, width = 6.5, height = 4.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



