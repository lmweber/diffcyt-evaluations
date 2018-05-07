##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: runtimes
# - method: all methods
# 
# - main results
# 
# Lukas Weber, May 2018
##########################################################################################


library(ggplot2)
library(cowplot)


# load saved results
DIR_RDATA_MAIN <- "../../../../RData/BCR_XL_sim/main"
DIR_RDATA_CITRUS <- "../../../../RData/BCR_XL_sim/comparisons_Citrus"
DIR_RDATA_CELLCNN <- "../../../../RData/BCR_XL_sim/comparisons_CellCnn"
DIR_RDATA_CYDAR <- "../../../../RData/BCR_XL_sim/comparisons_cydar"

load(file.path(DIR_RDATA_MAIN, "outputs_BCR_XL_sim_diffcyt_DS_limma_main.RData"))
load(file.path(DIR_RDATA_MAIN, "outputs_BCR_XL_sim_diffcyt_DS_LMM_main.RData"))
load(file.path(DIR_RDATA_CITRUS, "outputs_BCR_XL_sim_Citrus_main.RData"))
load(file.path(DIR_RDATA_CELLCNN, "outputs_BCR_XL_sim_CellCnn_main.RData"))
load(file.path(DIR_RDATA_CYDAR, "outputs_BCR_XL_sim_cydar_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/comparisons_runtimes"




################
# Generate plots
################

# create data frame for plotting
d_runtimes <- as.data.frame(c(
  Citrus = runtime_Citrus_main, 
  CellCnn = runtime_CellCnn_main, 
  cydar = runtime_cydar_main, 
  diffcyt_DS_limma = runtime_diffcyt_DS_limma_main, 
  diffcyt_DS_LMM = runtime_diffcyt_DS_LMM_main
))

colnames(d_runtimes) <- "runtime"

d_runtimes$method <- factor(rownames(d_runtimes), levels = rownames(d_runtimes))

# color scheme
colors <- c("gold2", "darkorange1", "brown4", "firebrick1", "darkviolet")

y_range <- c(10, 3000)

# create plot
p_runtimes <- 
  ggplot(d_runtimes, aes(x = method, y = runtime, color = method, label = sprintf("%.1f", runtime))) + 
  geom_point(shape = 4, size = 1.75, stroke = 1.5) + 
  geom_text(color = "black", vjust = -1.5, size = 3.4) + 
  scale_color_manual(values = colors) + 
  scale_y_log10(limits = y_range) + 
  ylab("runtime (sec, log10 scale)") + 
  ggtitle("BCR-XL-sim: runtimes") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))

p_runtimes

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_comparisons_runtimes.pdf")
ggsave(fn, width = 5.5, height = 4.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



