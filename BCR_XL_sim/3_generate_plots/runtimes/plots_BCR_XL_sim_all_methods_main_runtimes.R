##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: runtimes
# - method: all methods
# 
# - main results
# 
# Lukas Weber, November 2017
##########################################################################################


library(ggplot2)
library(cowplot)


# load saved results
DIR_RDATA_MAIN <- "../../../../RData/BCR_XL_sim/main"
DIR_RDATA_CYDAR <- "../../../../RData/BCR_XL_sim/comparisons_cydar"
DIR_RDATA_CELLCNN <- "../../../../RData/BCR_XL_sim/comparisons_CellCnn"
DIR_RDATA_CITRUS <- "../../../../RData/BCR_XL_sim/comparisons_Citrus"

load(file.path(DIR_RDATA_MAIN, "outputs_BCR_XL_sim_diffcyt_DS_med_main.RData"))
load(file.path(DIR_RDATA_CYDAR, "outputs_BCR_XL_sim_cydar_main.RData"))
load(file.path(DIR_RDATA_CELLCNN, "outputs_BCR_XL_sim_CellCnn_main.RData"))
load(file.path(DIR_RDATA_CITRUS, "outputs_BCR_XL_sim_Citrus_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/runtimes"




################
# Generate plots
################

# create data frame for plotting
d_runtimes <- as.data.frame(c(
  CellCnn = runtime_CellCnn_main, 
  Citrus = runtime_Citrus_main, 
  cydar = runtime_cydar_main, 
  diffcyt_DS_med = runtime_diffcyt_DS_med_main
))

colnames(d_runtimes) <- "runtime"

d_runtimes$method <- factor(rownames(d_runtimes))

# color scheme
colors <- c("mediumorchid3", "gold", "salmon", "darkblue")

y_range <- c(0, 520)

# create plot
p_runtimes <- 
  ggplot(d_runtimes, aes(x = method, y = runtime, color = method)) + 
  geom_point(shape = 4, size = 1.75, stroke = 1.5) + 
  scale_color_manual(values = colors) + 
  ylim(y_range) + 
  ylab("runtime (s)") + 
  ggtitle("Runtimes: BCR-XL-sim") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_all_methods_main_runtimes.pdf")
ggsave(fn, width = 5, height = 4)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



