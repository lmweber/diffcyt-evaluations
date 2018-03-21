##########################################################################################
# Generate plots
# 
# - data set: Anti-PD-1
# - plot type: runtimes
# - method: diffcyt methods
# 
# - main results
# 
# Lukas Weber, March 2018
##########################################################################################


library(ggplot2)
library(reshape2)


# load saved results
DIR_RDATA <- "../../../../RData/Anti_PD_1/main"

load(file.path(DIR_RDATA, "outputs_Anti_PD_1_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA, "outputs_Anti_PD_1_diffcyt_DA_GLMM_main.RData"))
load(file.path(DIR_RDATA, "outputs_Anti_PD_1_diffcyt_DA_limma_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/Anti_PD_1/main_runtimes"




################
# Generate plots
################

runtime <- c(
  diffcyt_DA_edgeR = runtime_diffcyt_DA_edgeR_main, 
  diffcyt_DA_GLMM = runtime_diffcyt_DA_GLMM_main, 
  diffcyt_DA_limma = runtime_diffcyt_DA_limma_main
)


# create data frame for plotting
d_plot <- as.data.frame(runtime)
d_plot$method <- factor(rownames(d_plot), levels = rownames(d_plot))

# color scheme
colors <- c("darkblue", "deepskyblue2", "darkslategray2")

# create plot
p <- 
  ggplot(d_plot, aes(x = method, y = runtime, group = method, color = method)) + 
  geom_point(shape = 4, size = 2.25, stroke = 1.125) + 
  scale_color_manual(values = colors) + 
  ylim(c(0, 100)) + 
  ylab("runtime (s)") + 
  ggtitle("Anti-PD-1, main results: runtimes") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))

# save plot
fn <- file.path(DIR_PLOTS, "results_Anti_PD_1_diffcyt_runtimes.pdf")
ggsave(file = fn, width = 6, height = 5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



