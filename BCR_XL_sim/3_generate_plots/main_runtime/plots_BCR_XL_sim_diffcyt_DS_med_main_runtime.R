##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: runtime
# - method: diffcyt-DS-med
# 
# - main results
# 
# Lukas Weber, November 2017
##########################################################################################


library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL_sim/main"

load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_med_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/main_runtime"




################
# Generate plots
################

runtime_list <- list(diffcyt_DS_med = runtime_diffcyt_DS_med_main)

# create data frame for plotting
d_plot <- as.data.frame(t(as.data.frame(runtime_list)))
colnames(d_plot) <- "runtime"
d_plot$method <- factor(rownames(d_plot), levels = rownames(d_plot))

# color scheme
colors <- "darkblue"

# shapes
shapes <- 4

# axis limits
y_max <- 20

# create plot
ggplot(d_plot, aes(x = method, y = runtime, group = method, color = method, shape = method)) + 
  geom_point(size = 2.25, stroke = 1.125) + 
  scale_color_manual(values = colors) + 
  scale_shape_manual(values = shapes) + 
  ylim(c(0, y_max)) + 
  ylab("runtime (s)") + 
  ggtitle("BCR-XL-sim: main results", subtitle = "Runtime") + 
  theme_bw() + 
  theme(axis.title.x = element_blank())#, 
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_med_main_runtime.pdf")
ggsave(fn, width = 5, height = 3.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



