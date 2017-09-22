##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: runtime
# - method: all methods
# 
# - main results
# 
# Lukas Weber, September 2017
##########################################################################################


# note: all methods except CellCnn use lineage markers only; CellCnn uses all markers


library(ggplot2)
library(reshape2)
library(cowplot)


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/main"

load(file.path(DIR_RDATA, "outputs_AML_sim_CellCnn_main.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_Citrus_main.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_cydar_main.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_GLMM_main.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_limma_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/main"




################
# Generate plots
################

# one panel per condition (j)


# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")


# store plots in list
plots_runtime <- vector("list", length(cond_names))


for (j in 1:length(cond_names)) {
  
  runtime_list <- vector("list", length(thresholds))
  names(runtime_list) <- thresholds
  
  for (th in 1:length(thresholds)) {
    runtime_list[[th]] <- c(CellCnn = runtime_CellCnn_main[[th]][[j]], 
                            Citrus = runtime_Citrus_main[[th]][[j]], 
                            cydar = runtime_cydar_main[[th]][[j]], 
                            diffcyt_DA_edgeR = runtime_diffcyt_DA_edgeR_main[[th]][[j]], 
                            diffcyt_DA_GLMM = runtime_diffcyt_DA_GLMM_main[[th]][[j]], 
                            diffcyt_DA_limma = runtime_diffcyt_DA_limma_main[[th]][[j]])
  }
  
  # create data frame for plotting
  d_plot <- as.data.frame(runtime_list)
  colnames(d_plot) <- thresholds
  d_plot$method <- factor(rownames(d_plot), levels = rownames(d_plot))
  
  d_plot <- melt(d_plot, id.vars = "method", variable.name = "threshold", value.name = "runtime")
  
  # color scheme
  colors <- c("mediumorchid3", "gold", "salmon", "darkblue", "deepskyblue2", "darkslategray2")
  
  # shapes
  shapes <- c(21, 24, 22, 23)
  names(shapes) <- thresholds
  
  # axis ranges
  y_range <- c(0, 425)
  
  # create plot
  p <- ggplot(d_plot, aes(x = method, y = runtime, group = method, color = method, shape = threshold)) + 
    geom_point(size = 2.25, stroke = 1.125) + 
    scale_color_manual(values = colors) + 
    #scale_fill_manual(values = colors) + 
    scale_shape_manual(values = shapes) + 
    #scale_shape_manual(values = shapes, guide = guide_legend(override.aes = list(fill = "black"))) + 
    #scale_y_log10() + 
    ylim(y_range) + 
    ylab("runtime (s)") + 
    ggtitle(cond_names[j]) + 
    theme_bw() + 
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))
    #guides(color = FALSE, fill = FALSE)
  
  plots_runtime[[j]] <- p
  
  # save individual panel plot
  p <- p + 
    ggtitle(paste0("AML-sim, main results, ", cond_names[j], ": runtime"))
  
  fn <- file.path(DIR_PLOTS, "panels", 
                  paste0("results_all_methods_main_runtime_", cond_names[j], ".pdf"))
  ggsave(fn, width = 7.5, height = 5)
  
}




########################
# Save multi-panel plots
########################

# remove duplicated annotation
plots_runtime <- lapply(plots_runtime, function(p) {
  p + theme(legend.position = "none", 
            axis.title.x = element_blank(), 
            axis.title.y = element_blank())
})

# format into grid
grid_runtime <- do.call(plot_grid, append(plots_runtime, list(labels = "AUTO", nrow = 1, ncol = 2, 
                                                              scale = 0.95, label_y = 0.975)))

# add combined axis titles
#xaxis_runtime <- ggdraw() + draw_label("", size = 12)
yaxis_runtime <- ggdraw() + draw_label("runtime (s)", size = 12, angle = 90, hjust = -0.1)

#grid_runtime <- plot_grid(grid_runtime, xaxis_runtime, ncol = 1, rel_heights = c(15, 1))
grid_runtime <- plot_grid(yaxis_runtime, grid_runtime, nrow = 1, rel_widths = c(1, 30))

# add combined legend
legend_runtime <- get_legend(plots_runtime[[1]] + theme(legend.position = "right"))
grid_legend <- plot_grid(legend_runtime, nrow = 2, rel_heights = c(10, 2.3))
grid_runtime <- plot_grid(grid_runtime, grid_legend, nrow = 1, rel_widths = c(10, 1.75))

# add combined title
title_runtime <- ggdraw() + draw_label("AML-sim, main results: runtime", fontface = "bold")
grid_runtime <- plot_grid(title_runtime, grid_runtime, ncol = 1, rel_heights = c(1, 16))

# save plots
fn_runtime <- file.path(DIR_PLOTS, "results_all_methods_main_runtime.pdf")
ggsave(fn_runtime, grid_runtime, width = 10, height = 5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



