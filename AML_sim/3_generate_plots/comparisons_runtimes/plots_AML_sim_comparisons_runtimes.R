##########################################################################################
# Generate plots
# 
# - data set: AML-sim
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
DIR_RDATA_MAIN <- "../../../../RData/AML_sim/main"
DIR_RDATA_CITRUS <- "../../../../RData/AML_sim/comparisons_Citrus"
DIR_RDATA_CELLCNN <- "../../../../RData/AML_sim/comparisons_CellCnn"
DIR_RDATA_CYDAR <- "../../../../RData/AML_sim/comparisons_cydar"

load(file.path(DIR_RDATA_MAIN, "outputs_AML_sim_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA_MAIN, "outputs_AML_sim_diffcyt_DA_voom_main.RData"))
load(file.path(DIR_RDATA_MAIN, "outputs_AML_sim_diffcyt_DA_GLMM_main.RData"))
load(file.path(DIR_RDATA_CITRUS, "outputs_AML_sim_Citrus_main.RData"))
load(file.path(DIR_RDATA_CELLCNN, "outputs_AML_sim_CellCnn_main.RData"))
load(file.path(DIR_RDATA_CYDAR, "outputs_AML_sim_cydar_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/comparisons_runtimes"




################
# Generate plots
################

# loop over thresholds (th) and conditions (j)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# condition names
cond_names <- c("CN", "CBF")


# store plots in list
plots_all <- vector("list", length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  for (j in 1:length(cond_names)) {
    
    # index to store plots sequentially in list
    ix <- (j * length(thresholds)) - (length(thresholds) - th)
    
    
    # create data frame for plotting
    d_runtimes <- as.data.frame(c(
      Citrus = runtime_Citrus_main[[th]][[j]], 
      CellCnn = runtime_CellCnn_main[[th]][[j]], 
      cydar = runtime_cydar_main[[th]][[j]], 
      diffcyt_DA_edgeR = runtime_diffcyt_DA_edgeR_main[[th]][[j]], 
      diffcyt_DA_voom = runtime_diffcyt_DA_voom_main[[th]][[j]], 
      diffcyt_DA_GLMM = runtime_diffcyt_DA_GLMM_main[[th]][[j]]
    ))
    
    colnames(d_runtimes) <- "runtime"
    
    d_runtimes$method <- factor(rownames(d_runtimes), levels = rownames(d_runtimes))
    
    # color scheme
    colors <- c("gold2", "darkorange1", "brown4", "darkblue", "dodgerblue", "darkslategray3")
    
    y_range <- c(10, 3000)
    
    # create plot
    p_runtimes <- 
      ggplot(d_runtimes, aes(x = method, y = runtime, color = method, label = sprintf("%.1f", runtime))) + 
      geom_point(shape = 4, size = 1.75, stroke = 1.5) + 
      geom_text(color = "black", vjust = -1.5, size = 3.4) + 
      scale_color_manual(values = colors) + 
      scale_y_log10(limits = y_range) + 
      ylab("runtime (s)") + 
      ggtitle(paste0("AML-sim: runtimes, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th]))) + 
      theme_bw() + 
      theme(axis.title.x = element_blank(), 
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))
    
    plots_all[[ix]] <- p_runtimes
    
    # save individual panel plot
    fn <- file.path(DIR_PLOTS, "panels", 
                    paste0("results_AML_sim_comparisons_runtimes_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 5.5, height = 4.5)
    
  }
}




########################
# Save multi-panel plots
########################

# modify plot elements
plots_multi <- lapply(plots_all, function(p) {
  p + 
    labs(title = gsub("^.*runtimes, ", "", p$labels$title)) + 
    theme(legend.position = "none")
})

# format into grid
plots_multi <- plot_grid(plotlist = plots_multi, nrow = 2, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_multi <- ggdraw() + draw_label(gsub(",.*$", "", plots_all[[1]]$labels$title), fontface = "bold")
grid_multi <- plot_grid(title_multi, plots_multi, ncol = 1, rel_heights = c(1, 25))

# add combined legend
legend_multi <- get_legend(plots_all[[1]] + theme(legend.position = "right", 
                                                  legend.title = element_text(size = 12, face = "bold"), 
                                                  legend.text = element_text(size = 12)))
grid_multi <- plot_grid(grid_multi, legend_multi, nrow = 1, rel_widths = c(5, 1))

# save plots
fn_multi <- file.path(DIR_PLOTS, paste0("results_AML_sim_comparisons_runtimes.pdf"))
ggsave(fn_multi, grid_multi, width = 11, height = 8)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



