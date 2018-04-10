##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: clustering performance
# - method: diffcyt methods
# 
# - main results
# 
# Lukas Weber, April 2018
##########################################################################################


# note: clustering step is the same for all 'diffcyt' methods (diffcyt-DA-edgeR, diffcyt-DA-voom, diffcyt-DA-GLMM)


library(SummarizedExperiment)
library(reshape2)
library(ggplot2)
library(cowplot)


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/main"

# note: only need to load one set of results, since clustering step is the same for all 'diffcyt' methods
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA, "out_objects_AML_sim_diffcyt_DA_edgeR_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/main_clustering_performance"




##################################
# Calculate clustering performance
##################################

# loop over thresholds (th) and conditions (j)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# condition names
cond_names <- c("CN", "CBF")
cond_names_all <- c("healthy", cond_names)

# lists to store clustering performance results
clustering_pr <- clustering_re <- clustering_F1 <- labels <- vector("list", length(thresholds))
names(clustering_pr) <- names(clustering_re) <- names(clustering_F1) <- names(labels) <- thresholds


for (th in 1:length(thresholds)) {
  
  clustering_pr[[th]] <- clustering_re[[th]] <- clustering_F1[[th]] <- labels[[th]] <- vector("list", length(cond_names))
  names(clustering_pr[[th]]) <- names(clustering_re[[th]]) <- names(clustering_F1[[th]]) <- names(labels[[th]]) <- cond_names
  
  
  for (j in 1:length(cond_names)) {
    
    # ------------------------------------------
    # load data objects and true spike-in status
    # ------------------------------------------
    
    # load data objects
    # note: clustering is performed once on all samples from both conditions together
    d_se <- out_objects_diffcyt_DA_edgeR_main[[th]]$d_se
    
    # load spike-in status at cell level (for condition j)
    spikein <- out_diffcyt_DA_edgeR_main[[th]][[j]]$spikein
    
    # add spike-in status to data object (for condition j)
    rowData(d_se)$spikein <- 0
    rowData(d_se)$spikein[rowData(d_se)$group_id %in% c("healthy", cond_names[j])] <- spikein
    
    
    # --------------------------------------------------------------------------------
    # calculate clustering performance for all clusters containing true spike-in cells
    # --------------------------------------------------------------------------------
    
    # find matching clusters (clusters containing true spike-in cells)
    
    # check no missing clusters
    stopifnot(all(names(table(rowData(d_se)[rowData(d_se)$spikein == 1, ]$cluster_id)) == levels(rowData(d_se)$cluster_id)))
    
    labels_matched <- unname(which(table(rowData(d_se)[rowData(d_se)$spikein == 1, ]$cluster_id) > 0))
    labels_matched
    
    # total number of cells in each matching cluster
    n_matched <- sapply(labels_matched, function(l) sum(rowData(d_se)$cluster_id == l))
    n_matched
    
    # number of true spike-in cells in each matching cluster
    n_correct <- sapply(labels_matched, function(l) sum(rowData(d_se)$cluster_id == l & rowData(d_se)$spikein == 1))
    n_correct
    
    # total number of true spike-in cells
    n_spikein <- sum(rowData(d_se)$spikein == 1)
    n_spikein
    
    # calculate precision, recall, F1 score for each matching cluster
    
    stopifnot(length(n_matched) == length(n_correct), 
              length(n_spikein) == 1)
    
    pr <- n_correct / n_matched
    re <- n_correct / n_spikein
    F1 <- 2 * (pr * re) / (pr + re)
    
    # store results
    labels[[th]][[j]] <- labels_matched
    clustering_pr[[th]][[j]] <- pr
    clustering_re[[th]][[j]] <- re
    clustering_F1[[th]][[j]] <- F1
    
  }
}




################
# Generate plots
################

# loop over thresholds (th) and conditions (j)

# store plots in list
plots_clustering <- vector("list", length(thresholds) * length(cond_names))

plot_widths <- rep(NA, length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  for (j in 1:length(cond_names)) {
    
    # index to store plots sequentially in list
    ix <- (j * length(thresholds)) - (length(thresholds) - th)
    
    
    # create data frame for plotting
    
    d_plot <- data.frame(
      cluster = labels[[th]][[j]], 
      precision = clustering_pr[[th]][[j]], 
      recall = clustering_re[[th]][[j]], 
      F1_score = clustering_F1[[th]][[j]]
    )
    
    plot_widths[ix] <- 2 + nrow(d_plot) / 7
    
    # sort by F1 score
    d_plot <- d_plot[rev(order(d_plot$F1_score)), ]
    
    d_plot$cluster <- factor(d_plot$cluster, levels = as.character(d_plot$cluster))
    d_plot <- melt(d_plot, id.vars = "cluster", variable.name = "measure")
    d_plot$measure <- factor(d_plot$measure, levels = c("F1_score", "precision", "recall"))
    
    # create plot
    
    colors <- c("firebrick1", "forestgreen", "deepskyblue")
    
    p <- 
      ggplot(d_plot, aes(x = cluster, y = value, color = measure)) + 
      geom_point(shape = 1, stroke = 1) + 
      scale_color_manual(values = colors) + 
      ylim(c(-0.025, 1.025)) + 
      ggtitle(paste0(cond_names[j], ", threshold ", gsub("pc$", "\\%", thresholds[th]))) + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8), 
            axis.title.y = element_blank())
    
    plots_clustering[[ix]] <- p
    
    # save individual panel plot
    fn <- file.path(DIR_PLOTS, "panels", 
                    paste0("results_AML_sim_diffcyt_main_clustering_performance_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = plot_widths[ix], height = 3)
    
  }
}




########################
# Save multi-panel plots
########################

# modify plot elements
plots_clustering <- lapply(plots_clustering, function(p) {
  p + theme(legend.position = "none")
})

# format into grid
plot_widths_avg <- c(7, 3.75, 2)
grid_clustering <- do.call(plot_grid, append(plots_clustering, list(
  nrow = 2, ncol = 3, align = "hv", axis = "bl", rel_widths = plot_widths_avg))
)

# add combined title
title_clustering <- ggdraw() + draw_label("AML-sim, diffcyt methods: clustering performance", fontface = "bold")
grid_clustering <- plot_grid(title_clustering, grid_clustering, ncol = 1, rel_heights = c(1, 25))

# add combined legend (one legend per row)
legend_clustering <- get_legend(plots_clustering[[1]] + theme(legend.position = "right", 
                                                              legend.title = element_text(size = 11, face = "bold"), 
                                                              legend.text = element_text(size = 11)))
legend_clustering <- plot_grid(legend_clustering, legend_clustering, ncol = 1)
grid_clustering <- plot_grid(grid_clustering, legend_clustering, nrow = 1, rel_widths = c(9, 1))

# save plots
fn_clustering <- file.path(DIR_PLOTS, paste0("results_AML_sim_diffcyt_main_clustering_performance.pdf"))
ggsave(fn_clustering, grid_clustering, width = 10, height = 5.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



