##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: tSNE
# - method: diffcyt-DA-limma
# 
# - main results
# 
# Lukas Weber, October 2017
##########################################################################################


library(SummarizedExperiment)
library(dplyr)
library(magrittr)
library(Rtsne)
library(ggplot2)
library(cowplot)


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/main"

load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_limma_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_limma_main.RData"))
load(file.path(DIR_RDATA, "out_objects_AML_sim_diffcyt_DA_limma_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/main_tSNE"




################
# Generate plots
################

# loop over thresholds (th) and conditions (j)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")
cond_names_all <- c("healthy", cond_names)


# store plots in list
plots_tSNE <- vector("list", length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  # load plot data objects (same for both conditions j)
  d_se <- out_objects_diffcyt_DA_limma_main[[th]]$d_se
  d_counts <- out_objects_diffcyt_DA_limma_main[[th]]$d_counts
  d_medians_all <- out_objects_diffcyt_DA_limma_main[[th]]$d_medians_all
  
  # run tSNE
  
  # note: using clustering markers only
  d_tsne <- assay(d_medians_all)[, colData(d_medians_all)$is_clustering_col]
  d_tsne <- as.matrix(d_tsne)
  
  # remove any duplicate rows (required by Rtsne)
  dups <- duplicated(d_tsne)
  d_tsne <- d_tsne[!dups, ]
  
  # run Rtsne
  # (note: initial PCA step not required, since we do not have too many dimensions)
  set.seed(123)
  out_tsne <- Rtsne(d_tsne, pca = FALSE, verbose = TRUE)
  
  tsne_coords_tmp <- as.data.frame(out_tsne$Y)
  colnames(tsne_coords_tmp) <- c("tSNE_1", "tSNE_2")
  
  # fill in any missing clusters due to duplicate rows
  tsne_coords <- as.data.frame(matrix(NA, nrow = nlevels(rowData(d_counts)$cluster), ncol = 2))
  rownames(tsne_coords) <- rowData(d_counts)$cluster
  colnames(tsne_coords) <- colnames(tsne_coords_tmp)
  tsne_coords[!dups, ] <- tsne_coords_tmp
  
  
  for (j in 1:length(cond_names)) {
    
    # index to store plots sequentially in list
    ix <- (th * length(cond_names)) - (length(cond_names) - j)
    
    
    # (i) from cluster-level results
    
    # load cluster-level results (for condition j)
    d_clus <- out_clusters_diffcyt_DA_limma_main[[th]][[j]]
    stopifnot(nrow(d_clus) == nrow(rowData(d_counts)))
    stopifnot(all(d_clus$cluster == rowData(d_counts)$cluster))
    
    # significant differential clusters
    cutoff_sig <- 0.1
    sig <- d_clus$adj.P.Val <= cutoff_sig
    # set filtered clusters to FALSE
    sig[is.na(sig)] <- FALSE
    
    # set up data frame for plotting
    d_plot <- data.frame(cluster = rowData(d_counts)$cluster, 
                         sig = as.numeric(sig), 
                         n_cells = rowData(d_counts)$n_cells)
    
    
    # (ii) from cell-level results
    
    # load spike-in status at cell level (for condition j)
    spikein <- out_diffcyt_DA_limma_main[[th]][[j]]$spikein
    
    n_cells_cond <- rowData(d_se) %>% as.data.frame %>% group_by(group) %>% tally
    n_cells_cond <- unname(unlist(n_cells_cond[, "n"]))
    
    # identify spike-in cells (for condition j) within full data set (all conditions)
    spikein_list <- vector("list", length(cond_names_all))
    for (s in 1:length(spikein_list)) {
      if (cond_names_all[s] == cond_names[j]) {
        spikein_list[[s]] <- spikein
      } else {
        spikein_list[[s]] <- rep(0, n_cells_cond[s])
      }
    }
    spikein_all <- unlist(spikein_list)
    
    # calculate proportion true spike-in cells (from condition j) for each cluster
    df_j <- as.data.frame(rowData(d_se))
    stopifnot(nrow(df_j) == length(spikein_all))
    
    df_j$spikein <- spikein_all
    
    d_true <- df_j %>% group_by(cluster) %>% summarize(prop_spikein = mean(spikein)) %>% as.data.frame
    
    # fill in any missing clusters (zero cells)
    if (nrow(d_true) < nlevels(rowData(d_se)$cluster)) {
      ix_missing <- which(!(levels(rowData(d_se)$cluster) %in% d_true$cluster))
      d_true_tmp <- data.frame(factor(ix_missing, levels = levels(rowData(d_se)$cluster)), 0)
      colnames(d_true_tmp) <- colnames(d_true)
      rownames(d_true_tmp) <- ix_missing
      d_true <- rbind(d_true, d_true_tmp)
      # re-order rows
      d_true <- d_true[order(d_true$cluster), ]
      rownames(d_true) <- d_true$cluster
    }
    
    stopifnot(nrow(d_true) == nlevels(rowData(d_se)$cluster))
    stopifnot(nrow(d_true) == nrow(d_plot))
    stopifnot(all(d_true$cluster == d_plot$cluster))
    stopifnot(nrow(d_plot) == nrow(tsne_coords))
    
    # identify clusters containing significant proportion of spike-in cells
    cutoff_prop <- 0.1
    d_true$spikein <- as.numeric(d_true$prop_spikein > cutoff_prop)
    
    # data frame for plotting
    d_plot$prop_spikein <- d_true$prop_spikein
    d_plot$spikein <- as.factor(d_true$spikein)
    d_plot$sig <- as.factor(d_plot$sig)
    d_plot <- cbind(d_plot, tsne_coords)
    
    
    # (iii) create plot
    
    # with short title for multi-panel plot
    p <- ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, size = n_cells, color = sig)) + 
      # first layer
      geom_point(alpha = 0.5) + 
      scale_size_area(max_size = 3) + 
      scale_color_manual(values = c("gray70", "red"), labels = c("FALSE", "TRUE")) + 
      # additional layer: outline clusters containing significant proportion spike-in cells
      geom_point(data = subset(d_plot, spikein == 1), shape = 1, color = "black", stroke = 0.75) + 
      # additional layer: emphasize significant differential clusters
      geom_point(data = subset(d_plot, sig == 1), color = "red", alpha = 0.5) + 
      ggtitle(paste0(cond_names[j], ", threshold ", gsub("pc$", "\\%", thresholds[th]))) + 
      theme_bw() + 
      theme(aspect.ratio = 1) + 
      guides(color = guide_legend("significant"), 
             size = guide_legend(override.aes = list(color = "gray70", stroke = 0.25)))
    
    plots_tSNE[[ix]] <- p
    
    # save individual panel plot
    p <- p + 
      ggtitle(paste0("AML-sim, main results: diffcyt-DA-limma: ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th]), ": tSNE"))
    
    fn <- file.path(DIR_PLOTS, "panels", 
                    paste0("results_AML_sim_diffcyt_DA_limma_main_tSNE_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 7.5, height = 6)
    
  }
}




########################
# Save multi-panel plots
########################

# modify plot elements
plots_tSNE <- lapply(plots_tSNE, function(p) {
  p + theme(legend.position = "none", 
            axis.title.x = element_blank(), axis.title.y = element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

# format into grid
grid_tSNE <- do.call(plot_grid, append(plots_tSNE, list(nrow = 4, ncol = 2, 
                                                        align = "hv", axis = "bl", 
                                                        scale = 0.975)))

# add combined axis titles
xaxis_tSNE <- ggdraw() + draw_label("tSNE dimension 1", size = 14)
yaxis_tSNE <- ggdraw() + draw_label("tSNE dimension 2", size = 14, angle = 90)

grid_tSNE <- plot_grid(grid_tSNE, xaxis_tSNE, ncol = 1, rel_heights = c(50, 1))
grid_tSNE <- plot_grid(yaxis_tSNE, grid_tSNE, nrow = 1, rel_widths = c(1, 30))

# add combined title
title_tSNE <- ggdraw() + draw_label("AML-sim, main results: diffcyt-DA-limma: tSNE", fontface = "bold")
grid_tSNE <- plot_grid(title_tSNE, grid_tSNE, ncol = 1, rel_heights = c(1, 32))

# add combined legend
legend_tSNE <- get_legend(plots_tSNE[[1]] + theme(legend.position = "right", 
                                                  legend.title = element_text(size = 12, face = "bold"), 
                                                  legend.text = element_text(size = 12)))
grid_tSNE <- plot_grid(grid_tSNE, legend_tSNE, nrow = 1, rel_widths = c(5, 1))

# save plots
fn_tSNE <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_DA_limma_main_tSNE.pdf")
ggsave(fn_tSNE, grid_tSNE, width = 10, height = 13)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



