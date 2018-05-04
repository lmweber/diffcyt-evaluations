##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: t-SNE
# - method: diffcyt-DA-GLMM
# 
# - main results
# 
# Lukas Weber, May 2018
##########################################################################################


library(SummarizedExperiment)
library(dplyr)
library(magrittr)
library(Rtsne)
library(ggplot2)
library(cowplot)


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/main"

load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_GLMM_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_GLMM_main.RData"))
load(file.path(DIR_RDATA, "out_objects_AML_sim_diffcyt_DA_GLMM_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/main_tSNE"




################
# Generate plots
################

# loop over thresholds (th) and conditions (j)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# condition names
cond_names <- c("CN", "CBF")
cond_names_all <- c("healthy", cond_names)


# store plots in list
plots_tSNE <- vector("list", length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  # load plot data objects (same for both conditions j)
  d_se <- out_objects_diffcyt_DA_GLMM_main[[th]]$d_se
  d_counts <- out_objects_diffcyt_DA_GLMM_main[[th]]$d_counts
  d_medians_by_cluster_marker <- out_objects_diffcyt_DA_GLMM_main[[th]]$d_medians_by_cluster_marker
  
  # run t-SNE
  
  # note: using cell type markers only
  d_tsne <- assay(d_medians_by_cluster_marker)[, colData(d_medians_by_cluster_marker)$marker_class == "type"]
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
  tsne_coords <- as.data.frame(matrix(NA, nrow = nlevels(rowData(d_counts)$cluster_id), ncol = 2))
  rownames(tsne_coords) <- rowData(d_counts)$cluster_id
  colnames(tsne_coords) <- colnames(tsne_coords_tmp)
  tsne_coords[!dups, ] <- tsne_coords_tmp
  
  
  for (j in 1:length(cond_names)) {
    
    # index to store plots sequentially in list
    ix <- (j * length(thresholds)) - (length(thresholds) - th)
    
    
    # (i) from cluster-level results
    
    # load cluster-level results (for condition j)
    d_clus <- out_clusters_diffcyt_DA_GLMM_main[[th]][[j]]
    
    stopifnot(nrow(d_clus) == nrow(rowData(d_counts)), 
              all(d_clus$cluster_id == rowData(d_counts)$cluster_id))
    
    # significant differential clusters
    cutoff_sig <- 0.1
    sig <- d_clus$p_adj <= cutoff_sig
    # set filtered clusters to FALSE
    sig[is.na(sig)] <- FALSE
    
    table(sig)
    
    # set up data frame for plotting
    d_plot <- data.frame(cluster = rowData(d_counts)$cluster_id, 
                         sig = as.numeric(sig), 
                         n_cells = rowData(d_counts)$n_cells)
    
    
    # (ii) from cell-level results
    
    # load spike-in status at cell level (for condition j)
    spikein <- out_diffcyt_DA_GLMM_main[[th]][[j]]$spikein
    
    n_cells_cond <- rowData(d_se) %>% as.data.frame %>% group_by(group_id) %>% tally
    n_cells_cond <- unname(unlist(n_cells_cond[, "n"]))
    
    # calculate proportion true spike-in cells (from condition j) for each cluster
    
    # note: select cells from this condition and healthy
    cond_keep <- rowData(d_se)$group_id %in% c(cond_names[j], "healthy")
    df_j <- as.data.frame(rowData(d_se)[cond_keep, ])
    stopifnot(nrow(df_j) == length(spikein))
    
    df_j$spikein <- spikein
    
    d_true <- df_j %>% group_by(cluster_id) %>% summarize(prop_spikein = mean(spikein)) %>% as.data.frame
    
    # fill in any missing clusters (zero cells)
    if (nrow(d_true) < nlevels(rowData(d_se)$cluster_id)) {
      ix_missing <- which(!(levels(rowData(d_se)$cluster_id) %in% d_true$cluster_id))
      d_true_tmp <- data.frame(factor(ix_missing, levels = levels(rowData(d_se)$cluster_id)), 0)
      colnames(d_true_tmp) <- colnames(d_true)
      rownames(d_true_tmp) <- ix_missing
      d_true <- rbind(d_true, d_true_tmp)
      # re-order rows
      d_true <- d_true[order(d_true$cluster_id), ]
      rownames(d_true) <- d_true$cluster_id
    }
    
    stopifnot(nrow(d_true) == nlevels(rowData(d_se)$cluster_id), 
              nrow(d_true) == nrow(d_plot), 
              all(d_true$cluster_id == d_plot$cluster_id), 
              nrow(d_plot) == nrow(tsne_coords))
    
    # identify clusters containing significant proportion of spike-in cells
    d_true$spikein <- as.numeric(d_true$prop_spikein > 0.5)
    
    # data frame for plotting
    d_plot$prop_spikein <- d_true$prop_spikein
    d_plot$spikein <- as.factor(d_true$spikein)
    d_plot$sig <- as.factor(d_plot$sig)
    d_plot <- cbind(d_plot, tsne_coords)
    
    
    # (iii) create plots
    
    # note: 'true' clusters are defined as containing at least x% true spike-in cells
    
    p <- 
      ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, size = n_cells, color = sig)) + 
      # first layer
      geom_point(alpha = 0.5) + 
      scale_size_area(max_size = 3) + 
      scale_color_manual(values = c("gray70", "red"), labels = c("no", "yes")) + 
      # additional layer: outline clusters containing significant proportion spike-in cells
      geom_point(data = subset(d_plot, spikein == 1), aes(shape = spikein), color = "black", stroke = 1.5) + 
      scale_shape_manual(values = 1, labels = ">50%") + 
      # additional layer: emphasize significant differential clusters
      geom_point(data = subset(d_plot, sig == 1), color = "red", alpha = 1) + 
      xlab("t-SNE 1") + 
      ylab("t-SNE 2") + 
      ggtitle(paste0(cond_names[j], ", threshold ", gsub("pc$", "\\%", thresholds[th]))) + 
      theme_bw() + 
      theme(aspect.ratio = 1) + 
      guides(color = guide_legend("significant", override.aes = list(alpha = 1, size = 3), order = 1), 
             shape = guide_legend("true spike-in cells", override.aes = list(size = 2, stroke = 1.25), order = 2), 
             size = guide_legend("no. cells", override.aes = list(color = "gray70", stroke = 0.25), order = 3))
    
    plots_tSNE[[ix]] <- p
    
    # save individual panel plot
    p <- p + 
      ggtitle(paste0("AML-sim, main results: diffcyt-DA-GLMM: ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th]), ": t-SNE"))
    
    fn <- file.path(DIR_PLOTS, "panels", 
                    paste0("results_AML_sim_diffcyt_DA_GLMM_main_tSNE_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 6.25, height = 5.25)
    
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
grid_tSNE <- do.call(plot_grid, append(plots_tSNE, list(
  nrow = 2, ncol = 3, align = "hv", axis = "bl", scale = 0.975))
)

# add combined axis titles
xaxis_tSNE <- ggdraw() + draw_label("t-SNE dimension 1")
yaxis_tSNE <- ggdraw() + draw_label("t-SNE dimension 2", angle = 90)

grid_tSNE <- plot_grid(grid_tSNE, xaxis_tSNE, ncol = 1, rel_heights = c(15, 1))
grid_tSNE <- plot_grid(yaxis_tSNE, grid_tSNE, nrow = 1, rel_widths = c(1, 20))

# add combined title
title_tSNE <- ggdraw() + draw_label("AML-sim, main results: diffcyt-DA-GLMM: t-SNE", fontface = "bold")
grid_tSNE <- plot_grid(title_tSNE, grid_tSNE, ncol = 1, rel_heights = c(1, 18))

# add combined legend
legend_tSNE <- get_legend(plots_tSNE[[1]] + theme(
  legend.position = "right", legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12))
)
grid_tSNE <- plot_grid(grid_tSNE, legend_tSNE, nrow = 1, rel_widths = c(4, 1))

# save plots
fn_tSNE <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_DA_GLMM_main_tSNE.pdf")
ggsave(fn_tSNE, grid_tSNE, width = 8, height = 5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



