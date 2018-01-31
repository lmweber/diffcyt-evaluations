##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: t-SNE
# - method: diffcyt-DS-LMM
# 
# - main results
# 
# Lukas Weber, January 2018
##########################################################################################


library(SummarizedExperiment)
library(dplyr)
library(magrittr)
library(Rtsne)
library(ggplot2)
library(cowplot)


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL_sim/main"

load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_LMM_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DS_LMM_main.RData"))
load(file.path(DIR_RDATA, "out_objects_BCR_XL_sim_diffcyt_DS_LMM_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/main_tSNE"




################
# Generate plots
################

# load plot data objects
d_se <- out_objects_diffcyt_DS_LMM_main$d_se
d_counts <- out_objects_diffcyt_DS_LMM_main$d_counts
d_medians_all <- out_objects_diffcyt_DS_LMM_main$d_medians_all

# run t-SNE

# note: using cell type markers only
d_tsne <- assay(d_medians_all)[, colData(d_medians_all)$is_celltype_marker]
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


# (i) from cluster-level results

# load cluster-level results (note: marker pS6 only)
d_clus <- out_clusters_diffcyt_DS_LMM_main[out_clusters_diffcyt_DS_LMM_main$marker == "pS6", ]

stopifnot(nrow(d_clus) == nrow(rowData(d_counts)), 
          all(d_clus$cluster == rowData(d_counts)$cluster))

# significant differential clusters
cutoff_sig <- 0.1
sig <- d_clus$p_adj <= cutoff_sig
# set filtered clusters to FALSE
sig[is.na(sig)] <- FALSE

table(sig)

# set up data frame for plotting
d_plot <- data.frame(cluster = rowData(d_counts)$cluster, 
                     sig = as.numeric(sig), 
                     n_cells = rowData(d_counts)$n_cells)


# (ii) from cell-level results

# load true B-cell status of each cell
B_cells <- out_diffcyt_DS_LMM_main$B_cell

# calculate proportion true B-cells for each cluster
df_tmp <- as.data.frame(rowData(d_se))
stopifnot(nrow(df_tmp) == length(B_cells))

df_tmp$B_cells <- B_cells

d_true <- df_tmp %>% group_by(cluster) %>% summarize(prop_B_cells = mean(B_cells)) %>% as.data.frame

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

stopifnot(nrow(d_true) == nlevels(rowData(d_se)$cluster), 
          nrow(d_true) == nrow(d_plot), 
          all(d_true$cluster == d_plot$cluster), 
          nrow(d_plot) == nrow(tsne_coords))

# identify clusters containing significant proportion of spike-in cells
d_true$B_cells_0.99 <- as.numeric(d_true$prop_B_cells > 0.99)
d_true$B_cells_0.9 <- as.numeric(d_true$prop_B_cells > 0.9)
d_true$B_cells_0.5 <- as.numeric(d_true$prop_B_cells > 0.5)
d_true$B_cells_0.1 <- as.numeric(d_true$prop_B_cells > 0.1)

# data frame for plotting
d_plot$prop_B_cells <- d_true$prop_B_cells
d_plot$B_cells_0.99 <- as.factor(d_true$B_cells_0.99)
d_plot$B_cells_0.9 <- as.factor(d_true$B_cells_0.9)
d_plot$B_cells_0.5 <- as.factor(d_true$B_cells_0.5)
d_plot$B_cells_0.1 <- as.factor(d_true$B_cells_0.1)
d_plot$sig <- as.factor(d_plot$sig)
d_plot <- cbind(d_plot, tsne_coords)


# (iii) create plots

# one plot for each threshold defining 'true' clusters ('true' clusters are defined as
# containing at least x% true B cells)


# threshold 99%

p_0.99 <- 
  ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, size = n_cells, color = sig)) + 
  # first layer
  geom_point(alpha = 0.5) + 
  scale_size_area(max_size = 3) + 
  scale_color_manual(values = c("gray70", "red"), labels = c("no", "yes")) + 
  # additional layer: outline clusters containing significant proportion true B cells
  geom_point(data = subset(d_plot, B_cells_0.99 == 1), aes(shape = B_cells_0.99), color = "black", stroke = 1.5) + 
  scale_shape_manual(values = 1, labels = ">99%") + 
  # additional layer: emphasize significant differential clusters
  geom_point(data = subset(d_plot, sig == 1), color = "red", alpha = 0.75) + 
  xlab("t-SNE 1") + 
  ylab("t-SNE 2") + 
  ggtitle("BCR-XL-sim, main results: diffcyt-DS-LMM: t-SNE") + 
  theme_bw() + 
  theme(aspect.ratio = 1) + 
  guides(color = guide_legend("significant", override.aes = list(alpha = 1, size = 3), order = 1), 
         shape = guide_legend("true B cells", override.aes = list(size = 2, stroke = 1.25), order = 2), 
         size = guide_legend("no. cells", override.aes = list(color = "gray70", stroke = 0.25), order = 3))

# save plot
fn <- file.path(DIR_PLOTS, "panels", "results_BCR_XL_sim_diffcyt_DS_LMM_main_tSNE_0.99.pdf")
ggsave(fn, width = 6, height = 5)


# threshold 90%

p_0.9 <- 
  ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, size = n_cells, color = sig)) + 
  # first layer
  geom_point(alpha = 0.5) + 
  scale_size_area(max_size = 3) + 
  scale_color_manual(values = c("gray70", "red"), labels = c("no", "yes")) + 
  # additional layer: outline clusters containing significant proportion true B cells
  geom_point(data = subset(d_plot, B_cells_0.9 == 1), aes(shape = B_cells_0.9), color = "black", stroke = 1.5) + 
  scale_shape_manual(values = 1, labels = ">90%") + 
  # additional layer: emphasize significant differential clusters
  geom_point(data = subset(d_plot, sig == 1), color = "red", alpha = 0.75) + 
  xlab("t-SNE 1") + 
  ylab("t-SNE 2") + 
  ggtitle("BCR-XL-sim, main results: diffcyt-DS-LMM: t-SNE") + 
  theme_bw() + 
  theme(aspect.ratio = 1) + 
  guides(color = guide_legend("significant", override.aes = list(alpha = 1, size = 3), order = 1), 
         shape = guide_legend("true B cells", override.aes = list(size = 2, stroke = 1.25), order = 2), 
         size = guide_legend("no. cells", override.aes = list(color = "gray70", stroke = 0.25), order = 3))

# save plot
fn <- file.path(DIR_PLOTS, "panels", "results_BCR_XL_sim_diffcyt_DS_LMM_main_tSNE_0.9.pdf")
ggsave(fn, width = 6, height = 5)


# threshold 50%

p_0.5 <- 
  ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, size = n_cells, color = sig)) + 
  # first layer
  geom_point(alpha = 0.5) + 
  scale_size_area(max_size = 3) + 
  scale_color_manual(values = c("gray70", "red"), labels = c("no", "yes")) + 
  # additional layer: outline clusters containing significant proportion true B cells
  geom_point(data = subset(d_plot, B_cells_0.5 == 1), aes(shape = B_cells_0.5), color = "black", stroke = 1.5) + 
  scale_shape_manual(values = 1, labels = ">50%") + 
  # additional layer: emphasize significant differential clusters
  geom_point(data = subset(d_plot, sig == 1), color = "red", alpha = 0.75) + 
  xlab("t-SNE 1") + 
  ylab("t-SNE 2") + 
  ggtitle("BCR-XL-sim, main results: diffcyt-DS-LMM: t-SNE") + 
  theme_bw() + 
  theme(aspect.ratio = 1) + 
  guides(color = guide_legend("significant", override.aes = list(alpha = 1, size = 3), order = 1), 
         shape = guide_legend("true B cells", override.aes = list(size = 2, stroke = 1.25), order = 2), 
         size = guide_legend("no. cells", override.aes = list(color = "gray70", stroke = 0.25), order = 3))

# save plot
fn <- file.path(DIR_PLOTS, "panels", "results_BCR_XL_sim_diffcyt_DS_LMM_main_tSNE_0.5.pdf")
ggsave(fn, width = 6, height = 5)


# threshold 10%

p_0.1 <- 
  ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, size = n_cells, color = sig)) + 
  # first layer
  geom_point(alpha = 0.5) + 
  scale_size_area(max_size = 3) + 
  scale_color_manual(values = c("gray70", "red"), labels = c("no", "yes")) + 
  # additional layer: outline clusters containing significant proportion true B cells
  geom_point(data = subset(d_plot, B_cells_0.1 == 1), aes(shape = B_cells_0.1), color = "black", stroke = 1.5) + 
  scale_shape_manual(values = 1, labels = ">10%") + 
  # additional layer: emphasize significant differential clusters
  geom_point(data = subset(d_plot, sig == 1), color = "red", alpha = 0.75) + 
  xlab("t-SNE 1") + 
  ylab("t-SNE 2") + 
  ggtitle("BCR-XL-sim, main results: diffcyt-DS-LMM: t-SNE") + 
  theme_bw() + 
  theme(aspect.ratio = 1) + 
  guides(color = guide_legend("significant", override.aes = list(alpha = 1, size = 3), order = 1), 
         shape = guide_legend("true B cells", override.aes = list(size = 2, stroke = 1.25), order = 2), 
         size = guide_legend("no. cells", override.aes = list(color = "gray70", stroke = 0.25), order = 3))

# save plot
fn <- file.path(DIR_PLOTS, "panels", "results_BCR_XL_sim_diffcyt_DS_LMM_main_tSNE_0.1.pdf")
ggsave(fn, width = 6, height = 5)




##################
# Multi-panel plot
##################

plots_list <- list(p_0.99, p_0.9, p_0.5, p_0.1)

# modify plot elements
plots_list <- lapply(plots_list, function(p) {
  p + 
    theme(plot.title = element_blank())
})

plots_multi <- plot_grid(plotlist = plots_list, 
                         nrow = 2, ncol = 2, align = "hv", axis = "bl", scale = 0.98, 
                         labels = "AUTO", label_x = 0.01, label_y = 0.99)

# add combined title
title_single <- p_0.99$labels$title
plots_title <- ggdraw() + draw_label(title_single)
plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 20))

# save multi-panel plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_LMM_main_tSNE_all.pdf")
ggsave(fn, width = 10, height = 7.5)





###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



