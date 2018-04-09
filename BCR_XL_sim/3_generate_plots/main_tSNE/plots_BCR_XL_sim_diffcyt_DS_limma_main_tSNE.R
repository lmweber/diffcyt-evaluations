##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: t-SNE
# - method: diffcyt-DS-limma
# 
# - main results
# 
# Lukas Weber, April 2018
##########################################################################################


library(SummarizedExperiment)
library(dplyr)
library(magrittr)
library(Rtsne)
library(ggplot2)
library(cowplot)


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL_sim/main"

load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_limma_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DS_limma_main.RData"))
load(file.path(DIR_RDATA, "out_objects_BCR_XL_sim_diffcyt_DS_limma_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/main_tSNE"




################
# Generate plots
################

# load plot data objects
d_se <- out_objects_diffcyt_DS_limma_main$d_se
d_counts <- out_objects_diffcyt_DS_limma_main$d_counts
d_medians_by_cluster_marker <- out_objects_diffcyt_DS_limma_main$d_medians_by_cluster_marker

# run t-SNE

# note: using cell type markers only
d_tsne <- assay(d_medians_by_cluster_marker)[, colData(d_medians_by_cluster_marker)$marker_class == "cell_type"]
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


# (i) from cluster-level results

# load cluster-level results (note: marker pS6 only)
d_clus <- out_clusters_diffcyt_DS_limma_main[out_clusters_diffcyt_DS_limma_main$marker == "pS6", ]

stopifnot(nrow(d_clus) == nrow(rowData(d_counts)), 
          all(d_clus$cluster_id == rowData(d_counts)$cluster_id))

# significant differential clusters
cutoff_sig <- 0.1
sig <- d_clus$adj.P.Val <= cutoff_sig
# set filtered clusters to FALSE
sig[is.na(sig)] <- FALSE

table(sig)

# set up data frame for plotting
d_plot <- data.frame(cluster = rowData(d_counts)$cluster_id, 
                     sig = as.numeric(sig), 
                     n_cells = rowData(d_counts)$n_cells)


# (ii) from cell-level results

# load true B-cell status of each cell
B_cells <- out_diffcyt_DS_limma_main$B_cell

# calculate proportion true B-cells for each cluster
df_tmp <- as.data.frame(rowData(d_se))
stopifnot(nrow(df_tmp) == length(B_cells))

df_tmp$B_cells <- B_cells

d_true <- df_tmp %>% group_by(cluster) %>% summarize(prop_B_cells = mean(B_cells)) %>% as.data.frame

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
d_true$B_cells <- as.numeric(d_true$prop_B_cells > 0.5)

# data frame for plotting
d_plot$prop_B_cells <- d_true$prop_B_cells
d_plot$B_cells <- as.factor(d_true$B_cells)
d_plot$sig <- as.factor(d_plot$sig)
d_plot <- cbind(d_plot, tsne_coords)


# (iii) create plots

# note: 'true' clusters are defined as containing at least x% true B cells

p <- 
  ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, size = n_cells, color = sig)) + 
  # first layer
  geom_point(alpha = 0.5) + 
  scale_size_area(max_size = 3) + 
  scale_color_manual(values = c("gray70", "darkorange1"), labels = c("no", "yes")) + 
  # additional layer: outline clusters containing significant proportion true B cells
  geom_point(data = subset(d_plot, B_cells == 1), aes(shape = B_cells), color = "black", stroke = 1.5) + 
  scale_shape_manual(values = 1, labels = ">50%") + 
  # additional layer: emphasize significant differential clusters
  geom_point(data = subset(d_plot, sig == 1), color = "darkorange1", alpha = 1) + 
  xlab("t-SNE 1") + 
  ylab("t-SNE 2") + 
  ggtitle("BCR-XL-sim, main results: diffcyt-DS-limma: t-SNE") + 
  theme_bw() + 
  theme(aspect.ratio = 1) + 
  guides(color = guide_legend("significant", override.aes = list(alpha = 1, size = 3), order = 1), 
         shape = guide_legend("true B cells", override.aes = list(size = 2, stroke = 1.25), order = 2), 
         size = guide_legend("no. cells", override.aes = list(color = "gray70", stroke = 0.25), order = 3))

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_limma_main_tSNE.pdf")
ggsave(fn, width = 5, height = 4)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



