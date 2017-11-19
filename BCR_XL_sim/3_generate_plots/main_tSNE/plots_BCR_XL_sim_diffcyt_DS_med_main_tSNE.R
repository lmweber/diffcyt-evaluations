##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: tSNE
# - method: diffcyt-DS-med
# 
# - main results
# 
# Lukas Weber, November 2017
##########################################################################################


library(SummarizedExperiment)
library(dplyr)
library(magrittr)
library(Rtsne)
library(ggplot2)


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL_sim/main"

load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_med_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DS_med_main.RData"))
load(file.path(DIR_RDATA, "out_objects_BCR_XL_sim_diffcyt_DS_med_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/main_tSNE"




################
# Generate plots
################

# load plot data objects
d_se <- out_objects_diffcyt_DS_med_main$d_se
d_counts <- out_objects_diffcyt_DS_med_main$d_counts
d_medians_all <- out_objects_diffcyt_DS_med_main$d_medians_all

# run tSNE

# note: using identity markers only
d_tsne <- assay(d_medians_all)[, colData(d_medians_all)$is_identity_col]
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
d_clus <- out_clusters_diffcyt_DS_med_main[out_clusters_diffcyt_DS_med_main$marker == "pS6(Yb172)Dd", ]

stopifnot(nrow(d_clus) == nrow(rowData(d_counts)), 
          all(d_clus$cluster == rowData(d_counts)$cluster))

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

# load spike-in status at cell level
spikein <- out_diffcyt_DS_med_main$spikein

# include cells from both 'base' and 'spike' conditions
spikein_all <- c(rep(0, sum(rowData(d_se)$group == "base")), spikein)

# calculate proportion true spike-in cells for each cluster
df_tmp <- as.data.frame(rowData(d_se))
stopifnot(nrow(df_tmp) == length(spikein_all))

df_tmp$spikein <- spikein_all

d_true <- df_tmp %>% group_by(cluster) %>% summarize(prop_spikein = mean(spikein)) %>% as.data.frame

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
cutoff_prop <- 0.1
d_true$spikein <- as.numeric(d_true$prop_spikein > cutoff_prop)

# data frame for plotting
d_plot$prop_spikein <- d_true$prop_spikein
d_plot$spikein <- as.factor(d_true$spikein)
d_plot$sig <- as.factor(d_plot$sig)
d_plot <- cbind(d_plot, tsne_coords)


# (iii) create plot

p <- 
  ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, size = n_cells, color = sig)) + 
  # first layer
  geom_point(alpha = 0.5) + 
  scale_size_area(max_size = 3) + 
  scale_color_manual(values = c("gray70", "red"), labels = c("no", "yes")) + 
  # additional layer: outline clusters containing significant proportion spike-in cells
  geom_point(data = subset(d_plot, spikein == 1), aes(shape = spikein), color = "black", stroke = 2) + 
  scale_shape_manual(values = 1, labels = ">10%") + 
  # additional layer: emphasize significant differential clusters
  geom_point(data = subset(d_plot, sig == 1), color = "red", alpha = 0.75) + 
  ggtitle("BCR-XL-sim, main results: diffcyt-DS-med: tSNE") + 
  theme_bw() + 
  theme(aspect.ratio = 1) + 
  guides(color = guide_legend("significant", override.aes = list(alpha = 1, size = 3), order = 1), 
         shape = guide_legend("spike-in", override.aes = list(size = 1.5), order = 2), 
         size = guide_legend("no. cells", override.aes = list(color = "gray70", stroke = 0.25), order = 3))

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_med_main_tSNE.pdf")
ggsave(fn, width = 6.25, height = 5.25)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



