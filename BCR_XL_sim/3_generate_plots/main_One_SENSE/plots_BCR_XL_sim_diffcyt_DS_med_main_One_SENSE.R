##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: One-SENSE
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
library(cowplot)


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL_sim/main"

load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_med_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DS_med_main.RData"))
load(file.path(DIR_RDATA, "out_objects_BCR_XL_sim_diffcyt_DS_med_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/main_One_SENSE"




################
# Generate plots
################

# load plot data objects
d_se <- out_objects_diffcyt_DS_med_main$d_se
d_counts <- out_objects_diffcyt_DS_med_main$d_counts
d_medians_all <- out_objects_diffcyt_DS_med_main$d_medians_all


# run 1-D tSNE separately on identity markers and functional markers

d_tsne_id <- assay(d_medians_all)[, colData(d_medians_all)$is_identity_col]
d_tsne_id <- as.matrix(d_tsne_id)

d_tsne_func <- assay(d_medians_all)[, colData(d_medians_all)$is_func_col]
d_tsne_func <- as.matrix(d_tsne_func)

# remove any duplicate rows (required by Rtsne)
dups_id <- duplicated(d_tsne_id)
d_tsne_id <- d_tsne_id[!dups_id, ]

dups_func <- duplicated(d_tsne_func)
d_tsne_func <- d_tsne_func[!dups_func, ]

# run 1-D Rtsne
# (note: initial PCA step not required, since we do not have too many dimensions)
set.seed(123)
out_tsne_id <- Rtsne(d_tsne_id, dims = 1, pca = FALSE, verbose = TRUE)

set.seed(123)
out_tsne_func <- Rtsne(d_tsne_func, dims = 1, pca = FALSE, verbose = TRUE)

stopifnot(nrow(out_tsne_id) == nrow(out_tsne_func), 
          nrow(out_tsne_id) == nrow(rowData(d_counts)))

# data frame of One-SENSE coordinates
d_coords <- data.frame(coord_id = out_tsne_id$Y[, 1], 
                       coord_func = out_tsne_func$Y[, 1], 
                       cluster = 1:nrow(out_tsne_id$Y))

stopifnot(all(d_coords$cluster == rowData(d_counts)$cluster))


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

table(sig)

# set up data frame for plotting
d_plot <- data.frame(cluster = rowData(d_counts)$cluster, 
                     sig = as.numeric(sig), 
                     n_cells = rowData(d_counts)$n_cells)


# (ii) from cell-level results

# load true B-cell status of each cell
B_cells <- out_diffcyt_DS_med_main$B_cell

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
          nrow(d_plot) == nrow(d_coords))

# identify clusters containing significant proportion of spike-in cells
d_true$B_cells_0.5 <- as.numeric(d_true$prop_B_cells > 0.5)

# data frame for plotting
d_plot$prop_B_cells <- d_true$prop_B_cells
d_plot$B_cells_0.5 <- as.factor(d_true$B_cells_0.5)
d_plot$sig <- as.factor(d_plot$sig)

d_plot <- cbind(d_plot, d_coords[, c("coord_id", "coord_func")])


# (iii) create plots

p <- 
  ggplot(d_plot, aes(x = coord_func, y = coord_id, size = n_cells, color = sig)) + 
  # first layer
  geom_point(alpha = 0.5) + 
  scale_size_area(max_size = 3) + 
  scale_color_manual(values = c("gray70", "red"), labels = c("no", "yes")) + 
  # additional layer: outline clusters containing significant proportion true B cells
  geom_point(data = subset(d_plot, B_cells_0.5 == 1), aes(shape = B_cells_0.5), color = "black", stroke = 1.5) + 
  scale_shape_manual(values = 1, labels = ">50%") + 
  # additional layer: emphasize significant differential clusters
  geom_point(data = subset(d_plot, sig == 1), color = "red", alpha = 0.75) + 
  xlab("functional (1D tSNE)") + 
  ylab("identity (1D tSNE)") + 
  ggtitle("BCR-XL-sim, main results: diffcyt-DS-med: One-SENSE") + 
  theme_bw() + 
  theme(aspect.ratio = 1) + 
  guides(color = guide_legend("significant", override.aes = list(alpha = 1, size = 3), order = 1), 
         shape = guide_legend("true B cells", override.aes = list(size = 2, stroke = 1.25), order = 2), 
         size = guide_legend("no. cells", override.aes = list(color = "gray70", stroke = 0.25), order = 3))

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_med_main_One_SENSE.pdf")
ggsave(fn, width = 6, height = 5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



