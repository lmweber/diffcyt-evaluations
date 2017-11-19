##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: heatmaps
# - method: diffcyt-DS-med
# 
# - main results
# 
# Lukas Weber, November 2017
##########################################################################################


library(SummarizedExperiment)
library(dplyr)
library(magrittr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL_sim/main"

load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_med_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DS_med_main.RData"))
load(file.path(DIR_RDATA, "out_objects_BCR_XL_sim_diffcyt_DS_med_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/main_heatmaps"




################
# Generate plots
################

# load plot data objects
d_se <- out_objects_diffcyt_DS_med_main$d_se
d_counts <- out_objects_diffcyt_DS_med_main$d_counts
d_medians_all <- out_objects_diffcyt_DS_med_main$d_medians_all


# -----------------------------------------------------------
# create heatmap: main panel showing marker expression values
# -----------------------------------------------------------

# note: using all markers
# note: no additional scaling (using asinh-transformed values directly)

d_heatmap <- assay(d_medians_all)[, colData(d_medians_all)$is_marker_col]

colnames(d_heatmap) <- gsub("\\(.*$", "", colnames(d_heatmap))

# arrange columns (identity and functional markers)
d_heatmap <- cbind(d_heatmap[, metadata(d_medians_all)$id_identity_markers], 
                   d_heatmap[, metadata(d_medians_all)$id_func_markers])

# arrange each group (identity and functional) alphabetically
n_identity <- sum(metadata(d_medians_all)$id_identity_markers)
n_func <- sum(metadata(d_medians_all)$id_func_markers)
d_heatmap <- cbind(d_heatmap[, seq_len(n_identity)][, order(colnames(d_heatmap)[seq_len(n_identity)])], 
                   d_heatmap[, (seq_len(n_func)) + n_identity][, order(colnames(d_heatmap)[(seq_len(n_func)) + n_identity])])

# column annotation
col_annot <- data.frame(
  marker_type = factor(c(rep("identity", n_identity), rep("functional", n_func)), levels = c("identity", "functional"))
)

ha_col <- columnAnnotation(df = col_annot, 
                           col = list(marker_type = c("identity" = "gold", "functional" = "darkgreen")), 
                           colname = anno_text(colnames(d_heatmap), rot = 90, just = "right", offset = unit(1, "npc") - unit(2, "mm")), 
                           annotation_height = unit.c(unit(4, "mm"), max_text_width(colnames(d_heatmap)) + unit(2, "mm")), 
                           annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)))

# use 1% and 99% percentiles for color scale
colors <- colorRamp2(quantile(d_heatmap, c(0.01, 0.5, 0.99)), 
                   c("royalblue", "white", "tomato"))

ht <- Heatmap(d_heatmap, col = colors, name = "expression", 
              row_title = "clusters", row_title_gp = gpar(fontsize = 14), 
              column_title = "markers", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
              show_column_names = FALSE, 
              column_names_gp = gpar(fontsize = 12), 
              heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
              cluster_columns = FALSE, show_row_names = FALSE, 
              bottom_annotation = ha_col)


# -------------------------------------------------------------
# row annotation for significant and true differential clusters
# -------------------------------------------------------------

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

# set up data frame
d_sig <- data.frame(cluster = rowData(d_counts)$cluster, 
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
          nrow(d_true) == nrow(d_sig), 
          all(d_true$cluster == d_sig$cluster), 
          nrow(d_sig) == nrow(d_heatmap))

# identify clusters containing significant proportion of spike-in cells
cutoff_prop <- 0.1
d_true$spikein <- as.numeric(d_true$prop_spikein > cutoff_prop)


# (iii) add row annotation and title

row_annot <- data.frame(
  significant = factor(d_sig$sig, levels = c(0, 1), labels = c("no", "yes")), 
  spikein = factor(d_true$spikein, levels = c(0, 1), labels = c("no", "yes"))
)

ha_row <- rowAnnotation(df = row_annot, 
                        col = list(significant = c("no" = "gray90", "yes" = "red"), 
                                   spikein = c("no" = "gray90", "yes" = "black")), 
                        annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
                        width = unit(1, "cm"))

ht_title <- "BCR-XL-sim: diffcyt-DS-med"


# (iv) save individual plot

fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_med_main_heatmap.pdf")
pdf(fn, width = 8, height = 8)
draw(ht + ha_row, newpage = FALSE, 
     column_title = ht_title, column_title_gp = gpar(fontface = "bold", fontsize = 14))
dev.off()




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



