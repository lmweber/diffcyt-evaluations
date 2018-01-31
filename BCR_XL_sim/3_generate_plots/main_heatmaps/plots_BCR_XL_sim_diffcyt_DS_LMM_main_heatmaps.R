##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: heatmaps
# - method: diffcyt-DS-LMM
# 
# - main results
# 
# Lukas Weber, January 2018
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

load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_LMM_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DS_LMM_main.RData"))
load(file.path(DIR_RDATA, "out_objects_BCR_XL_sim_diffcyt_DS_LMM_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/main_heatmaps"




################
# Generate plots
################

# load plot data objects
d_se <- out_objects_diffcyt_DS_LMM_main$d_se
d_counts <- out_objects_diffcyt_DS_LMM_main$d_counts
d_medians_all <- out_objects_diffcyt_DS_LMM_main$d_medians_all


# -----------------------------------------------------------
# create heatmap: main panel showing marker expression values
# -----------------------------------------------------------

# note: using all markers
# note: no additional scaling (using asinh-transformed values directly)

d_heatmap <- assay(d_medians_all)[, colData(d_medians_all)$is_marker]

# arrange columns (cell type and state markers)
d_heatmap <- cbind(d_heatmap[, metadata(d_medians_all)$id_celltype_markers], 
                   d_heatmap[, metadata(d_medians_all)$id_state_markers])

# arrange each group (cell type and state) alphabetically
n_celltype <- sum(metadata(d_medians_all)$id_celltype_markers)
n_state <- sum(metadata(d_medians_all)$id_state_markers)
d_heatmap <- cbind(d_heatmap[, seq_len(n_celltype)][, order(colnames(d_heatmap)[seq_len(n_celltype)])], 
                   d_heatmap[, (seq_len(n_state)) + n_celltype][, order(colnames(d_heatmap)[(seq_len(n_state)) + n_celltype])])

# column annotation
col_annot <- data.frame(
  "marker type" = factor(c(rep("cell type", n_celltype), rep("state", n_state)), levels = c("cell type", "state")), 
  check.names = FALSE
)

ha_col <- columnAnnotation(df = col_annot, 
                           col = list("marker type" = c("cell type" = "gold", "state" = "darkgreen")), 
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
              clustering_distance_rows = "euclidean", clustering_method_rows = "median", 
              bottom_annotation = ha_col)


# -------------------------------------------------------------
# row annotation for significant and true differential clusters
# -------------------------------------------------------------

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

# set up data frame
d_sig <- data.frame(cluster = rowData(d_counts)$cluster, 
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
          nrow(d_true) == nrow(d_sig), 
          all(d_true$cluster == d_sig$cluster), 
          nrow(d_sig) == nrow(d_heatmap))

# identify clusters containing significant proportion of spike-in cells
d_true$B_cells <- as.numeric(d_true$prop_B_cells > 0.5)


# (iii) add row annotation and title

row_annot <- data.frame(
  "significant" = factor(d_sig$sig, levels = c(0, 1), labels = c("no", "yes")), 
  "true B cells" = factor(d_true$B_cells, levels = c(0, 1), labels = c("no", "yes")), 
  check.names = FALSE
)

ha_row <- rowAnnotation(df = row_annot, 
                        col = list("significant" = c("no" = "gray90", "yes" = "red"), 
                                   "true B cells" = c("no" = "gray90", "yes" = "black")), 
                        annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
                        width = unit(1, "cm"))

ht_title <- "BCR-XL-sim: diffcyt-DS-LMM"


# (iv) save individual plot

fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_LMM_main_heatmap.pdf")
pdf(fn, width = 8.25, height = 8)
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



