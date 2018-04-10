##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: heatmaps
# - method: diffcyt-DS-LMM
# 
# - main results
# 
# Lukas Weber, April 2018
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
d_medians <- out_objects_diffcyt_DS_LMM_main$d_medians
d_medians_by_cluster_marker <- out_objects_diffcyt_DS_LMM_main$d_medians_by_cluster_marker


is_marker <- colData(d_medians_by_cluster_marker)$marker_class != "none"
is_celltype_marker <- colData(d_medians_by_cluster_marker)$marker_class == "cell_type"
is_state_marker <- colData(d_medians_by_cluster_marker)$marker_class == "cell_state"


# -------------------------------------------------------
# heatmap: main panel - expression of 'cell type' markers
# -------------------------------------------------------

# note: show top 'n' clusters only (otherwise heatmaps are too small on multi-panel plot)
# note: no additional scaling (using asinh-transformed values directly)

d_heatmap <- assay(d_medians_by_cluster_marker)[, is_marker]

d_heatmap_celltype <- assay(d_medians_by_cluster_marker)[, is_celltype_marker]

# arrange alphabetically
d_heatmap_celltype <- d_heatmap_celltype[, order(colnames(d_heatmap_celltype))]

# load cluster-level results
d_clus <- out_clusters_diffcyt_DS_LMM_main[out_clusters_diffcyt_DS_LMM_main$marker == "pS6", ]
stopifnot(nrow(d_clus) == nrow(d_heatmap_celltype))
# select top 'n' clusters
n <- 20
top_n <- order(d_clus$p_adj)[1:n]
d_heatmap_celltype <- d_heatmap_celltype[top_n, ]

# column annotation for cell type and cell state markers
n_celltype <- sum(metadata(d_medians_by_cluster_marker)$id_type_markers)
n_state <- sum(metadata(d_medians_by_cluster_marker)$id_state_markers)

col_annot_celltype <- data.frame(
  "marker type" = factor(c(rep("cell type", n_celltype), rep("state", 0)), levels = c("cell type", "state")), 
  check.names = FALSE
)

ha_col_celltype <- columnAnnotation(
  df = col_annot_celltype, show_legend = FALSE, 
  col = list("marker type" = c("cell type" = "gold", "state" = "darkgreen")), 
  colname = anno_text(colnames(d_heatmap_celltype), rot = 90, just = "right", offset = unit(1, "npc") - unit(2, "mm")), 
  annotation_height = unit.c(unit(5, "mm"), max_text_width(colnames(d_heatmap_celltype)) + unit(2, "mm")), 
  annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12))
)

# color scale: 1%, 50%, 99% percentiles across all medians and all markers
colors <- colorRamp2(
  quantile(assay(d_medians_by_cluster_marker)[, is_marker], c(0.01, 0.5, 0.99)), 
  c("royalblue3", "white", "tomato2")
)

ht_main <- Heatmap(
  d_heatmap_celltype, col = colors, name = "expression", 
  row_title = "clusters", row_title_gp = gpar(fontsize = 14), 
  column_title = "markers (cell type)", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
  show_column_names = FALSE, 
  column_names_gp = gpar(fontsize = 12), 
  heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
  cluster_columns = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 11), 
  clustering_distance_rows = "euclidean", clustering_method_rows = "median", 
  bottom_annotation = ha_col_celltype
)


# -----------------------------------------------------
# heatmap: second panel - expression of 'state' markers
# -----------------------------------------------------

d_heatmap_state <- assay(d_medians_by_cluster_marker)[, is_state_marker]

# arrange alphabetically
d_heatmap_state <- d_heatmap_state[, order(colnames(d_heatmap_state))]

# select top 'n' clusters
d_heatmap_state <- d_heatmap_state[top_n, ]

# column annotation for cell type and cell state markers
col_annot_state <- data.frame(
  "marker type" = factor(c(rep("cell type", 0), rep("state", n_state)), levels = c("cell type", "state")), 
  check.names = FALSE
)

ha_col_state <- columnAnnotation(
  df = col_annot_state, show_legend = FALSE, 
  col = list("marker type" = c("cell type" = "gold", "state" = "darkgreen")), 
  colname = anno_text(colnames(d_heatmap_state), rot = 90, just = "right", offset = unit(1, "npc") - unit(2, "mm")), 
  annotation_height = unit.c(unit(5, "mm"), max_text_width(colnames(d_heatmap_state)) + unit(2, "mm")), 
  annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12))
)

ht_state <- Heatmap(
  d_heatmap_state, col = colors, 
  show_heatmap_legend = FALSE, 
  column_title = "markers (state)", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
  show_column_names = FALSE, 
  column_names_gp = gpar(fontsize = 12), 
  cluster_columns = FALSE, show_row_names = FALSE, 
  bottom_annotation = ha_col_state
)


# --------------------------------------------------
# heatmap: third panel - expression of pS6 by sample
# --------------------------------------------------

cnd_which <- c(which(colData(d_counts)$group_id == "base"), 
               which(colData(d_counts)$group_id == "spike"))

d_pS6 <- assays(d_medians)[["pS6"]][top_n, cnd_which, drop = FALSE]

stopifnot(all(rownames(d_heatmap_celltype) == rownames(d_pS6)), 
          nrow(d_heatmap_celltype) == nrow(d_pS6))

# column annotation (empty; need for formatting purposes)
col_annot_pS6 <- data.frame(
  "state marker" = factor(rep("pS6", ncol(d_pS6)), levels = "pS6"), 
  check.names = FALSE
)

ha_col_pS6 <- columnAnnotation(
  df = col_annot_pS6, show_legend = FALSE, 
  col = list("state marker" = c("pS6" = "gray")), 
  colname = anno_text(colnames(d_pS6), rot = 90, just = "right", offset = unit(1, "npc") - unit(2, "mm")), 
  # use zero height
  annotation_height = unit.c(unit(0, "mm"), max_text_width(colnames(d_pS6)) + unit(2, "mm")), 
  annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12))
)

# use full range for color scale
colors_pS6 <- colorRamp2(range(d_pS6), 
                         c("navy", "yellow"))

# note: row ordering is automatically matched when multiple heatmaps are combined
ht_pS6 <- Heatmap(
  d_pS6, col = colors_pS6, name = "expression: pS6",
  column_title = "samples", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
  show_column_names = FALSE, 
  column_names_gp = gpar(fontsize = 12), 
  heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
  cluster_columns = FALSE, show_row_names = FALSE, 
  bottom_annotation = ha_col_pS6
)


# -------------------------------------------------------------
# row annotation for significant and true differential clusters
# -------------------------------------------------------------

# (i) from cluster-level results

# load cluster-level results (note: marker pS6 only)
d_clus <- out_clusters_diffcyt_DS_LMM_main[out_clusters_diffcyt_DS_LMM_main$marker == "pS6", ]

stopifnot(nrow(d_clus) == nrow(rowData(d_counts)), 
          all(d_clus$cluster_id == rowData(d_counts)$cluster_id))

# significant differential clusters
cutoff_sig <- 0.1
sig <- d_clus$p_adj <= cutoff_sig
# set filtered clusters to FALSE
sig[is.na(sig)] <- FALSE

table(sig)

# set up data frame
d_sig <- data.frame(cluster = rowData(d_counts)$cluster_id, 
                    sig = as.numeric(sig), 
                    n_cells = rowData(d_counts)$n_cells)


# select top n
d_sig <- d_sig[top_n, ]

stopifnot(all(rownames(d_sig) == rownames(d_heatmap_celltype)), 
          nrow(d_sig) == nrow(d_heatmap_celltype))


# (ii) from cell-level results

# load true B-cell status of each cell
B_cells <- out_diffcyt_DS_LMM_main$B_cell

# calculate proportion true B-cells for each cluster
df_tmp <- as.data.frame(rowData(d_se))
stopifnot(nrow(df_tmp) == length(B_cells))

df_tmp$B_cells <- B_cells

d_true <- df_tmp %>% group_by(cluster_id) %>% summarize(prop_B_cells = mean(B_cells)) %>% as.data.frame

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
          all(d_true$cluster_id == rowData(d_counts)$cluster_id))

# identify clusters containing significant proportion of spike-in cells
d_true$B_cells <- as.numeric(d_true$prop_B_cells > 0.5)


# select top n
d_true <- d_true[top_n, ]

stopifnot(all(rownames(d_true) == rownames(d_heatmap_celltype)), 
          nrow(d_true) == nrow(d_heatmap_celltype))


# (iii) add row annotation and title

row_annot <- data.frame(
  "significant" = factor(d_sig$sig, levels = c(0, 1), labels = c("no", "yes")), 
  "true B cells" = factor(d_true$B_cells, levels = c(0, 1), labels = c("no", "yes")), 
  check.names = FALSE
)

ha_row <- rowAnnotation(
  df = row_annot, 
  col = list("significant" = c("no" = "gray90", "yes" = "darkorange1"), 
             "true B cells" = c("no" = "gray90", "yes" = "black")), 
  annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
  width = unit(1.2, "cm")
)

ht_title <- "BCR-XL-sim: diffcyt-DS-LMM"


# (iv) save individual plot

fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_LMM_main_heatmap.pdf")
pdf(fn, width = 12.5, height = 6.5)
draw(ht_main + ht_state + ht_pS6 + ha_row, newpage = FALSE, 
     column_title = ht_title, column_title_gp = gpar(fontface = "bold", fontsize = 14))
dev.off()




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



