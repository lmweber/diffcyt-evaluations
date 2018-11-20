##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: heatmaps
# - method: diffcyt-DA-edgeR
# 
# - supplementary results: treating all markers as 'cell type' markers; using DA methods
# 
# Lukas Weber, November 2018
##########################################################################################


library(SummarizedExperiment)
library(dplyr)
library(magrittr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL_sim/supp_cell_type_markers"

load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DA_edgeR_supp_cell_type_markers.RData"))
load(file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DA_edgeR_supp_cell_type_markers.RData"))
load(file.path(DIR_RDATA, "out_objects_BCR_XL_sim_diffcyt_DA_edgeR_supp_cell_type_markers.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/supp_cell_type_markers"




################
# Generate plots
################

# load plot data objects
d_se <- out_objects_diffcyt_DA_edgeR_supp_cell_type_markers$d_se
d_counts <- out_objects_diffcyt_DA_edgeR_supp_cell_type_markers$d_counts
d_medians <- out_objects_diffcyt_DA_edgeR_supp_cell_type_markers$d_medians
d_medians_by_cluster_marker <- out_objects_diffcyt_DA_edgeR_supp_cell_type_markers$d_medians_by_cluster_marker


is_marker <- colData(d_medians_by_cluster_marker)$marker_class != "none"
is_celltype_marker <- colData(d_medians_by_cluster_marker)$marker_class == "type"
is_state_marker <- colData(d_medians_by_cluster_marker)$marker_class == "state"


# -------------------------------------------------------------
# heatmap: main panel showing expression of 'cell type' markers
# -------------------------------------------------------------

# note: using clustering markers only
# note: show top 'n' clusters only (otherwise heatmaps are too small on multi-panel plot)
# note: no additional scaling (using asinh-transformed values directly)

d_heatmap <- assay(d_medians_by_cluster_marker)[, is_marker]

d_heatmap_celltype <- assay(d_medians_by_cluster_marker)[, is_celltype_marker]

# arrange alphabetically
d_heatmap_celltype <- d_heatmap_celltype[, order(colnames(d_heatmap_celltype))]

# load cluster-level results
d_clus <- out_clusters_diffcyt_DA_edgeR_supp_cell_type_markers
stopifnot(nrow(d_clus) == nrow(d_heatmap_celltype))
# select top 'n' clusters
n <- 20
top_n <- order(d_clus$p_adj)[1:n]
d_heatmap_celltype <- d_heatmap_celltype[top_n, ]

# color scale: 1%, 50%, 99% percentiles across all medians and all markers
colors <- colorRamp2(
  quantile(assay(d_medians_by_cluster_marker)[, is_marker], c(0.01, 0.5, 0.99)), 
  c("royalblue3", "white", "tomato2")
)

ht_main <- Heatmap(
  d_heatmap_celltype, col = colors, name = "expression", 
  row_title = "clusters", row_title_gp = gpar(fontsize = 14), 
  column_title = "markers (cell type)", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
  show_column_names = TRUE, 
  column_names_gp = gpar(fontsize = 12), 
  heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
  cluster_columns = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 11), 
  clustering_distance_rows = "euclidean", clustering_method_rows = "average"
)


# ----------------------------------------------------------
# heatmap: second panel showing cluster abundances by sample
# ----------------------------------------------------------

cnd_which <- c(which(colData(d_counts)$group_id == "base"), 
               which(colData(d_counts)$group_id == "spike"))

d_abundance <- assay(d_counts)[top_n, cnd_which, drop = FALSE]

# scale to show relative counts per cluster
d_abundance <- d_abundance / apply(d_abundance, 1, max)

stopifnot(all(rownames(d_heatmap_celltype) == rownames(d_abundance)), 
          nrow(d_heatmap_celltype) == nrow(d_abundance))

# use full range for color scale
colors_counts <- colorRamp2(range(d_abundance), 
                            c("#132a13", "yellow"))

# note: row ordering is automatically matched when multiple heatmaps are combined
ht_abundance <- Heatmap(
  d_abundance, col = colors_counts, name = "prop_cells", 
  column_title = "samples", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
  column_names_gp = gpar(fontsize = 12), 
  heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
  cluster_columns = FALSE, show_row_names = FALSE
)


# --------------------------------------------------------
# heatmap: third panel showing expression of pS6 by sample
# --------------------------------------------------------

cnd_which <- c(which(colData(d_counts)$group_id == "base"), 
               which(colData(d_counts)$group_id == "spike"))

d_pS6 <- assays(d_medians)[["pS6"]][top_n, cnd_which, drop = FALSE]

stopifnot(all(rownames(d_heatmap_celltype) == rownames(d_pS6)), 
          nrow(d_heatmap_celltype) == nrow(d_pS6))

# use full range for color scale
colors_pS6 <- colorRamp2(range(d_pS6, na.rm = TRUE), 
                         c("navy", "yellow"))

# note: row ordering is automatically matched when multiple heatmaps are combined
ht_pS6 <- Heatmap(
  d_pS6, col = colors_pS6, name = "expression: pS6",
  column_title = "samples", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
  show_column_names = TRUE, 
  column_names_gp = gpar(fontsize = 12), 
  heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
  cluster_columns = FALSE, show_row_names = FALSE
)


# -------------------------------------------------------------
# row annotation for significant and true differential clusters
# -------------------------------------------------------------

# (i) from cluster-level results

# load cluster-level results (note: marker pS6 only)
d_clus

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
B_cells <- out_diffcyt_DA_edgeR_supp_cell_type_markers$B_cell

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
  "detected" = factor(d_sig$sig, levels = c(0, 1), labels = c("no", "yes")), 
  "true B cells" = factor(d_true$B_cells, levels = c(0, 1), labels = c("no", "yes")), 
  check.names = FALSE
)

ha_row <- rowAnnotation(
  df = row_annot, 
  col = list("detected" = c("no" = "gray90", "yes" = "red"), 
             "true B cells" = c("no" = "gray90", "yes" = "black")), 
  annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
  width = unit(1.2, "cm")
)

ht_title <- "BCR-XL-sim: diffcyt-DA-edgeR"


# (iv) save individual plot

fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DA_edgeR_supp_cell_type_markers_heatmap.pdf")
pdf(fn, width = 15, height = 6.5)
draw(ht_main + ht_abundance + ht_pS6 + ha_row, newpage = FALSE, 
     column_title = ht_title, column_title_gp = gpar(fontface = "bold", fontsize = 14))
dev.off()




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



