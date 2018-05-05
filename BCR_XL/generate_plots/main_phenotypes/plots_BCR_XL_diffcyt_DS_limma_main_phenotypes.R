##########################################################################################
# Generate plots
# 
# - data set: BCR-XL
# - plot type: phenotypes and expression of cell state markers for top detected
# cluster-marker combinations
# - method: diffcyt-DS-limma
# 
# - main results
# 
# Lukas Weber, May 2018
##########################################################################################


library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL/main"

load(file.path(DIR_RDATA, "out_clusters_BCR_XL_diffcyt_DS_limma_main.RData"))
load(file.path(DIR_RDATA, "out_objects_BCR_XL_diffcyt_DS_limma_main.RData"))
load(file.path(DIR_RDATA, "outputs_BCR_XL_diffcyt_DS_limma_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL/main_phenotypes"




################
# Create heatmap
################

# load plot data objects
d_se <- out_objects_diffcyt_DS_limma_main$d_se
d_counts <- out_objects_diffcyt_DS_limma_main$d_counts
d_medians <- out_objects_diffcyt_DS_limma_main$d_medians
d_medians_by_cluster_marker <- out_objects_diffcyt_DS_limma_main$d_medians_by_cluster_marker


is_marker <- colData(d_medians_by_cluster_marker)$marker_class != "none"
is_celltype_marker <- colData(d_medians_by_cluster_marker)$marker_class == "type"


# -------------------------------------------------------------------
# First panel: phenotypes of top detected cluster-marker combinations
# -------------------------------------------------------------------

# note: show top 'n' clusters only
# note: no additional scaling (using asinh-transformed values directly)

d_heatmap <- assay(d_medians_by_cluster_marker)[, is_marker]

d_heatmap_celltype <- assay(d_medians_by_cluster_marker)[, is_celltype_marker]

# arrange alphabetically
d_heatmap_celltype <- d_heatmap_celltype[, order(colnames(d_heatmap_celltype))]

# load cluster-level results
d_clus <- out_clusters_diffcyt_DS_limma_main
# select top 'n' cluster-marker combinations
n <- 40
d_top <- d_clus[order(d_clus$p_adj)[seq_len(n)], , drop = FALSE]

d_heatmap_celltype <- d_heatmap_celltype[d_top$cluster_id, , drop = FALSE]

# color scale: 1%, 50%, 99% percentiles across all medians and all markers
colors <- colorRamp2(
  quantile(assay(d_medians_by_cluster_marker)[, is_marker], c(0.01, 0.5, 0.99)), 
  c("royalblue3", "white", "tomato2")
)

ht_phenotype <- Heatmap(
  d_heatmap_celltype, col = colors, name = "expression", 
  row_title = "clusters", row_title_gp = gpar(fontsize = 14), 
  column_title = "markers (cell type)", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
  show_column_names = TRUE, 
  column_names_gp = gpar(fontsize = 12), 
  heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
  cluster_columns = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 11), 
  cluster_rows = FALSE
)


# --------------------------------------------------------
# Second panel: expression of cell state markers by sample
# --------------------------------------------------------

# create data frame of expression values of top 'n' cluster-marker combinations by sample

top_clusters <- as.list(d_top$cluster_id)
top_markers <- as.character(d_top$marker)

assays_ordered <- assays(d_medians)[top_markers]

d_markers <- mapply(function(a, cl) {
  a[cl, , drop = FALSE]
}, assays_ordered, top_clusters, SIMPLIFY = FALSE)

d_markers <- do.call("rbind", d_markers)

stopifnot(all(rownames(d_markers) == rownames(d_heatmap_celltype)), 
          nrow(d_markers) == nrow(d_heatmap_celltype), 
          all(top_markers == d_top$marker))

rownames(d_markers) <- paste0(top_markers, " (", formatC(d_top$p_adj, format = "e", digits = 1), ")")

# re-order columns (by group)
cols_order <- c(seq(2, 16, by = 2), seq(1, 16, by = 2))
d_markers <- d_markers[, cols_order]

# color scale: full range
colors_markers <- colorRamp2(range(d_markers, na.rm = TRUE), c("navy", "yellow"))

# note: row ordering is automatically matched when multiple heatmaps are combined
ht_markers <- Heatmap(
  d_markers, col = colors_markers, name = "expression\n(by sample)",
  column_title = "samples", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
  column_names_gp = gpar(fontsize = 12), 
  heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
  cluster_columns = FALSE, 
  show_row_names = TRUE, row_names_side = "right", row_names_gp = gpar(fontsize = 11)
)


# --------------
# Combine panels
# --------------

ht_title <- "BCR-XL, diffcyt-DS-limma: top detected cluster-marker combinations"

fn <- file.path(DIR_PLOTS, "diffcyt_DS_limma", "results_BCR_XL_diffcyt_DS_limma_main_phenotypes.pdf")
pdf(fn, width = 9, height = 10)
draw(ht_phenotype + ht_markers, newpage = FALSE, 
     column_title = ht_title, column_title_gp = gpar(fontface = "bold", fontsize = 14))
dev.off()




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



