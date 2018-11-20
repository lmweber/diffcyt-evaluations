##########################################################################################
# Generate plots
# 
# - data set: BCR-XL
# - plot type: proportions of true populations and detected differential cell state
# markers for each cluster (using cutoff for adjusted p-values)
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
DIR_PLOTS <- "../../../../plots/BCR_XL/main_proportions"




################
# Create heatmap
################

# -------------------------------------------------------------
# First panel: proportions of true populations for each cluster
# -------------------------------------------------------------

# calculate proportions of true populations for each cluster

d_se <- out_objects_diffcyt_DS_limma_main$d_se

population_IDs <- rowData(d_se)$population_id
length(population_IDs)
table(population_IDs)

cluster_IDs <- rowData(d_se)$cluster_id

d_heatmap_props <- as.data.frame.matrix(t(table(population_IDs, cluster_IDs)))
stopifnot(nrow(d_heatmap_props) == nlevels(cluster_IDs))

# convert to proportions per cluster
n_cells_clus <- rowSums(d_heatmap_props)
stopifnot(length(n_cells_clus) == nlevels(cluster_IDs))

d_heatmap_props <- d_heatmap_props / n_cells_clus

# color scale: black and white
colors <- colorRamp2(range(d_heatmap_props, na.rm = TRUE), c("gray95", "black"))

ht_props <- Heatmap(
  d_heatmap_props, col = colors, name = "proportion", 
  row_title = "clusters", row_title_gp = gpar(fontsize = 14), 
  row_names_side = "left", row_names_gp = gpar(fontsize = 8), 
  column_title = "populations", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
  column_names_gp = gpar(fontsize = 12), 
  heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
  cluster_columns = FALSE, 
  clustering_distance_rows = "euclidean", clustering_method_rows = "average"
)


# -----------------------------------------------------------------------
# Second panel: detected differential cell state markers for each cluster
# -----------------------------------------------------------------------

# detected differential cell state markers for each cluster

res <- out_clusters_diffcyt_DS_limma_main
table(res$p_adj <= 0.1)
dim(res)

n_clus <- nlevels(cluster_IDs)
n_state_markers <- sum(metadata(out_objects_diffcyt_DS_limma_main$d_medians)$id_state_markers)

# note: showing detected markers at a given significance cutoff
cutoff <- 0.1

d_detected_markers <- matrix(FALSE, nrow = n_clus, ncol = n_state_markers)
marker_names <- sort(as.character(unique(res$marker)))
stopifnot(length(marker_names) == sum(n_state_markers))
colnames(d_detected_markers) <- marker_names
rownames(d_detected_markers) <- levels(cluster_IDs)

for (i in seq_along(marker_names)) {
  res_i <- res[res$marker == marker_names[i], ]
  res_i_detected <- res_i[res_i$p_adj <= cutoff, ]
  clus_i_detected <- res_i_detected$cluster_id
  d_detected_markers[clus_i_detected, marker_names[i]] <- TRUE
}


# color scale: red and white
colors <- colorRamp2(range(d_heatmap_props, na.rm = TRUE), c("gray95", "red"))

ht_markers <- Heatmap(
  d_detected_markers, col = colors, name = "detected", 
  row_title = "clusters", row_title_gp = gpar(fontsize = 14), 
  show_row_names = FALSE, 
  column_title = "markers (cell state)", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
  column_names_gp = gpar(fontsize = 12), 
  heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12), 
                              color_bar = "discrete", at = c(0, 1), labels = c("no", "yes")), 
  cluster_columns = FALSE, cluster_rows = FALSE
)


# --------------
# Combine panels
# --------------

ht_title <- "BCR-XL, diffcyt-DS-limma: cluster proportions and detected markers"

fn <- file.path(DIR_PLOTS, "diffcyt_DS_limma", "results_BCR_XL_diffcyt_DS_limma_main_proportions_cutoff.pdf")
pdf(fn, width = 10, height = 11)
draw(ht_props + ht_markers, newpage = FALSE, 
     column_title = ht_title, column_title_gp = gpar(fontface = "bold", fontsize = 14))
dev.off()




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



