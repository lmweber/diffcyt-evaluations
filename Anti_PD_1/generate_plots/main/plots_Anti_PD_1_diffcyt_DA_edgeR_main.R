##########################################################################################
# Generate plots
# 
# - data set: Anti-PD-1
# - plot type: phenotype and abundance of detected clusters
# - method: diffcyt-DA-edgeR
# 
# - main results
# 
# Lukas Weber, April 2018
##########################################################################################


library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)


# load saved results
DIR_RDATA <- "../../../../RData/Anti_PD_1/main"

load(file.path(DIR_RDATA, "out_clusters_Anti_PD_1_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA, "out_objects_Anti_PD_1_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA, "outputs_Anti_PD_1_diffcyt_DA_edgeR_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/Anti_PD_1/main"




#######################################################
# Heatmap: phenotype and abundance of detected clusters
#######################################################

d_se <- out_objects_diffcyt_DA_edgeR_main$d_se
d_counts <- out_objects_diffcyt_DA_edgeR_main$d_counts
d_medians_by_cluster_marker <- out_objects_diffcyt_DA_edgeR_main$d_medians_by_cluster_marker


# get detected clusters

res <- out_clusters_diffcyt_DA_edgeR_main

res <- res[!is.na(res$FDR), , drop = FALSE]

res_sorted <- res[order(res$FDR), , drop = FALSE]
head(res_sorted, 10)

res_detected <- res_sorted[res_sorted$FDR <= 0.1, , drop = FALSE]
res_detected

clus_detected <- res_detected$cluster_id


# medians for plotting

meds_detected <- assay(d_medians_by_cluster_marker)[clus_detected, , drop = FALSE]
if (length(clus_detected) > 0) {
  rownames(meds_detected) <- paste0("cluster ", rownames(meds_detected))
}

meds_all <- data.frame(t(colMeans(assay(d_se)[, colData(d_se)$marker_class != "none"])), 
                       check.names = FALSE)
rownames(meds_all) <- "all cells"

stopifnot(ncol(meds_detected) == ncol(meds_all), 
          all(colnames(meds_detected) == colnames(meds_all)))

meds <- rbind(meds_detected, meds_all)
meds <- meds[, order(colnames(meds))]


# --------------
# create heatmap
# --------------

# main panel showing cluster phenotype

colors_phenotype <- colorRamp2(quantile(assay(d_medians_by_cluster_marker), (c(0, 0.5, 1))), 
                               c("royalblue3", "white", "tomato2"))

split <- factor(c(rep("clusters", length(clus_detected)), "all cells"), levels = c("clusters", "all cells"))

ht_phenotype <- Heatmap(
  meds, col = colors_phenotype, name = "expression", 
  split = split, combined_name_fun = NULL, 
  column_title = "markers", column_title_side = "bottom", column_title_gp = gpar(fontsize = 12), 
  column_names_gp = gpar(fontsize = 11), 
  cluster_rows = TRUE, clustering_distance_rows = "euclidean", clustering_method_rows = "median", 
  cluster_columns = FALSE, 
  row_names_side = "left", row_names_gp = gpar(fontsize = 12), 
  heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12))
)


# second panel showning cluster abundance by sample

n_cells_detected <- assay(d_counts)[clus_detected, , drop = FALSE]
rownames(n_cells_detected) <- paste0("cluster ", rownames(n_cells_detected))

n_cells_plot <- rbind(n_cells_detected, "all cells" = NA)
ix_ord <- c(1:5, 11:14, 6:10, 15:20)
n_cells_plot <- n_cells_plot[, ix_ord]

colors_abundance <- colorRamp2(range(n_cells_detected, na.rm = TRUE), 
                               c("#132a13", "yellow"))

ht_abundance <- Heatmap(
  n_cells_plot, col = colors_abundance, name = "n_cells", 
  na_col = "white", 
  column_title = "samples", column_title_side = "bottom", column_title_gp = gpar(fontsize = 12), 
  column_names_gp = gpar(fontsize = 11), 
  cluster_rows = FALSE, cluster_columns = FALSE, 
  show_row_names = FALSE, 
  heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12))
)


# combine heatmaps and save

ht_title <- "Anti-PD-1, diffcyt-DA-edgeR: detected clusters"

fn <- file.path(DIR_PLOTS, "results_Anti_PD_1_diffcyt_DA_edgeR_main_heatmap.pdf")
pdf(fn, width = 12, height = 3)
draw(ht_phenotype + ht_abundance, newpage = FALSE, 
     column_title = ht_title, column_title_gp = gpar(fontface = "bold", fontsize = 12))
dev.off()




######################################################
# Boxplots: abundance of detected clusters (by sample)
######################################################

# combined abundance of detected clusters matching the population of interest

# select clusters that match the population of interest (from previous plot)
rownames_keep <- paste("cluster", c(358, 317, 380))

# calculate number of cells in selected clusters
n_cells_keep <- n_cells_detected[rownames_keep, ]
stopifnot(all(colnames(n_cells_keep) == colnames(assay(d_counts))), 
          all(colnames(n_cells_keep) == colData(d_counts)$sample_id))
# as percentages
perc <- colSums(n_cells_keep) / colSums(assay(d_counts)) * 100
perc
# average percentage in each group
perc_NR <- mean(perc[colData(d_counts)$group_id == "NR"])
perc_R <- mean(perc[colData(d_counts)$group_id == "R"])
perc_NR
perc_R

d_boxplots <- data.frame(
  group = colData(d_counts)$group_id, 
  percent = perc
)


# boxplots
ggplot(d_boxplots, aes(x = group, y = percent, color = group)) + 
  geom_boxplot(alpha = 0, width = 0.25) + 
  geom_point() + 
  labs(title = "Anti-PD-1, diffcyt-DA-edgeR", 
       subtitle = "Total abundance of detected clusters, by sample") + 
  theme_bw()

fn <- file.path(DIR_PLOTS, "results_Anti_PD_1_diffcyt_DA_edgeR_main_boxplots.pdf")
ggsave(file = fn, width = 4.75, height = 4)




#########
# Runtime
#########

d_runtimes <- as.data.frame(c(
  diffcyt_DA_edgeR = runtime_diffcyt_DA_edgeR_main
))

colnames(d_runtimes) <- "runtime"

d_runtimes$method <- factor(rownames(d_runtimes), levels = rownames(d_runtimes))

# color scheme
colors <- "darkblue"

y_range <- c(1, 1000)

# create plot
p_runtimes <- 
  ggplot(d_runtimes, aes(x = method, y = runtime, color = method, label = round(runtime, 1))) + 
  geom_point(shape = 4, size = 1.75, stroke = 1.5) + 
  geom_text(color = "black", vjust = -1.5, size = 3.4) + 
  scale_color_manual(values = colors) + 
  scale_y_log10(limits = y_range) + 
  ylab("runtime (sec, log10 scale)") + 
  ggtitle("Anti-PD-1, main results: runtime") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank())

# save plot
fn <- file.path(DIR_PLOTS, "results_Anti_PD_1_diffcyt_DA_edgeR_main_runtime.pdf")
ggsave(file = fn, width = 5, height = 3.25)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



