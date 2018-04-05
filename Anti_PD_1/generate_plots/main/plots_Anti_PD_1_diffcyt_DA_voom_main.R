##########################################################################################
# Generate plots
# 
# - data set: Anti-PD-1
# - plot type: phenotype and abundance of detected clusters
# - method: diffcyt-DA-voom
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

load(file.path(DIR_RDATA, "out_clusters_Anti_PD_1_diffcyt_DA_voom_main.RData"))
load(file.path(DIR_RDATA, "out_objects_Anti_PD_1_diffcyt_DA_voom_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/Anti_PD_1/main"




#########################################
# Heatmap: phenotype of detected clusters
#########################################

d_se <- out_objects_diffcyt_DA_voom_main$d_se
d_counts <- out_objects_diffcyt_DA_voom_main$d_counts
d_medians_by_cluster_marker <- out_objects_diffcyt_DA_voom_main$d_medians_by_cluster_marker


# get detected clusters

res <- out_clusters_diffcyt_DA_voom_main

res <- res[!is.na(res$adj.P.Val), , drop = FALSE]

res_sorted <- res[order(res$adj.P.Val), , drop = FALSE]
head(res_sorted, 10)

res_detected <- res_sorted[res_sorted$adj.P.Val <= 0.1, , drop = FALSE]
res_detected

clus_detected <- res_detected$cluster


# medians for plotting

meds_detected <- assay(d_medians_by_cluster_marker)[clus_detected, , drop = FALSE]
rownames(meds_detected) <- paste0("cluster ", rownames(meds_detected))

meds_all <- colMeans(assay(d_se)[, colData(d_se)$is_marker])

stopifnot(ncol(meds_detected) == ncol(meds_all))

meds <- rbind(meds_detected, "all cells" = meds_all)
meds <- meds[, order(colnames(meds))]


# create heatmap
colors <- colorRamp2(quantile(assay(d_medians_by_cluster_marker), (c(0, 0.5, 1))), 
                     c("royalblue3", "white", "tomato2"))

split <- factor(c(rep("clusters", length(clus_detected)), "all cells"), levels = c("clusters", "all cells"))

ht_main <- Heatmap(
  meds, col = colors, name = "expression", 
  split = split, combined_name_fun = NULL, 
  column_title = "markers", column_title_side = "bottom", column_title_gp = gpar(fontsize = 12), 
  column_names_gp = gpar(fontsize = 12), 
  cluster_rows = FALSE, cluster_columns = FALSE, 
  row_names_side = "left", row_names_gp = gpar(fontsize = 12), 
  heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12))
)

ht_main_title <- "Anti-PD-1, diffcyt-DA-voom: phenotype of detected clusters"


# save plot
fn <- file.path(DIR_PLOTS, "results_Anti_PD_1_diffcyt_DA_voom_main_heatmap_phenotype.pdf")
pdf(fn, width = 10, height = 7)
draw(ht_main, newpage = FALSE, 
     column_title = ht_main_title, column_title_gp = gpar(fontface = "bold", fontsize = 12))
dev.off()




#########################################
# Heatmap: abundance of detected clusters
#########################################

n_cells_detected <- assay(d_counts)[clus_detected, , drop = FALSE]
rownames(n_cells_detected) <- paste0("cluster ", rownames(n_cells_detected))


# create heatmap
colors_counts <- colorRamp2(range(n_cells_detected), 
                            c("#132a13", "yellow"))

ht_abundance <- Heatmap(
  n_cells_detected, col = colors_counts, name = "n_cells", 
  column_title = "samples", column_title_side = "bottom", column_title_gp = gpar(fontsize = 12), 
  column_names_gp = gpar(fontsize = 12), 
  cluster_rows = FALSE, cluster_columns = FALSE, 
  row_names_side = "left", row_names_gp = gpar(fontsize = 12), 
  heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12))
)

ht_abundance_title <- "Anti-PD-1, diffcyt-DA-voom: abundance of detected clusters"


# save plot
fn <- file.path(DIR_PLOTS, "results_Anti_PD_1_diffcyt_DA_voom_main_heatmap_abundance.pdf")
pdf(fn, width = 9, height = 7)
draw(ht_abundance, newpage = FALSE, 
     column_title = ht_abundance_title, column_title_gp = gpar(fontface = "bold", fontsize = 12))
dev.off()




##########################################
# Boxplots: abundance of detected clusters
##########################################

# abundance of all detected clusters combined


# number of cells in detected clusters
colSums(n_cells_detected)
colSums(assay(d_counts))
# as percentage
perc <- colSums(n_cells_detected) / colSums(assay(d_counts)) * 100
perc
# average percentage in each group
perc_NR <- mean(perc[colData(d_counts)$group == "NR"])
perc_R <- mean(perc[colData(d_counts)$group == "R"])
perc_NR
perc_R

d_boxplots <- data.frame(
  group = colData(d_counts)$group, 
  percent = perc
)


# boxplots
ggplot(d_boxplots, aes(x = group, y = percent, color = group)) + 
  geom_boxplot(alpha = 0, width = 0.25) + 
  geom_point() + 
  labs(title = "Anti-PD-1, diffcyt-DA-voom: combined abundance") + 
  theme_bw()

fn <- file.path(DIR_PLOTS, "results_Anti_PD_1_diffcyt_DA_voom_main_boxplots.pdf")
ggsave(file = fn, width = 5, height = 4)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



