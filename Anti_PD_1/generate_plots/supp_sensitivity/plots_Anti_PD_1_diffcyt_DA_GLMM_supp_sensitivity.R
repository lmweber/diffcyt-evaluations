##########################################################################################
# Generate plots
# 
# - data set: Anti-PD-1
# - plot type: phenotype and abundance of detected clusters
# - method: diffcyt-DA-GLMM
# 
# - supplementary results: sensitivity to random seeds
# 
# Lukas Weber, May 2018
##########################################################################################


library(SummarizedExperiment)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)


# load saved results
DIR_RDATA <- "../../../../RData/Anti_PD_1/supp_sensitivity"

load(file.path(DIR_RDATA, "out_clusters_Anti_PD_1_diffcyt_DA_GLMM_supp_sensitivity.RData"))
load(file.path(DIR_RDATA, "out_objects_Anti_PD_1_diffcyt_DA_GLMM_supp_sensitivity.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/Anti_PD_1/supp_sensitivity"


# varying random seeds
seeds <- as.numeric(names(out_objects_diffcyt_DA_GLMM_supp_sensitivity))




#######################################################
# Heatmap: phenotype and abundance of detected clusters
#######################################################

# note: one heatmap for each random seed

plot_widths <- c(14, 13, 14, 14, 14)
plot_heights <- c(12, 12, 12, 12, 12)


for (s in 1:length(seeds)) {
  
  d_se <- out_objects_diffcyt_DA_GLMM_supp_sensitivity[[s]]$d_se
  d_counts <- out_objects_diffcyt_DA_GLMM_supp_sensitivity[[s]]$d_counts
  d_medians_by_cluster_marker <- out_objects_diffcyt_DA_GLMM_supp_sensitivity[[s]]$d_medians_by_cluster_marker
  
  
  # get detected clusters
  
  res <- out_clusters_diffcyt_DA_GLMM_supp_sensitivity[[s]]
  
  res <- res[!is.na(res$p_adj), , drop = FALSE]
  
  res_sorted <- res[order(res$p_adj), , drop = FALSE]
  head(res_sorted, 10)
  
  res_detected <- res_sorted[res_sorted$p_adj <= 0.1, , drop = FALSE]
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
    cluster_rows = TRUE, clustering_distance_rows = "euclidean", clustering_method_rows = "average", 
    cluster_columns = FALSE, 
    row_names_side = "left", row_names_gp = gpar(fontsize = 12), 
    heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12))
  )
  
  
  # second panel showning cluster abundance by sample
  
  n_cells_detected <- assay(d_counts)[clus_detected, , drop = FALSE]
  if (length(n_cells_detected) > 0) {
    rownames(n_cells_detected) <- paste0("cluster ", rownames(n_cells_detected))
  }
  
  n_cells_plot <- rbind(n_cells_detected, "all cells" = NA)
  if (!is.data.frame(n_cells_plot)) {
    n_cells_plot <- as.data.frame(n_cells_plot)
  }
  ix_ord <- c(1:5, 11:14, 6:10, 15:20)
  n_cells_plot <- n_cells_plot[, ix_ord]
  
  n_cells_plot <- n_cells_plot / apply(n_cells_plot, 1, max)
  
  if (length(n_cells_detected) > 0) {
    colors_abundance <- colorRamp2(range(n_cells_plot, na.rm = TRUE), c("#132a13", "yellow"))
  } else {
    colors_abundance <- colorRamp2(c(0, 1), c("#132a13", "yellow"))
  }
  
  ht_abundance <- Heatmap(
    n_cells_plot, col = colors_abundance, name = "prop_cells", 
    na_col = "white", 
    column_title = "samples", column_title_side = "bottom", column_title_gp = gpar(fontsize = 12), 
    column_names_gp = gpar(fontsize = 11), 
    cluster_rows = FALSE, cluster_columns = FALSE, 
    show_row_names = FALSE, 
    heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12))
  )
  
  
  # combine heatmaps and save
  
  ht_title <- paste0("Anti-PD-1, diffcyt-DA-GLMM: detected clusters: random seed ", s)
  
  fn <- file.path(DIR_PLOTS, "diffcyt_DA_GLMM", paste0("results_Anti_PD_1_diffcyt_DA_GLMM_supp_sensitivity_heatmap_seed_", s, ".pdf"))
  pdf(fn, width = plot_widths[s], height = plot_heights[s])
  draw(ht_phenotype + ht_abundance, newpage = FALSE, 
       column_title = ht_title, column_title_gp = gpar(fontface = "bold", fontsize = 12))
  dev.off()
  
}




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



