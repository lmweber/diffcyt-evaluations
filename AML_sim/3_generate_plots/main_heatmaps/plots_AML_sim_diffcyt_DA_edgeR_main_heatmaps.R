##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: heatmaps
# - method: diffcyt-DA-edgeR
# 
# - main results
# 
# Lukas Weber, March 2018
##########################################################################################


library(SummarizedExperiment)
library(dplyr)
library(magrittr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/main"

load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA, "out_objects_AML_sim_diffcyt_DA_edgeR_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/main_heatmaps"




################
# Generate plots
################

# loop over thresholds (th) and conditions (j)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# condition names
cond_names <- c("CN", "CBF")
cond_names_all <- c("healthy", cond_names)


# store plots in list
plots_heatmaps <- vector("list", length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  # load plot data objects (same for both conditions j)
  d_se <- out_objects_diffcyt_DA_edgeR_main[[th]]$d_se
  d_counts <- out_objects_diffcyt_DA_edgeR_main[[th]]$d_counts
  d_medians_all <- out_objects_diffcyt_DA_edgeR_main[[th]]$d_medians_all
  
  
  for (j in 1:length(cond_names)) {
    
    # index to store plots sequentially in list
    ix <- (j * length(thresholds)) - (length(thresholds) - th)
    
    
    # -----------------------------------------------------------
    # create heatmap: main panel showing marker expression values
    # -----------------------------------------------------------
    
    # note: using clustering markers only
    # note: show top 'n' clusters only (otherwise heatmaps are too small on multi-panel plot)
    # note: no additional scaling (using asinh-transformed values directly)
    
    d_heatmap <- assay(d_medians_all)[, colData(d_medians_all)$is_celltype_marker]
    
    # load cluster-level results (for condition j)
    d_clus <- out_clusters_diffcyt_DA_edgeR_main[[th]][[j]]
    
    stopifnot(nrow(d_clus) == nrow(rowData(d_counts)), 
              all(d_clus$cluster == rowData(d_counts)$cluster))
    
    # select top 'n' clusters
    n <- 20
    top_n <- order(d_clus$FDR)[1:n]
    
    d_heatmap <- d_heatmap[top_n, ]
    
    # use 1% and 99% percentiles for color scale
    colors <- colorRamp2(quantile(d_heatmap, c(0.01, 0.5, 0.99)), 
                         c("royalblue3", "white", "tomato2"))
    
    ht_main <- Heatmap(
      d_heatmap, col = colors, name = "expression", 
      row_title = "clusters", row_title_gp = gpar(fontsize = 14), 
      column_title = "markers (cell type)", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
      column_names_gp = gpar(fontsize = 12), 
      heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
      cluster_columns = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 11), 
      clustering_distance_rows = "euclidean", clustering_method_rows = "median"
    )
    
    
    # --------------------------------------------
    # second heatmap: cluster abundances by sample
    # --------------------------------------------
    
    cnd_which <- c(which(colData(d_counts)$group_IDs == "healthy"), 
                   which(colData(d_counts)$group_IDs == cond_names[j]))
    
    d_abundance <- assay(d_counts)[top_n, cnd_which, drop = FALSE]
    
    stopifnot(all(rownames(d_heatmap) == rownames(d_abundance)), 
              nrow(d_heatmap) == nrow(d_abundance))
    
    # use full range for color scale
    colors_counts <- colorRamp2(range(d_abundance), 
                                c("#132a13", "yellow"))
    
    # note: row ordering is automatically matched when multiple heatmaps are combined
    ht_abundance <- Heatmap(
      d_abundance, col = colors_counts, name = "n_cells", 
      column_title = "samples", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
      column_names_gp = gpar(fontsize = 12), 
      heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
      cluster_columns = FALSE, show_row_names = FALSE
    )
    
    
    # -------------------------------------------------------------
    # row annotation for significant and true differential clusters
    # -------------------------------------------------------------
    
    # (i) from cluster-level results
    
    # significant differential clusters
    cutoff_sig <- 0.1
    sig <- d_clus$FDR <= cutoff_sig
    # set filtered clusters to FALSE
    sig[is.na(sig)] <- FALSE
    
    table(sig)
    
    # set up data frame for plotting
    d_plot <- data.frame(cluster = rowData(d_counts)$cluster, 
                         sig = as.numeric(sig), 
                         n_cells = rowData(d_counts)$n_cells)
    
    # select top 'n' clusters
    d_plot <- d_plot[top_n, ]
    
    
    # (ii) from cell-level results
    
    # load spike-in status at cell level (for condition j)
    spikein <- out_diffcyt_DA_edgeR_main[[th]][[j]]$spikein
    
    n_cells_cond <- rowData(d_se) %>% as.data.frame %>% group_by(group_IDs) %>% tally
    n_cells_cond <- unname(unlist(n_cells_cond[, "n"]))
    
    # calculate proportion true spike-in cells (from condition j) for each cluster
    
    # note: select cells from this condition and healthy
    cond_keep <- rowData(d_se)$group_IDs %in% c(cond_names[j], "healthy")
    df_j <- as.data.frame(rowData(d_se)[cond_keep, ])
    stopifnot(nrow(df_j) == length(spikein))
    
    df_j$spikein <- spikein
    
    d_true <- df_j %>% group_by(cluster) %>% summarize(prop_spikein = mean(spikein)) %>% as.data.frame
    
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
    
    stopifnot(nrow(d_true) == nlevels(rowData(d_se)$cluster))
    
    # select top 'n' clusters
    d_true <- d_true[top_n, ]
    
    stopifnot(nrow(d_true) == nrow(d_plot), 
              all(d_true$cluster == d_plot$cluster), 
              nrow(d_plot) == nrow(d_heatmap))
    
    # identify clusters containing significant proportion of spike-in cells
    d_true$spikein <- as.numeric(d_true$prop_spikein > 0.5)
    
    # data frame for plotting
    d_plot$prop_spikein <- d_true$prop_spikein
    d_plot$spikein <- as.factor(d_true$spikein)
    d_plot$sig <- as.factor(d_plot$sig)
    d_plot <- cbind(d_plot, d_heatmap)
    
    
    # (iii) add row annotation and title
    
    row_annot <- data.frame(
      "significant" = factor(d_plot$sig, levels = c(0, 1), labels = c("no", "yes")), 
      "spike-in" = factor(d_true$spikein, levels = c(0, 1), labels = c("no", "yes")), 
      check.names = FALSE
    )
    
    ha_row <- rowAnnotation(
      df = row_annot, 
      col = list("significant" = c("no" = "gray90", "yes" = "darkorange1"), 
                 "spike-in" = c("no" = "gray90", "yes" = "black")), 
      annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
      width = unit(1.2, "cm")
    )
    
    ht_title <- paste0("AML-sim, ", cond_names[j], ", threshold ", gsub("pc$", "\\%", thresholds[th]), ": diffcyt-DA-edgeR")
    
    
    # (iv) save individual plot
    
    fn <- file.path(DIR_PLOTS, paste0("panels/results_AML_sim_diffcyt_DA_edgeR_main_heatmap_AML_sim_", thresholds[th], "_", cond_names[j], ".pdf"))
    pdf(fn, width = 9, height = 6.5)
    plots_heatmaps[[ix]] <- draw(ht_main + ht_abundance + ha_row, newpage = FALSE, 
                                 column_title = ht_title, column_title_gp = gpar(fontface = "bold", fontsize = 14))
    dev.off()
    
  }
}




########################
# Save multi-panel plots
########################

fn <- file.path(DIR_PLOTS, paste0("results_AML_sim_diffcyt_DA_edgeR_main_heatmaps.pdf"))
pdf(fn, width = 25, height = 10)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 3)))

for (th in 1:length(thresholds)) {
  for (j in 1:length(cond_names)) {
    ix <- (j * length(thresholds)) - (length(thresholds) - th)
    
    pushViewport(viewport(layout.pos.row = j, layout.pos.col = th))
    draw(plots_heatmaps[[ix]], newpage = FALSE)
    upViewport()
  }
}

upViewport()

dev.off()




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



