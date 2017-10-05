##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: heatmaps (showing significant and true spike-in clusters only)
# - method: diffcyt-DA-GLMM
# 
# - main results
# 
# Lukas Weber, October 2017
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

load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_GLMM_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_GLMM_main.RData"))
load(file.path(DIR_RDATA, "out_objects_AML_sim_diffcyt_DA_GLMM_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/main_heatmaps"




################
# Generate plots
################

# loop over thresholds (th) and conditions (j)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")
cond_names_all <- c("healthy", cond_names)


# store plots in list
plots_heatmaps <- plots_heights <- vector("list", length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  # load plot data objects (same for both conditions j)
  d_se <- out_objects_diffcyt_DA_GLMM_main[[th]]$d_se
  d_counts <- out_objects_diffcyt_DA_GLMM_main[[th]]$d_counts
  d_medians_all <- out_objects_diffcyt_DA_GLMM_main[[th]]$d_medians_all
  
  
  for (j in 1:length(cond_names)) {
    
    # index to store plots sequentially in list
    ix <- (th * length(cond_names)) - (length(cond_names) - j)
    
    
    # ------------
    # heatmap data
    # ------------
    
    # note: using clustering markers only
    # note: no additional scaling (using asinh-transformed values directly)
    d_heatmap <- assay(d_medians_all)[, colData(d_medians_all)$is_clustering_col]
    
    colnames(d_heatmap) <- gsub("\\(.*$", "", colnames(d_heatmap))
    
    
    # -------------------------------------------------------------
    # row annotation for significant and true differential clusters
    # -------------------------------------------------------------
    
    # (i) from cluster-level results
    
    # load cluster-level results (for condition j)
    d_clus <- out_clusters_diffcyt_DA_GLMM_main[[th]][[j]]
    stopifnot(nrow(d_clus) == nrow(rowData(d_counts)))
    stopifnot(all(d_clus$cluster == rowData(d_counts)$cluster))
    
    # significant differential clusters
    cutoff_sig <- 0.1
    sig <- d_clus$p_adj <= cutoff_sig
    # set filtered clusters to FALSE
    sig[is.na(sig)] <- FALSE
    
    # set up data frame
    d_sig <- data.frame(cluster = rowData(d_counts)$cluster, 
                        sig = as.numeric(sig), 
                        n_cells = rowData(d_counts)$n_cells)
    
    
    # (ii) from cell-level results
    
    # load spike-in status at cell level (for condition j)
    spikein <- out_diffcyt_DA_GLMM_main[[th]][[j]]$spikein
    
    n_cells_cond <- rowData(d_se) %>% as.data.frame %>% group_by(group) %>% tally
    n_cells_cond <- unname(unlist(n_cells_cond[, "n"]))
    
    # identify spike-in cells (for condition j) within full data set (all conditions)
    spikein_list <- vector("list", length(cond_names_all))
    for (s in 1:length(spikein_list)) {
      if (cond_names_all[s] == cond_names[j]) {
        spikein_list[[s]] <- spikein
      } else {
        spikein_list[[s]] <- rep(0, n_cells_cond[s])
      }
    }
    spikein_all <- unlist(spikein_list)
    
    # calculate proportion true spike-in cells (from condition j) for each cluster
    df_j <- as.data.frame(rowData(d_se))
    stopifnot(nrow(df_j) == length(spikein_all))
    
    df_j$spikein <- spikein_all
    
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
    stopifnot(nrow(d_true) == nrow(d_sig))
    stopifnot(all(d_true$cluster == d_sig$cluster))
    stopifnot(nrow(d_sig) == nrow(d_heatmap))
    
    # identify clusters containing significant proportion of spike-in cells
    cutoff_prop <- 0.1
    d_true$spikein <- as.numeric(d_true$prop_spikein > cutoff_prop)
    
    
    # (iii) select clusters of interest: significant differential and true spike-in
    
    clus_keep <- (d_sig$sig == 1) | (d_true$spikein == 1)
    
    # skip this iteration if there are no clusters of interest
    if (sum(clus_keep) == 0) {
      plots_heatmaps[[ix]] <- textGrob("NA")
      plots_heights[[ix]] <- 1.5
      next
    }
    
    d_heatmap <- d_heatmap[clus_keep, , drop = FALSE]
    d_sig <- d_sig[clus_keep, , drop = FALSE]
    d_true <- d_true[clus_keep, , drop = FALSE]
    
    
    # (iv) row annotation: significant and true spike-in clusters
    
    d_annot <- data.frame(significant = factor(d_sig$sig, levels = c(0, 1), labels = c("no", "yes")), 
                          spikein = factor(d_true$spikein, levels = c(0, 1), labels = c("no", "yes")))
    
    ha_sig <- rowAnnotation(df = d_annot, 
                            col = list(significant = c("no" = "gray90", "yes" = "red"), 
                                       spikein = c("no" = "gray90", "yes" = "purple4")), 
                            annotation_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
                            width = unit(0.75, "cm"))
    
    
    # (v) row annotation: number of cells per cluster
    
    cnd_which <- colData(d_counts)$group %in% c("healthy", cond_names[j])
    cnd_healthy <- colData(d_counts)$group[cnd_which] == "healthy"
    cnd_j <- colData(d_counts)$group[cnd_which] == cond_names[j]
    
    d_abundance <- assay(d_counts)[clus_keep, cnd_which, drop = FALSE]
    
    col <- rep("navy", sum(cnd_which))
    col[cnd_j] <- "deepskyblue"
    
    # see examples in ComplexHeatmap package vignette 'Heatmap Annotations' on Bioconductor website
    func_anno <- function(index) {
      n <- length(index)
      pushViewport(viewport(xscale = range(d_abundance) + c(-1 * max(d_abundance) / 10, max(d_abundance) / 8), 
                            yscale = c(0.5, n + 0.5)))
      for (i in 1:ncol(d_abundance)) {
        # note 'index' argument matches row order to row clustering order
        grid.points(d_abundance[index, i], rev(seq_along(index)), 
                    pch = 16, size = unit(0.65, "char"), gp = gpar(col = col[i], alpha = 0.75))
      }
      grid.xaxis(at = round(seq(min(d_abundance), max(d_abundance), length.out = 4), digits = 0), 
                 gp = gpar(fontsize = 12, col = "gray25"))
      upViewport()
    }
    
    ha_abundance <- rowAnnotation(counts = func_anno, show_annotation_name = TRUE, 
                                  annotation_name_rot = 0, annotation_name_offset = unit(1.4, "cm"), 
                                  annotation_name_gp = gpar(fontsize = 14), 
                                  width = unit(3.5, "cm"))
    
    # legend
    lgd <- Legend(at = c("healthy", cond_names[j]), title = "group", type = "points", size = unit(0.65, "char"), 
                  title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12), 
                  legend_gp = gpar(col = c("navy", "deepskyblue")))
    
    
    # (vi) create heatmap
    
    # use 1% and 99% percentiles for color scale
    colors <- colorRamp2(quantile(d_heatmap, c(0.01, 0.5, 0.99)), 
                         c("dodgerblue", "white", "darkorange"))
    
    ht_title <- paste0("AML-sim, ", cond_names[j], ", threshold ", gsub("pc$", "\\%", thresholds[th]), ": diffcyt-DA-GLMM")
    
    ht <- Heatmap(d_heatmap, col = colors, name = "expression", 
                  row_title = "clusters", row_title_gp = gpar(fontsize = 14), 
                  column_title = "markers", column_title_side = "bottom", column_title_gp = gpar(fontsize = 14), 
                  column_names_gp = gpar(fontsize = 12), 
                  heatmap_legend_param = list(title_gp = gpar(fontface = "bold", fontsize = 12), labels_gp = gpar(fontsize = 12)), 
                  cluster_columns = FALSE, show_row_names = FALSE)
    
    
    # (vii) add annotation and save individual plot
    
    # save individual plot
    fn <- file.path(DIR_PLOTS, paste0("panels/results_AML_sim_diffcyt_DA_GLMM_main_heatmap_sub_AML_sim_", thresholds[th], "_", cond_names[j], ".pdf"))
    pdf(fn, width = 8.5, height = max(2.5, 1.5 + 0.25 * nrow(d_abundance)))
    draw(ht + ha_sig + ha_abundance, 
         annotation_legend_list = list(lgd), 
         column_title = ht_title, column_title_gp = gpar(fontface = "bold", fontsize = 16), 
         newpage = FALSE)
    dev.off()
    
    
    # store plot object for multi-panel plot
    plots_heatmaps[[ix]] <- grid.grabExpr(draw(ht + ha_sig + ha_abundance, 
                                               annotation_legend_list = list(lgd), 
                                               column_title = ht_title, column_title_gp = gpar(fontface = "bold", fontsize = 16), 
                                               newpage = FALSE))
    
    # store plot heights for multi-panel plot
    plots_heights[[ix]] <- max(1.8, 1.6 + 0.19 * nrow(d_abundance))
    
  }
}




########################
# Save multi-panel plots
########################

heights <- sapply(split(unlist(plots_heights), rep(1:length(thresholds), each = length(cond_names))), max)
widths <- rep(7.7, length(cond_names))

fn <- file.path(DIR_PLOTS, paste0("results_AML_sim_diffcyt_DA_GLMM_main_heatmaps_sub.pdf"))

pdf(fn, width = 16, height = sum(heights) + 0.5)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = 2, just = "top", 
                                           heights = unit(heights, "inches"), widths = unit(widths, "inches"))))

for (th in 1:length(thresholds)) {
  for (j in 1:length(cond_names)) {
    ix <- (th * length(cond_names)) - (length(cond_names) - j)
    
    pushViewport(viewport(layout.pos.row = th, layout.pos.col = j))
    grid.draw(plots_heatmaps[[ix]])
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



