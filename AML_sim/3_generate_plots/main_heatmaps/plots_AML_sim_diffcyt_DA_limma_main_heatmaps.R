##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: heatmaps
# - method: diffcyt-DA-limma
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

load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_limma_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_limma_main.RData"))
load(file.path(DIR_RDATA, "out_objects_AML_sim_diffcyt_DA_limma_main.RData"))


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
plots_heatmaps <- vector("list", length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  # load plot data objects (same for both conditions j)
  d_se <- out_objects_diffcyt_DA_limma_main[[th]]$d_se
  d_counts <- out_objects_diffcyt_DA_limma_main[[th]]$d_counts
  d_medians_all <- out_objects_diffcyt_DA_limma_main[[th]]$d_medians_all
  
  
  # -----------------------------------------------------------
  # create heatmap: main panel showing marker expression values
  # -----------------------------------------------------------
  
  # note: using clustering markers only
  # note: no additional scaling (using asinh-transformed values directly)
  d_heatmap <- assay(d_medians_all)[, colData(d_medians_all)$is_clustering_col]
  
  colnames(d_heatmap) <- gsub("\\(.*$", "", colnames(d_heatmap))
  
  # use 1% and 99% percentiles for color scale
  colors <- colorRamp2(quantile(d_heatmap, c(0.01, 0.5, 0.99)), 
                       c("dodgerblue", "white", "darkorange"))
  
  ht <- Heatmap(d_heatmap, col = colors, 
                name = "expression", row_title = "clusters", 
                column_title = "markers", column_title_side = "bottom", 
                cluster_columns = FALSE, show_row_names = FALSE)
  
  
  for (j in 1:length(cond_names)) {
    
    # index to store plots sequentially in list
    ix <- (th * length(cond_names)) - (length(cond_names) - j)
    
    
    # -------------------------------------------------------------
    # row annotation for significant and true differential clusters
    # -------------------------------------------------------------
    
    # (i) from cluster-level results
    
    # load cluster-level results (for condition j)
    d_clus <- out_clusters_diffcyt_DA_limma_main[[th]][[j]]
    stopifnot(nrow(d_clus) == nrow(rowData(d_counts)))
    stopifnot(all(d_clus$cluster == rowData(d_counts)$cluster))
    
    # significant differential clusters
    cutoff_sig <- 0.1
    sig <- d_clus$FDR <= cutoff_sig
    # set filtered clusters to FALSE
    sig[is.na(sig)] <- FALSE
    
    # set up data frame
    d_sig <- data.frame(cluster = rowData(d_counts)$cluster, 
                        sig = as.numeric(sig), 
                        n_cells = rowData(d_counts)$n_cells)
    
    
    # (ii) from cell-level results
    
    # load spike-in status at cell level (for condition j)
    spikein <- out_diffcyt_DA_limma_main[[th]][[j]]$spikein
    
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
    
    
    # (iii) add row annotation and title
    
    d_annot <- data.frame(significant = factor(d_sig$sig, levels = c(0, 1), labels = c("no", "yes")), 
                          spikein = factor(d_true$spikein, levels = c(0, 1), labels = c("no", "yes")))
    
    ha_bar <- rowAnnotation(df = d_annot, 
                            col = list(significant = c("no" = "gray90", "yes" = "red"), 
                                       spikein = c("no" = "gray90", "yes" = "black")), 
                            width = unit(1, "cm"))
    
    ht_title <- paste0("AML-sim, ", cond_names[j], ", threshold ", gsub("pc$", "\\%", thresholds[th]), ": diffcyt-DA-limma")
    
    
    # (iv) save plot
    
    fn <- file.path(DIR_PLOTS, paste0("results_diffcyt_DA_limma_main_heatmap_AML_sim_", thresholds[th], "_", cond_names[j], ".pdf"))
    pdf(fn, width = 7, height = 5.5)
    ht_list <- draw(ht + ha_bar, column_title = ht_title)
    dev.off()
    
  }
}




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



