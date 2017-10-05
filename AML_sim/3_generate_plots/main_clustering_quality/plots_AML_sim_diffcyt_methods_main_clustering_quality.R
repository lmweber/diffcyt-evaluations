##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: clustering quality (bar plots)
# - method: diffcyt methods
# 
# - main results
# 
# Lukas Weber, October 2017
##########################################################################################


# note: clustering step is the same for all three 'diffcyt' methods (diffcyt-DA-edgeR,
# diffcyt-DA-GLMM, diffcyt-DA-limma)


library(SummarizedExperiment)
library(ggplot2)
library(cowplot)


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/main"

# note: only need to load one set of results, since clustering step is the same for all 'diffcyt' methods
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA, "out_objects_AML_sim_diffcyt_DA_edgeR_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/main_clustering_quality"




##################################
# Calculate clustering performance
##################################

# loop over thresholds (th) and conditions (j)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")

# lists to store clustering performance results
clustering_pr <- clustering_re <- clustering_F1 <- vector("list", length(thresholds))
names(clustering_pr) <- names(clustering_re) <- names(clustering_F1) <- thresholds


for (th in 1:length(thresholds)) {
  
  clustering_pr[[th]] <- clustering_re[[th]] <- clustering_F1[[th]] <- vector("list", length(cond_names))
  names(clustering_pr[[th]]) <- names(clustering_re[[th]]) <- names(clustering_F1[[th]]) <- cond_names
  
  
  for (j in 1:length(cond_names)) {
    
    # -------------------------------------
    # load data objects and spike-in status
    # -------------------------------------
    
    # load data objects
    # note: clustering is performed once on all samples from both conditions together
    d_se <- out_objects_diffcyt_DA_edgeR_main[[th]]$d_se
    
    group_IDs <- unname(sapply(split(rowData(d_se)$group, rowData(d_se)$sample), unique))
    
    # load spike-in status for condition j
    spikein <- out_diffcyt_DA_edgeR_main[[th]][[j]]$spikein
    
    # add spike-in status to data object
    rowData(d_se)$spikein <- 0
    rowData(d_se)$spikein[rowData(d_se)$group == cond_names[j]] <- spikein
    
    
    # --------------------------------------------------------------------------------
    # calculate clustering performance for all clusters containing true spike-in cells
    # --------------------------------------------------------------------------------
    
    # find all matching clusters (clusters containing true spike-in cells) for each sample
    
    d_split <- split(rowData(d_se), rowData(d_se)$sample)
    
    # check no missing clusters
    stopifnot(all(sapply(d_split, function(d) all(names(table(d[d$spikein == 1, ]$cluster)) == levels(rowData(d_se)$cluster)))))
    
    # find matching clusters
    labels_matched <- lapply(d_split, function(d) which(table(d[d$spikein == 1, ]$cluster) > 0))
    
    labels_matched[(group_IDs == "healthy") | !(group_IDs == cond_names[j])] <- NA
    labels_matched
    
    # total number of cells in each matching cluster
    
    n_matched_samp <- mapply(function(d, l) {
      sapply(l, function(l) sum(d$cluster == l))
    }, d_split, labels_matched)
    n_matched_samp
    
    # number of true spike-in cells in each matching cluster
    
    n_matched_samp_correct <- mapply(function(d, l) {
      sapply(l, function(l) sum(d$cluster == l & d$spikein == 1))
    }, d_split, labels_matched)
    n_matched_samp_correct
    
    # total number of spike-in cells per sample
    
    n_spikein_samp <- sapply(d_split, function(d) sum(d$spikein == 1))
    n_spikein_samp
    
    # calculate precision, recall, F1 score for each matching cluster
    
    pr <- mapply(function(n_smp, n_cor) {
      n_cor / n_smp
    }, n_matched_samp, n_matched_samp_correct)
    
    re <- mapply(function(n_spk, n_cor) {
      n_cor / n_spk
    }, n_spikein_samp, n_matched_samp_correct)
    
    F1 <- mapply(function(pr, re) {
      2 * (pr * re) / (pr + re)
    }, pr, re)
    
    pr[(group_IDs == "healthy") | !(group_IDs == cond_names[j])] <- NA
    re[(group_IDs == "healthy") | !(group_IDs == cond_names[j])] <- NA
    F1[(group_IDs == "healthy") | !(group_IDs == cond_names[j])] <- NA
    
    # store results
    clustering_pr[[th]][[j]] <- pr
    clustering_re[[th]][[j]] <- re
    clustering_F1[[th]][[j]] <- F1
    
  }
}




################
# Generate plots
################

# loop over thresholds (th) and conditions (j)

# store plots in list
plots_clustering <- vector("list", length(thresholds) * length(cond_names))

plot_widths <- rep(NA, length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  for (j in 1:length(cond_names)) {
    
    # index to store plots sequentially in list
    ix <- (th * length(cond_names)) - (length(cond_names) - j)
    
    
    # create data frame for plotting
    
    pr <- clustering_pr[[th]][[j]][group_IDs == cond_names[j]]
    re <- clustering_re[[th]][[j]][group_IDs == cond_names[j]]
    F1 <- clustering_F1[[th]][[j]][group_IDs == cond_names[j]]
    
    pr_all <- lapply(pr, function(vals) data.frame(cluster = names(vals), value = vals))
    re_all <- lapply(re, function(vals) data.frame(cluster = names(vals), value = vals))
    F1_all <- lapply(F1, function(vals) data.frame(cluster = names(vals), value = vals))
    
    for (i in 1:length(pr_all)) {
      pr_all[[i]] <- cbind(pr_all[[i]], sample = gsub("^.*_", "", names(pr_all)[i]), measure = "precision")
      re_all[[i]] <- cbind(re_all[[i]], sample = gsub("^.*_", "", names(re_all)[i]), measure = "recall")
      F1_all[[i]] <- cbind(F1_all[[i]], sample = gsub("^.*_", "", names(F1_all)[i]), measure = "F1_score")
    }
    
    pr_all <- do.call("rbind", pr_all)
    re_all <- do.call("rbind", re_all)
    F1_all <- do.call("rbind", F1_all)
    
    d_plot <- rbind(pr_all, re_all, F1_all)
    
    
    # vector of clusters for x-axis labels
    
    # sort clusters by sum of F1 scores (use sum instead of average to account for missing clusters)
    sum_F1 <- sapply(split(F1_all, F1_all$cluster), function(s) sum(s$value))
    x_clusters_ord <- as.numeric(names(rev(sort(sum_F1))))
    
    # check
    x_clusters_check <- sort(unique(as.numeric(unname(unlist(sapply(pr, names))))))
    stopifnot(all(x_clusters_ord %in% x_clusters_check) & all(x_clusters_check %in% x_clusters_ord))
    
    
    d_plot$cluster <- factor(d_plot$cluster, levels = x_clusters_ord, labels = x_clusters_ord)
    
    
    # create plot
    
    plot_widths[ix] <- 2 + length(x_clusters_ord) / 7
    
    colors <- c("navy", "deepskyblue", "red")
    #colors <- c("#009E73", "#56B4E9", "#F8766D")
    
    shapes <- c(21, 24, 22, 3, 4)
    
    p <- ggplot(d_plot, aes(x = cluster, y = value, color = measure, shape = sample)) + 
      geom_point(stroke = 1) + 
      scale_color_manual(values = colors) + 
      scale_shape_manual(values = shapes) + 
      ggtitle(paste0(cond_names[j], ", threshold ", gsub("pc$", "\\%", thresholds[th]))) + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
            axis.title.y = element_blank())
    
    plots_clustering[[ix]] <- p
    
    # save individual panel plot
    fn <- file.path(DIR_PLOTS, "panels", 
                    paste0("results_AML_sim_diffcyt_methods_main_clustering_quality_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = plot_widths[ix], height = 3.5)
    
  }
}




########################
# Save multi-panel plots
########################

# re-order plots sequentially
ord <- c(2 * (1:4) - 1, 2 * (1:4))
plots_clustering <- plots_clustering[ord]
plot_widths <- plot_widths[ord]

# modify plot elements
plots_clustering <- lapply(plots_clustering, function(p) {
  p + theme(legend.position = "none")
})


for (j in 1:length(cond_names)) {
  
  # format into grid
  plots_rows <- vector("list", length(thresholds))
  for (th in 1:length(thresholds)) {
    # divide each row into 2 columns to allow varying widths
    ix_rows <- th + (j - 1) * length(thresholds)
    plots_rows[[th]] <- plot_grid(plots_clustering[[ix_rows]], NULL, 
                                  nrow = 1, ncol = 2, 
                                  rel_widths = c(plot_widths[ix_rows], max(plot_widths) - plot_widths[ix_rows]))
  }
  # combine rows
  grid_clustering <- do.call(plot_grid, append(plots_rows, list(nrow = 4, ncol = 1, scale = 0.975)))
  
  # add combined title
  title_clustering <- ggdraw() + draw_label("AML-sim, diffcyt methods: clustering quality measures", fontface = "bold")
  grid_clustering <- plot_grid(title_clustering, grid_clustering, ncol = 1, rel_heights = c(1, 25))
  
  # add combined legend
  legend_clustering <- get_legend(plots_clustering[[1]] + theme(legend.position = "right", 
                                                                legend.title = element_text(size = 12, face = "bold"), 
                                                                legend.text = element_text(size = 12)))
  grid_clustering <- plot_grid(grid_clustering, legend_clustering, nrow = 1, rel_widths = c(6, 1))
  
  # save plots
  fn_clustering <- file.path(DIR_PLOTS, paste0("results_AML_sim_diffcyt_methods_main_clustering_quality_", cond_names[j], ".pdf"))
  ggsave(fn_clustering, grid_clustering, width = 9, height = 11)
  
}




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



