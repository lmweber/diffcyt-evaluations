##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: ROC curves and pAUC values
# - method: diffcyt-DA-edgeR
# 
# - supplementary results: varying clustering resolution
# 
# Lukas Weber, February 2018
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2
library(viridis)


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/supp_clustering_resolution"

load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_edgeR_supp_clustering_resolution.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/supp_clustering_resolution"




################
# Generate plots
################

# loop over thresholds (th) and conditions (j)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# condition names
cond_names <- c("CN", "CBF")


# store plots in list
plots_ROC <- plots_pAUC <- vector("list", length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  for (j in 1:length(cond_names)) {
    
    
    # index to store plots in list
    ix <- (j * length(thresholds)) - (length(thresholds) - th)
    
    
    # -------------------------------------
    # Pre-processing steps for iCOBRA plots
    # -------------------------------------
    
    # create 'COBRAData' object
    data <- list(k_9   = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_9"]][[th]][[j]], 
                 k_25  = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_25"]][[th]][[j]], 
                 k_49  = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_49"]][[th]][[j]], 
                 k_100 = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_100"]][[th]][[j]], 
                 k_196 = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_196"]][[th]][[j]], 
                 k_400 = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_400"]][[th]][[j]], 
                 k_900 = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_900"]][[th]][[j]], 
                 k_1600 = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_1600"]][[th]][[j]])
    
    # check
    stopifnot(all(sapply(data, function(d) all(d$spikein == data[[1]]$spikein))))
    
    # note: provide all available values
    # 'padj' is required for threshold points on TPR-FDR curves
    # depending on availability, plotting functions use 'score', then 'pval', then 'padj'
    cobradata <- COBRAData(pval = data.frame(k_9   = data[["k_9"]][, "p_vals"], 
                                             k_25  = data[["k_25"]][, "p_vals"], 
                                             k_49  = data[["k_49"]][, "p_vals"], 
                                             k_100 = data[["k_100"]][, "p_vals"], 
                                             k_196 = data[["k_196"]][, "p_vals"], 
                                             k_400 = data[["k_400"]][, "p_vals"], 
                                             k_900 = data[["k_900"]][, "p_vals"], 
                                             k_1600 = data[["k_1600"]][, "p_vals"]), 
                           padj = data.frame(k_9   = data[["k_9"]][, "p_adj"], 
                                             k_25  = data[["k_25"]][, "p_adj"], 
                                             k_49  = data[["k_49"]][, "p_adj"], 
                                             k_100 = data[["k_100"]][, "p_adj"], 
                                             k_196 = data[["k_196"]][, "p_adj"], 
                                             k_400 = data[["k_400"]][, "p_adj"], 
                                             k_900 = data[["k_900"]][, "p_adj"], 
                                             k_1600 = data[["k_1600"]][, "p_adj"]), 
                           truth = data.frame(spikein = data[["k_9"]][, "spikein"]))
    
    # calculate performance scores
    # (note: can ignore warning messages when 'padj' not available)
    cobraperf <- calculate_performance(cobradata, 
                                       binary_truth = "spikein", 
                                       aspects = "roc")
    
    # color scheme
    colors <- rev(viridis(8))
    
    colors <- colors[1:length(data)]
    names(colors) <- names(data)
    
    # prepare plotting object
    cobraplot <- prepare_data_for_plot(cobraperf, 
                                       colorscheme = colors, 
                                       conditionalfill = FALSE)
    
    # re-order legend
    cobraplot <- reorder_levels(cobraplot, levels = names(data))
    
    
    
    # ----------
    # ROC curves
    # ----------
    
    # create plot
    p_ROC <- 
      plot_roc(cobraplot, linewidth = 0.75) + 
      scale_color_manual(labels = gsub("k_", "", names(colors)), values = colors) + 
      coord_fixed() + 
      xlab("False positive rate") + 
      ylab("True positive rate") + 
      ggtitle(paste0("AML-sim: clustering resolution, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th])), subtitle = "ROC curves, diffcyt-DA-edgeR") + 
      theme_bw() + 
      theme(strip.text.x = element_blank()) + 
      guides(color = guide_legend("no. clusters"))
    
    plots_ROC[[ix]] <- p_ROC
    
    # save individual panel plot
    fn <- file.path(DIR_PLOTS, "panels", 
                    paste0("results_AML_sim_diffcyt_DA_edgeR_supp_clustering_resolution_ROC_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 4.75, height = 3.5)
    
    
    
    # -----------
    # pAUC values
    # -----------
    
    # calculate pAUC values
    pAUC <- rep(NA, length(data))
    names(pAUC) <- names(data)
    
    thresh <- 0.2
    
    for (i in 1:length(pAUC)) {
      roc_i <- roc(cobraperf)[roc(cobraperf)$method == names(pAUC)[i], ]
      roc_i_less <- roc_i[roc_i$FPR < thresh, ]
      roc_i_more <- roc_i[roc_i$FPR >= thresh, ]
      
      row_insert <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(roc_i)))
      colnames(row_insert) <- colnames(roc_i)
      row_insert[1, "FPR"] <- thresh
      row_insert[1, "TPR"] <- roc_i_less[nrow(roc_i_less), "TPR"]
      
      roc_i <- rbind(roc_i_less, row_insert, roc_i_more)
      rownames(roc_i) <- NULL
      
      #roc_i$FPR <- c(0, FPR[-1])
      roc_i$dFPR <- c(0, diff(roc_i$FPR))
      roc_i$dTPR <- c(0, diff(roc_i$TPR))
      roc_i$TPRs <- c(0, roc_i$TPR[-length(roc_i$TPR)])
      roc_i$AUC <- cumsum(roc_i$dFPR * roc_i$dTPR/2 + roc_i$dFPR * roc_i$TPRs)
      
      # extract pAUC
      pAUC[i] <- roc_i$AUC[roc_i$FPR == thresh]
    }
    
    # normalize pAUC values
    pAUC <- pAUC / thresh
    
    # create data frame for plotting
    names(pAUC) <- gsub("k_", "", names(pAUC))
    p_data <- as.data.frame(pAUC)
    p_data$n_clusters <- factor(names(pAUC), levels = names(pAUC))
    p_data$group <- "pAUC"
    
    names(colors) <- gsub("k_", "", names(colors))
    
    
    # create plot
    p_pAUC <- 
      ggplot(p_data, aes(x = n_clusters, y = pAUC, group = group, color = n_clusters)) + 
      geom_line(lwd = 0.7) + 
      geom_point(size = 2.25) + 
      scale_color_manual(values = colors) + 
      ylim(c(0, 1)) + 
      xlab("Number of clusters") + 
      ylab(paste0("pAUC (FPR < ", thresh, ")")) + 
      ggtitle(paste0("AML-sim: clustering resolution, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th])), subtitle = "pAUC, diffcyt-DA-edgeR") + 
      theme_bw() + 
      guides(color = guide_legend("no. clusters"))
    
    plots_pAUC[[ix]] <- p_pAUC
    
    # save individual panel plot
    fn <- file.path(DIR_PLOTS, "panels", 
                    paste0("results_AML_sim_diffcyt_DA_edgeR_supp_clustering_resolution_pAUC_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 4.75, height = 3.5)
    
  }
}




########################
# Save multi-panel plots
########################

# ----------
# ROC curves
# ----------

# modify plot elements
plots_multi <- lapply(plots_ROC, function(p) {
  p + 
    labs(title = gsub("^.*resolution, ", "", p$labels$title)) + 
    theme(legend.position = "none")
})

# format into grid
plots_multi <- plot_grid(plotlist = plots_multi, nrow = 2, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_multi <- ggdraw() + draw_label(gsub(",.*$", ", diffcyt-DA-edgeR", plots_ROC[[1]]$labels$title), fontface = "bold")
grid_multi <- plot_grid(title_multi, plots_multi, ncol = 1, rel_heights = c(1, 20))

# add combined legend
legend_multi <- get_legend(plots_ROC[[1]] + theme(legend.position = "right", 
                                                   legend.title = element_text(size = 12, face = "bold"), 
                                                   legend.text = element_text(size = 12)))
grid_multi <- plot_grid(grid_multi, legend_multi, nrow = 1, rel_widths = c(5.5, 1))

# save plots
fn_multi <- file.path(DIR_PLOTS, paste0("results_AML_sim_diffcyt_DA_edgeR_supp_clustering_resolution_ROC.pdf"))
ggsave(fn_multi, grid_multi, width = 9.5, height = 6)



# -----------
# pAUC values
# -----------

# modify plot elements
plots_multi <- lapply(plots_pAUC, function(p) {
  p + 
    labs(title = gsub("^.*resolution, ", "", p$labels$title)) + 
    theme(legend.position = "none")
})

# format into grid
plots_multi <- plot_grid(plotlist = plots_multi, nrow = 2, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_multi <- ggdraw() + draw_label(gsub(",.*$", ", diffcyt-DA-edgeR", plots_pAUC[[1]]$labels$title), fontface = "bold")
grid_multi <- plot_grid(title_multi, plots_multi, ncol = 1, rel_heights = c(1, 20))

# add combined legend
legend_multi <- get_legend(plots_pAUC[[1]] + theme(legend.position = "right", 
                                                  legend.title = element_text(size = 12, face = "bold"), 
                                                  legend.text = element_text(size = 12)))
grid_multi <- plot_grid(grid_multi, legend_multi, nrow = 1, rel_widths = c(5.5, 1))

# save plots
fn_multi <- file.path(DIR_PLOTS, paste0("results_AML_sim_diffcyt_DA_edgeR_supp_clustering_resolution_pAUC.pdf"))
ggsave(fn_multi, grid_multi, width = 9.5, height = 6)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



