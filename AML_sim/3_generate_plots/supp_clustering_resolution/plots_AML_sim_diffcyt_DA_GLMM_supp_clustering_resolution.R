##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: ROC curves and pAUC values
# - method: diffcyt-DA-GLMM
# 
# - supplementary results: clustering resolution
# 
# Lukas Weber, October 2017
##########################################################################################


library(iCOBRA)
library(dplyr)
library(magrittr)
library(reshape2)
library(ggplot2)
library(cowplot)


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/supp_clustering_resolution"

load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_GLMM_supp_clustering_resolution.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/supp_clustering_resolution"




############################
# Generate plots: ROC curves
############################

# one panel per threshold (th) and condition (j)


# clustering resolution
resolution <- c(3, 5, 10, 20, 30, 40, 50)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")


# store plots in list
plots_ROC <- vector("list", length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  for (j in 1:length(cond_names)) {
    
    # -------------------------------------
    # Pre-processing steps for iCOBRA plots
    # -------------------------------------
    
    # index to store plots sequentially in list
    ix <- (th * length(cond_names)) - (length(cond_names) - j)
    
    # create 'COBRAData' object
    data <- list(k_9    = out_diffcyt_DA_GLMM_supp_resolution[[1]][[th]][[j]], 
                 k_25   = out_diffcyt_DA_GLMM_supp_resolution[[2]][[th]][[j]], 
                 k_100  = out_diffcyt_DA_GLMM_supp_resolution[[3]][[th]][[j]], 
                 k_400  = out_diffcyt_DA_GLMM_supp_resolution[[4]][[th]][[j]], 
                 k_900  = out_diffcyt_DA_GLMM_supp_resolution[[5]][[th]][[j]], 
                 k_1600 = out_diffcyt_DA_GLMM_supp_resolution[[6]][[th]][[j]], 
                 k_2500 = out_diffcyt_DA_GLMM_supp_resolution[[7]][[th]][[j]])
    
    # check
    stopifnot(all(sapply(data, function(d) all(d$spikein == data[[1]]$spikein))))
    
    # note: provide all available values
    # - 'padj' is required for threshold points on TPR-FDR curves
    # - depending on availability, plotting functions use 'score', then 'pval', then 'padj'
    cobradata <- COBRAData(pval = data.frame(k_9    = data[["k_9"]][, "p_vals"], 
                                             k_25   = data[["k_25"]][, "p_vals"], 
                                             k_100  = data[["k_100"]][, "p_vals"], 
                                             k_400  = data[["k_400"]][, "p_vals"], 
                                             k_900  = data[["k_900"]][, "p_vals"], 
                                             k_1600 = data[["k_1600"]][, "p_vals"], 
                                             k_2500 = data[["k_2500"]][, "p_vals"]), 
                           padj = data.frame(k_9    = data[["k_9"]][, "p_adj"], 
                                             k_25   = data[["k_25"]][, "p_adj"], 
                                             k_100  = data[["k_100"]][, "p_adj"], 
                                             k_400  = data[["k_400"]][, "p_adj"], 
                                             k_900  = data[["k_900"]][, "p_adj"], 
                                             k_1600 = data[["k_1600"]][, "p_adj"], 
                                             k_2500 = data[["k_2500"]][, "p_adj"]), 
                           truth = data.frame(spikein = data[["k_9"]][, "spikein"]))
    
    # calculate performance scores
    # (note: can ignore warning messages when 'padj' not available)
    cobraperf <- calculate_performance(cobradata, 
                                       binary_truth = "spikein", 
                                       aspects = "roc")
    
    # color scheme
    colors <- colorRampPalette(c("#c7e9c0", "darkgreen"))(7)
    
    colors <- colors[1:length(data)]
    names(colors) <- names(data)
    
    # axis ranges
    x_range <- c(0, 1)
    y_range <- c(0, 1)
    
    # prepare plotting object
    cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colors)
    
    # re-order legend
    cobraplot <- reorder_levels(cobraplot, levels = names(data))
    
    
    # ----------
    # ROC curves
    # ----------
    
    # create plot
    p <- plot_roc(cobraplot, linewidth = 0.75)
    
    # with short title for multi-panel plot
    p <- p + 
      coord_fixed(xlim = x_range, ylim = y_range) + 
      xlab("False positive rate") + 
      ylab("True positive rate") + 
      ggtitle(paste0(cond_names[j], ", threshold ", gsub("pc$", "\\%", thresholds[th]))) + 
      theme_bw() + 
      theme(strip.text.x = element_blank()) + 
      guides(color = guide_legend("n_clusters"))
    
    plots_ROC[[ix]] <- p
    
    # save individual panel plot
    p <- p + 
      ggtitle(paste0("AML-sim, clustering resolution: diffcyt-DA-GLMM: ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th]), ": ROC curves"))
    
    fn <- file.path(DIR_PLOTS, "panels", 
                    paste0("results_AML_sim_diffcyt_DA_GLMM_supp_clustering_resolution_ROC_curves_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 7, height = 5.5)
    
  }
}




#############################
# Generate plots: pAUC values
#############################

# one panel per condition (j)


# clustering resolution
resolution <- c(3, 5, 10, 20, 30, 40, 50)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")


# store plots in list
plots_pAUC <- vector("list", length(cond_names))


for (j in 1:length(cond_names)) {
  
  # -------------------------------------
  # Pre-processing steps for iCOBRA plots
  # -------------------------------------
  
  # loop over thresholds (to create one line per threshold)
  
  pAUC_list <- vector("list", length(thresholds))
  names(pAUC_list) <- thresholds
  
  for (th in 1:length(thresholds)) {
    
    # create 'COBRAData' object
    data <- list(k_9    = out_diffcyt_DA_GLMM_supp_resolution[[1]][[th]][[j]], 
                 k_25   = out_diffcyt_DA_GLMM_supp_resolution[[2]][[th]][[j]], 
                 k_100  = out_diffcyt_DA_GLMM_supp_resolution[[3]][[th]][[j]], 
                 k_400  = out_diffcyt_DA_GLMM_supp_resolution[[4]][[th]][[j]], 
                 k_900  = out_diffcyt_DA_GLMM_supp_resolution[[5]][[th]][[j]], 
                 k_1600 = out_diffcyt_DA_GLMM_supp_resolution[[6]][[th]][[j]], 
                 k_2500 = out_diffcyt_DA_GLMM_supp_resolution[[7]][[th]][[j]])
    
    # check
    stopifnot(all(sapply(data, function(d) all(d$spikein == data[[1]]$spikein))))
    
    # note: provide all available values
    # - 'padj' is required for threshold points on TPR-FDR curves
    # - depending on availability, plotting functions use 'score', then 'pval', then 'padj'
    cobradata <- COBRAData(pval = data.frame(k_9    = data[["k_9"]][, "p_vals"], 
                                             k_25   = data[["k_25"]][, "p_vals"], 
                                             k_100  = data[["k_100"]][, "p_vals"], 
                                             k_400  = data[["k_400"]][, "p_vals"], 
                                             k_900  = data[["k_900"]][, "p_vals"], 
                                             k_1600 = data[["k_1600"]][, "p_vals"], 
                                             k_2500 = data[["k_2500"]][, "p_vals"]), 
                           padj = data.frame(k_9    = data[["k_9"]][, "p_adj"], 
                                             k_25   = data[["k_25"]][, "p_adj"], 
                                             k_100  = data[["k_100"]][, "p_adj"], 
                                             k_400  = data[["k_400"]][, "p_adj"], 
                                             k_900  = data[["k_900"]][, "p_adj"], 
                                             k_1600 = data[["k_1600"]][, "p_adj"], 
                                             k_2500 = data[["k_2500"]][, "p_adj"]), 
                           truth = data.frame(spikein = data[["k_9"]][, "spikein"]))
    
    # calculate performance scores
    # (note: can ignore warning messages when 'padj' not available)
    cobraperf <- calculate_performance(cobradata, 
                                       binary_truth = "spikein", 
                                       aspects = "roc")
    
    # calculate pAUC values (store in separate vector)
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
    
    # calculate AUC values (store in main object) (code from Charlotte)
    # roc(cobraperf) <- 
    #   roc(cobraperf) %>% 
    #   dplyr::group_by(method) %>% 
    #   dplyr::mutate(FPR = c(0, FPR[-1])) %>% 
    #   dplyr::mutate(dFPR = c(0, diff(FPR)), 
    #                 dTPR = c(0, diff(TPR)), 
    #                 TPRs = c(0, TPR[-length(TPR)])) %>% 
    #   dplyr::mutate(AUC = cumsum(dFPR * dTPR/2 + dFPR * TPRs)) %>% 
    #   dplyr::ungroup() %>% 
    #   as.data.frame()
    
    # store pAUC values
    pAUC_list[[th]] <- pAUC
  }
  
  # create data frame for plotting
  for (z in 1:length(pAUC_list)) stopifnot(all(names(pAUC_list[[z]]) == names(pAUC_list[[1]])))
  
  p_data <- as.data.frame(do.call("cbind", pAUC_list))
  p_data$k_clusters <- factor(names(pAUC_list[[1]]), levels = names(pAUC_list[[1]]))
  p_data <- melt(p_data, id.vars = "k_clusters", variable.name = "threshold", value.name = "pAUC")
  
  # color scheme
  colors <- colorRampPalette(c("#023858", "#d0d1e6"))(4)
  
  colors <- colors[1:length(pAUC_list)]
  names(colors) <- names(pAUC_list)
  
  # axis ranges
  y_range_pAUC <- c(0, 1)
  
  
  # ----------
  # pAUC plots
  # ----------
  
  # create plot
  p <- 
    ggplot(p_data, aes(x = k_clusters, y = pAUC, group = threshold, color = threshold)) + 
    geom_line(lwd = 0.7) + 
    geom_point(size = 2.25) + 
    scale_color_manual(values = colors) + 
    xlab("Number of clusters") + 
    ylim(y_range_pAUC) + 
    ylab(paste0("pAUC (FPR < ", thresh, ")")) + 
    theme_bw() + 
    ggtitle(cond_names[j])
  
  plots_pAUC[[j]] <- p
  
  # save individual panel plot
  p <- p + 
    ggtitle(paste0("AML-sim, clustering resolution: diffcyt-DA-GLMM: ", cond_names[j], ": pAUC (FPR < ", thresh, ")"))
  
  fn <- file.path(DIR_PLOTS, "panels", 
                  paste0("results_AML_sim_diffcyt_DA_GLMM_supp_clustering_resolution_pAUC_", cond_names[j], ".pdf"))
  ggsave(fn, width = 7, height = 5.5)
  
}




########################
# Save multi-panel plots
########################

# ----------
# ROC curves
# ----------

# re-order plots to fill each condition by row
ord <- c(2 * (1:4) - 1, 2 * (1:4))
plots_ROC <- plots_ROC[ord]

# modify plot elements
plots_ROC <- lapply(plots_ROC, function(p) {
  p + theme(legend.position = "none", 
            axis.title.x = element_blank(), axis.title.y = element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

# format into grid
grid_ROC <- do.call(plot_grid, append(plots_ROC, list(nrow = 2, ncol = 4, 
                                                      align = "hv", axis = "bl", 
                                                      scale = 0.975)))

# add combined axis titles
xaxis_ROC <- ggdraw() + draw_label("False positive rate", size = 14)
yaxis_ROC <- ggdraw() + draw_label("True positive rate", size = 14, angle = 90)

grid_ROC <- plot_grid(grid_ROC, xaxis_ROC, ncol = 1, rel_heights = c(12, 1))
grid_ROC <- plot_grid(yaxis_ROC, grid_ROC, nrow = 1, rel_widths = c(1, 30))

# add combined title
title_ROC <- ggdraw() + draw_label("AML-sim, clustering resolution: diffcyt-DA-GLMM: ROC curves", fontface = "bold")
grid_ROC <- plot_grid(title_ROC, grid_ROC, ncol = 1, rel_heights = c(1, 13))

# add combined legend
legend_ROC <- get_legend(plots_ROC[[1]] + theme(legend.position = "right", 
                                                legend.title = element_text(size = 12, face = "bold"), 
                                                legend.text = element_text(size = 12)))
grid_ROC <- plot_grid(grid_ROC, legend_ROC, nrow = 1, rel_widths = c(4.8, 1))

# save plots
fn_ROC <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_DA_GLMM_supp_clustering_resolution_ROC_curves.pdf")
ggsave(fn_ROC, grid_ROC, width = 11, height = 5.4)



# -----------
# pAUC values
# -----------

# modify plot elements
plots_pAUC <- lapply(plots_pAUC, function(p) {
  p + theme(legend.position = "none", 
            axis.title.x = element_blank(), axis.title.y = element_blank())
})

# format into grid
grid_pAUC <- do.call(plot_grid, append(plots_pAUC, list(nrow = 1, ncol = 2, 
                                                        align = "h", axis = "b", 
                                                        scale = 0.95)))

# add combined axis titles
xaxis_pAUC <- ggdraw() + draw_label("Number of clusters", size = 12)
yaxis_pAUC <- ggdraw() + draw_label(paste0("pAUC (FPR < ", thresh, ")"), size = 12, angle = 90)

grid_pAUC <- plot_grid(grid_pAUC, xaxis_pAUC, ncol = 1, rel_heights = c(15, 1))
grid_pAUC <- plot_grid(yaxis_pAUC, grid_pAUC, nrow = 1, rel_widths = c(1, 30))

# add combined title
title_pAUC <- ggdraw() + draw_label(paste0("AML-sim, clustering resolution: diffcyt-DA-GLMM: pAUC (FPR < ", thresh, ")"), fontface = "bold")
grid_pAUC <- plot_grid(title_pAUC, grid_pAUC, ncol = 1, rel_heights = c(1, 12))

# add combined legend
legend_pAUC <- get_legend(plots_pAUC[[1]] + theme(legend.position = "right", 
                                                  legend.title = element_text(size = 11, face = "bold"), 
                                                  legend.text = element_text(size = 11)))
grid_pAUC <- plot_grid(grid_pAUC, legend_pAUC, nrow = 1, rel_widths = c(10, 1))

# save plots
fn_pAUC <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_DA_GLMM_supp_clustering_resolution_pAUC.pdf")
ggsave(fn_pAUC, grid_pAUC, width = 9, height = 3.6)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



