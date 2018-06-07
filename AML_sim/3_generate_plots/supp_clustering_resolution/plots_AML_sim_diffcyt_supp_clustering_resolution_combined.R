##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: pAUC values
# - method: diffcyt methods (combined)
# 
# - supplementary results: varying clustering resolution
# 
# Lukas Weber, June 2018
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2
library(viridis)


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/supp_clustering_resolution"

load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_edgeR_supp_clustering_resolution.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_voom_supp_clustering_resolution.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_GLMM_supp_clustering_resolution.RData"))


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
plots_pAUC <- vector("list", length(thresholds) * length(cond_names))


for (th in 1:length(thresholds)) {
  
  for (j in 1:length(cond_names)) {
    
    
    # index to store plots in list
    ix <- (j * length(thresholds)) - (length(thresholds) - th)
    
    
    # -------------------------------------
    # Pre-processing steps for iCOBRA plots
    # -------------------------------------
    
    # create 'COBRAData' objects (one object per method)
    data_edgeR <- list(k_9    = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_9"]][[th]][[j]], 
                       k_25   = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_25"]][[th]][[j]], 
                       k_49   = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_49"]][[th]][[j]], 
                       k_100  = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_100"]][[th]][[j]], 
                       k_196  = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_196"]][[th]][[j]], 
                       k_400  = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_400"]][[th]][[j]], 
                       k_900  = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_900"]][[th]][[j]], 
                       k_1600 = out_diffcyt_DA_edgeR_supp_clustering_resolution[["k_1600"]][[th]][[j]])
    
    data_voom <- list(k_9    = out_diffcyt_DA_voom_supp_clustering_resolution[["k_9"]][[th]][[j]], 
                      k_25   = out_diffcyt_DA_voom_supp_clustering_resolution[["k_25"]][[th]][[j]], 
                      k_49   = out_diffcyt_DA_voom_supp_clustering_resolution[["k_49"]][[th]][[j]], 
                      k_100  = out_diffcyt_DA_voom_supp_clustering_resolution[["k_100"]][[th]][[j]], 
                      k_196  = out_diffcyt_DA_voom_supp_clustering_resolution[["k_196"]][[th]][[j]], 
                      k_400  = out_diffcyt_DA_voom_supp_clustering_resolution[["k_400"]][[th]][[j]], 
                      k_900  = out_diffcyt_DA_voom_supp_clustering_resolution[["k_900"]][[th]][[j]], 
                      k_1600 = out_diffcyt_DA_voom_supp_clustering_resolution[["k_1600"]][[th]][[j]])
    
    data_GLMM <- list(k_9    = out_diffcyt_DA_GLMM_supp_clustering_resolution[["k_9"]][[th]][[j]], 
                      k_25   = out_diffcyt_DA_GLMM_supp_clustering_resolution[["k_25"]][[th]][[j]], 
                      k_49   = out_diffcyt_DA_GLMM_supp_clustering_resolution[["k_49"]][[th]][[j]], 
                      k_100  = out_diffcyt_DA_GLMM_supp_clustering_resolution[["k_100"]][[th]][[j]], 
                      k_196  = out_diffcyt_DA_GLMM_supp_clustering_resolution[["k_196"]][[th]][[j]], 
                      k_400  = out_diffcyt_DA_GLMM_supp_clustering_resolution[["k_400"]][[th]][[j]], 
                      k_900  = out_diffcyt_DA_GLMM_supp_clustering_resolution[["k_900"]][[th]][[j]], 
                      k_1600 = out_diffcyt_DA_GLMM_supp_clustering_resolution[["k_1600"]][[th]][[j]])
    
    # check
    stopifnot(all(sapply(data_edgeR, function(d) all(d$spikein == data_edgeR[[1]]$spikein))))
    stopifnot(all(sapply(data_voom, function(d) all(d$spikein == data_voom[[1]]$spikein))))
    stopifnot(all(sapply(data_GLMM, function(d) all(d$spikein == data_GLMM[[1]]$spikein))))
    
    # note: provide all available values
    # 'padj' is required for threshold points on TPR-FDR curves
    # depending on availability, plotting functions use 'score', then 'pval', then 'padj'
    cobradata_edgeR <- COBRAData(pval = data.frame(k_9    = data_edgeR[["k_9"]][, "p_val"], 
                                                   k_25   = data_edgeR[["k_25"]][, "p_val"], 
                                                   k_49   = data_edgeR[["k_49"]][, "p_val"], 
                                                   k_100  = data_edgeR[["k_100"]][, "p_val"], 
                                                   k_196  = data_edgeR[["k_196"]][, "p_val"], 
                                                   k_400  = data_edgeR[["k_400"]][, "p_val"], 
                                                   k_900  = data_edgeR[["k_900"]][, "p_val"], 
                                                   k_1600 = data_edgeR[["k_1600"]][, "p_val"]), 
                                 padj = data.frame(k_9    = data_edgeR[["k_9"]][, "p_adj"], 
                                                   k_25   = data_edgeR[["k_25"]][, "p_adj"], 
                                                   k_49   = data_edgeR[["k_49"]][, "p_adj"], 
                                                   k_100  = data_edgeR[["k_100"]][, "p_adj"], 
                                                   k_196  = data_edgeR[["k_196"]][, "p_adj"], 
                                                   k_400  = data_edgeR[["k_400"]][, "p_adj"], 
                                                   k_900  = data_edgeR[["k_900"]][, "p_adj"], 
                                                   k_1600 = data_edgeR[["k_1600"]][, "p_adj"]), 
                                 truth = data.frame(spikein = data_edgeR[["k_9"]][, "spikein"]))
    
    cobradata_voom <- COBRAData(pval = data.frame(k_9    = data_voom[["k_9"]][, "p_val"], 
                                                   k_25   = data_voom[["k_25"]][, "p_val"], 
                                                   k_49   = data_voom[["k_49"]][, "p_val"], 
                                                   k_100  = data_voom[["k_100"]][, "p_val"], 
                                                   k_196  = data_voom[["k_196"]][, "p_val"], 
                                                   k_400  = data_voom[["k_400"]][, "p_val"], 
                                                   k_900  = data_voom[["k_900"]][, "p_val"], 
                                                   k_1600 = data_voom[["k_1600"]][, "p_val"]), 
                                 padj = data.frame(k_9    = data_voom[["k_9"]][, "p_adj"], 
                                                   k_25   = data_voom[["k_25"]][, "p_adj"], 
                                                   k_49   = data_voom[["k_49"]][, "p_adj"], 
                                                   k_100  = data_voom[["k_100"]][, "p_adj"], 
                                                   k_196  = data_voom[["k_196"]][, "p_adj"], 
                                                   k_400  = data_voom[["k_400"]][, "p_adj"], 
                                                   k_900  = data_voom[["k_900"]][, "p_adj"], 
                                                   k_1600 = data_voom[["k_1600"]][, "p_adj"]), 
                                 truth = data.frame(spikein = data_voom[["k_9"]][, "spikein"]))
    
    cobradata_GLMM <- COBRAData(pval = data.frame(k_9    = data_GLMM[["k_9"]][, "p_val"], 
                                                   k_25   = data_GLMM[["k_25"]][, "p_val"], 
                                                   k_49   = data_GLMM[["k_49"]][, "p_val"], 
                                                   k_100  = data_GLMM[["k_100"]][, "p_val"], 
                                                   k_196  = data_GLMM[["k_196"]][, "p_val"], 
                                                   k_400  = data_GLMM[["k_400"]][, "p_val"], 
                                                   k_900  = data_GLMM[["k_900"]][, "p_val"], 
                                                   k_1600 = data_GLMM[["k_1600"]][, "p_val"]), 
                                 padj = data.frame(k_9    = data_GLMM[["k_9"]][, "p_adj"], 
                                                   k_25   = data_GLMM[["k_25"]][, "p_adj"], 
                                                   k_49   = data_GLMM[["k_49"]][, "p_adj"], 
                                                   k_100  = data_GLMM[["k_100"]][, "p_adj"], 
                                                   k_196  = data_GLMM[["k_196"]][, "p_adj"], 
                                                   k_400  = data_GLMM[["k_400"]][, "p_adj"], 
                                                   k_900  = data_GLMM[["k_900"]][, "p_adj"], 
                                                   k_1600 = data_GLMM[["k_1600"]][, "p_adj"]), 
                                 truth = data.frame(spikein = data_GLMM[["k_9"]][, "spikein"]))
    
    # calculate performance scores
    # (note: can ignore warning messages when 'padj' not available)
    cobraperf_edgeR <- calculate_performance(cobradata_edgeR, binary_truth = "spikein", aspects = "roc")
    cobraperf_voom <- calculate_performance(cobradata_voom, binary_truth = "spikein", aspects = "roc")
    cobraperf_GLMM <- calculate_performance(cobradata_GLMM, binary_truth = "spikein", aspects = "roc")
    
    # color scheme
    colors <- rev(viridis(8))
    
    colors <- colors[1:length(data_edgeR)]
    names(colors) <- names(data_edgeR)
    
    # prepare plotting object
    cobraplot_edgeR <- prepare_data_for_plot(cobraperf_edgeR, colorscheme = colors, conditionalfill = FALSE)
    cobraplot_voom <- prepare_data_for_plot(cobraperf_voom, colorscheme = colors, conditionalfill = FALSE)
    cobraplot_GLMM <- prepare_data_for_plot(cobraperf_GLMM, colorscheme = colors, conditionalfill = FALSE)
    
    # re-order legend
    cobraplot_edgeR <- reorder_levels(cobraplot_edgeR, levels = names(data_edgeR))
    cobraplot_voom <- reorder_levels(cobraplot_voom, levels = names(data_edgeR))
    cobraplot_GLMM <- reorder_levels(cobraplot_GLMM, levels = names(data_edgeR))
    
    
    
    # -----------
    # pAUC values
    # -----------
    
    # calculate pAUC values
    
    pAUC_edgeR <- pAUC_voom <- pAUC_GLMM <- rep(NA, length(data_edgeR))
    names(pAUC_edgeR) <- names(pAUC_voom) <- names(pAUC_GLMM) <- names(data_edgeR)
    
    pAUC_list <- list(pAUC_edgeR, pAUC_voom, pAUC_GLMM)
    cobraperf_list <- list(cobraperf_edgeR, cobraperf_voom, cobraperf_GLMM)
    
    thresh <- 0.2
    
    
    for (m in 1:3) {
      
      pAUC <- pAUC_list[[m]]
      
      for (i in 1:length(pAUC)) {
        roc_i <- roc(cobraperf_list[[m]])[roc(cobraperf_list[[m]])$method == names(pAUC)[i], ]
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
      
      names(pAUC) <- gsub("k_", "", names(pAUC))
      
      # store pAUC values for each method
      pAUC_list[[m]] <- pAUC
      
    }
    
    
    # create data frame for plotting
    
    p_data <- vector("list", length = 3)
    method_names <- c("diffcyt-DA-edgeR", "diffcyt-DA-voom", "diffcyt-DA-GLMM")
    
    for (m in 1:length(p_data)) {
      p_data_m <- data.frame(pAUC = pAUC_list[[m]])
      p_data_m$n_clusters <- factor(names(pAUC_list[[m]]), levels = names(pAUC_list[[m]]))
      p_data_m$method <- method_names[m]
      p_data[[m]] <- p_data_m
    }
    
    p_data <- do.call("rbind", p_data)
    rownames(p_data) <- NULL
    
    p_data$method <- factor(p_data$method, levels = method_names)
    
    
    # colors for methods
    colors <- c("darkblue", "dodgerblue", "darkslategray3")
    names(colors) <- method_names
    
    # create plot
    p_pAUC <- 
      ggplot(p_data, aes(x = n_clusters, y = pAUC, group = method, color = method)) + 
      geom_line(lwd = 0.7) + 
      geom_point(size = 2.25) + 
      scale_color_manual(values = colors) + 
      ylim(c(0, 1)) + 
      xlab("Number of clusters") + 
      ylab(paste0("pAUC (FPR < ", thresh, ")")) + 
      ggtitle(paste0("AML-sim: clustering resolution, ", cond_names[j], ", ", gsub("pc$", "\\%", thresholds[th])), subtitle = "pAUC") + 
      theme_bw()
    
    plots_pAUC[[ix]] <- p_pAUC
    
    # save individual panel plot
    fn <- file.path(DIR_PLOTS, "panels", 
                    paste0("results_AML_sim_diffcyt_supp_clustering_resolution_pAUC_", thresholds[th], "_", cond_names[j], ".pdf"))
    ggsave(fn, width = 5.25, height = 3.5)
    
  }
}




########################
# Save multi-panel plots
########################

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
title_multi <- ggdraw() + draw_label(gsub(",.*$", "", plots_pAUC[[1]]$labels$title), fontface = "bold")
grid_multi <- plot_grid(title_multi, plots_multi, ncol = 1, rel_heights = c(1, 20))

# add combined legend
legend_multi <- get_legend(plots_pAUC[[1]] + theme(legend.position = "right", 
                                                   legend.title = element_text(size = 12, face = "bold"), 
                                                   legend.text = element_text(size = 12)))
grid_multi <- plot_grid(grid_multi, legend_multi, nrow = 1, rel_widths = c(5, 1))

# save plots
fn_multi <- file.path(DIR_PLOTS, paste0("results_AML_sim_diffcyt_supp_clustering_resolution_pAUC.pdf"))
ggsave(fn_multi, grid_multi, width = 10, height = 6)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



