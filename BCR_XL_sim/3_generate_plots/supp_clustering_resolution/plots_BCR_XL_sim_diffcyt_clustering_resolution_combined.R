##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: pAUC values
# - method: diffcyt methods (combined)
# 
# - supplementary results: varying clustering resolution
# 
# Lukas Weber, May 2019
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2
library(viridis)


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL_sim/supp_clustering_resolution"

load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_limma_supp_clustering_resolution.RData"))
load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_LMM_supp_clustering_resolution.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/supp_clustering_resolution"




################
# Generate plots
################

# -------------------------------------
# Pre-processing steps for iCOBRA plots
# -------------------------------------

# create 'COBRAData' objects (one object per method)
data_limma <- list(k_9    = out_diffcyt_DS_limma_supp_clustering_resolution[["k_9"]], 
                   k_25   = out_diffcyt_DS_limma_supp_clustering_resolution[["k_25"]], 
                   k_49   = out_diffcyt_DS_limma_supp_clustering_resolution[["k_49"]], 
                   k_100  = out_diffcyt_DS_limma_supp_clustering_resolution[["k_100"]], 
                   k_196  = out_diffcyt_DS_limma_supp_clustering_resolution[["k_196"]], 
                   k_400  = out_diffcyt_DS_limma_supp_clustering_resolution[["k_400"]], 
                   k_900  = out_diffcyt_DS_limma_supp_clustering_resolution[["k_900"]], 
                   k_1600 = out_diffcyt_DS_limma_supp_clustering_resolution[["k_1600"]])

data_LMM <- list(k_9    = out_diffcyt_DS_LMM_supp_clustering_resolution[["k_9"]], 
                 k_25   = out_diffcyt_DS_LMM_supp_clustering_resolution[["k_25"]], 
                 k_49   = out_diffcyt_DS_LMM_supp_clustering_resolution[["k_49"]], 
                 k_100  = out_diffcyt_DS_LMM_supp_clustering_resolution[["k_100"]], 
                 k_196  = out_diffcyt_DS_LMM_supp_clustering_resolution[["k_196"]], 
                 k_400  = out_diffcyt_DS_LMM_supp_clustering_resolution[["k_400"]], 
                 k_900  = out_diffcyt_DS_LMM_supp_clustering_resolution[["k_900"]], 
                 k_1600 = out_diffcyt_DS_LMM_supp_clustering_resolution[["k_1600"]])

# check
stopifnot(all(sapply(data_limma, function(d) all(d$B_cell == data_limma[[1]]$B_cell))))
stopifnot(all(sapply(data_LMM, function(d) all(d$B_cell == data_LMM[[1]]$B_cell))))

# note: provide all available values
# 'padj' is required for threshold points on TPR-FDR curves
# depending on availability, plotting functions use 'score', then 'pval', then 'padj'
cobradata_limma <- COBRAData(pval = data.frame(k_9    = data_limma[["k_9"]][, "p_val"], 
                                               k_25   = data_limma[["k_25"]][, "p_val"], 
                                               k_49   = data_limma[["k_49"]][, "p_val"], 
                                               k_100  = data_limma[["k_100"]][, "p_val"], 
                                               k_196  = data_limma[["k_196"]][, "p_val"], 
                                               k_400  = data_limma[["k_400"]][, "p_val"], 
                                               k_900  = data_limma[["k_900"]][, "p_val"], 
                                               k_1600 = data_limma[["k_1600"]][, "p_val"]), 
                             padj = data.frame(k_9    = data_limma[["k_9"]][, "p_adj"], 
                                               k_25   = data_limma[["k_25"]][, "p_adj"], 
                                               k_49   = data_limma[["k_49"]][, "p_adj"], 
                                               k_100  = data_limma[["k_100"]][, "p_adj"], 
                                               k_196  = data_limma[["k_196"]][, "p_adj"], 
                                               k_400  = data_limma[["k_400"]][, "p_adj"], 
                                               k_900  = data_limma[["k_900"]][, "p_adj"], 
                                               k_1600 = data_limma[["k_1600"]][, "p_adj"]), 
                             truth = data.frame(B_cell = data_limma[["k_9"]][, "B_cell"]))

cobradata_LMM <- COBRAData(pval = data.frame(k_9    = data_LMM[["k_9"]][, "p_val"], 
                                             k_25   = data_LMM[["k_25"]][, "p_val"], 
                                             k_49   = data_LMM[["k_49"]][, "p_val"], 
                                             k_100  = data_LMM[["k_100"]][, "p_val"], 
                                             k_196  = data_LMM[["k_196"]][, "p_val"], 
                                             k_400  = data_LMM[["k_400"]][, "p_val"], 
                                             k_900  = data_LMM[["k_900"]][, "p_val"], 
                                             k_1600 = data_LMM[["k_1600"]][, "p_val"]), 
                           padj = data.frame(k_9    = data_LMM[["k_9"]][, "p_adj"], 
                                             k_25   = data_LMM[["k_25"]][, "p_adj"], 
                                             k_49   = data_LMM[["k_49"]][, "p_adj"], 
                                             k_100  = data_LMM[["k_100"]][, "p_adj"], 
                                             k_196  = data_LMM[["k_196"]][, "p_adj"], 
                                             k_400  = data_LMM[["k_400"]][, "p_adj"], 
                                             k_900  = data_LMM[["k_900"]][, "p_adj"], 
                                             k_1600 = data_LMM[["k_1600"]][, "p_adj"]), 
                           truth = data.frame(B_cell = data_LMM[["k_9"]][, "B_cell"]))

# calculate performance scores
# (note: can ignore warning messages when 'padj' not available)
cobraperf_limma <- calculate_performance(cobradata_limma, binary_truth = "B_cell", aspects = "roc")
cobraperf_LMM <- calculate_performance(cobradata_LMM, binary_truth = "B_cell", aspects = "roc")

# color scheme
colors <- rev(viridis(8))

colors <- colors[1:length(data_limma)]
names(colors) <- names(data_limma)

# prepare plotting object
cobraplot_limma <- prepare_data_for_plot(cobraperf_limma, colorscheme = colors, conditionalfill = FALSE)
cobraplot_LMM <- prepare_data_for_plot(cobraperf_LMM, colorscheme = colors, conditionalfill = FALSE)

# re-order legend
cobraplot_limma <- reorder_levels(cobraplot_limma, levels = names(data_limma))
cobraplot_LMM <- reorder_levels(cobraplot_LMM, levels = names(data_LMM))



# -----------
# pAUC values
# -----------

# calculate pAUC values

pAUC_limma <- pAUC_LMM <- rep(NA, length(data_limma))
names(pAUC_limma) <- names(pAUC_LMM) <- names(data_limma)

pAUC_list <- list(pAUC_limma, pAUC_LMM)
cobraperf_list <- list(cobraperf_limma, cobraperf_LMM)

thresh <- 0.2


for (m in 1:2) {
  
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

p_data <- vector("list", length = 2)
method_names <- c("diffcyt-DS-limma", "diffcyt-DS-LMM")

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
colors <- c("firebrick1", "darkviolet")
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
  ggtitle("BCR-XL-sim: clustering resolution", subtitle = "pAUC") + 
  theme_bw()

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_supp_clustering_resolution_pAUC.pdf")
ggsave(fn, width = 5.25, height = 3.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



