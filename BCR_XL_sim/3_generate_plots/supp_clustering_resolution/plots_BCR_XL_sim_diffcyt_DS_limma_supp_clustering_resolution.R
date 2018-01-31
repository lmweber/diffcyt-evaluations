##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: performance metrics
# - method: diffcyt-DS-limma
# 
# - supplementary results: varying clustering resolution
# 
# Lukas Weber, January 2018
##########################################################################################


library(iCOBRA)
library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2
library(viridis)


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL_sim/supp_clustering_resolution"

load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_limma_supp_clustering_resolution.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/supp_clustering_resolution"




################
# Generate plots
################

# -------------------------------------
# Pre-processing steps for iCOBRA plots
# -------------------------------------

# create 'COBRAData' object
data <- list(k_9   = out_diffcyt_DS_limma_supp_clustering_resolution[["k_9"]], 
             k_25  = out_diffcyt_DS_limma_supp_clustering_resolution[["k_25"]], 
             k_49  = out_diffcyt_DS_limma_supp_clustering_resolution[["k_49"]], 
             k_100 = out_diffcyt_DS_limma_supp_clustering_resolution[["k_100"]], 
             k_196 = out_diffcyt_DS_limma_supp_clustering_resolution[["k_196"]], 
             k_400 = out_diffcyt_DS_limma_supp_clustering_resolution[["k_400"]], 
             k_900 = out_diffcyt_DS_limma_supp_clustering_resolution[["k_900"]])

# check
stopifnot(all(sapply(data, function(d) all(d$B_cell == data[[1]]$B_cell))))

# note: provide all available values
# 'padj' is required for threshold points on TPR-FDR curves
# depending on availability, plotting functions use 'score', then 'pval', then 'padj'
cobradata <- COBRAData(pval = data.frame(k_9   = data[["k_9"]][, "p_vals"], 
                                         k_25  = data[["k_25"]][, "p_vals"], 
                                         k_49  = data[["k_49"]][, "p_vals"], 
                                         k_100 = data[["k_100"]][, "p_vals"], 
                                         k_196 = data[["k_196"]][, "p_vals"], 
                                         k_400 = data[["k_400"]][, "p_vals"], 
                                         k_900 = data[["k_900"]][, "p_vals"]), 
                       padj = data.frame(k_9   = data[["k_9"]][, "p_adj"], 
                                         k_25  = data[["k_25"]][, "p_adj"], 
                                         k_49  = data[["k_49"]][, "p_adj"], 
                                         k_100 = data[["k_100"]][, "p_adj"], 
                                         k_196 = data[["k_196"]][, "p_adj"], 
                                         k_400 = data[["k_400"]][, "p_adj"], 
                                         k_900 = data[["k_900"]][, "p_adj"]), 
                       truth = data.frame(B_cell = data[["k_9"]][, "B_cell"]))

# calculate performance scores
# (note: can ignore warning messages when 'padj' not available)
cobraperf <- calculate_performance(cobradata, 
                                   binary_truth = "B_cell", 
                                   aspects = "roc")

# color scheme
colors <- rev(viridis(7))

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
  coord_fixed() + 
  xlab("False positive rate") + 
  ylab("True positive rate") + 
  ggtitle("BCR-XL-sim: clustering resolution", subtitle = "ROC curves, diffcyt-DS-limma") + 
  theme_bw() + 
  theme(strip.text.x = element_blank()) + 
  guides(color = guide_legend("no. clusters"))

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_limma_supp_clustering_resolution_ROC.pdf")
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
p_data <- as.data.frame(pAUC)
p_data$n_clusters <- factor(names(pAUC), levels = names(pAUC))
p_data$group <- "pAUC"


# create plot
p_pAUC <- 
  ggplot(p_data, aes(x = n_clusters, y = pAUC, group = group, color = n_clusters)) + 
  geom_line(lwd = 0.7) + 
  geom_point(size = 2.25) + 
  scale_color_manual(values = colors) + 
  ylim(c(0, 1)) + 
  xlab("Number of clusters") + 
  ylab(paste0("pAUC (FPR < ", thresh, ")")) + 
  ggtitle("BCR-XL-sim: clustering resolution", subtitle = "pAUC, diffcyt-DS-limma") + 
  theme_bw() + 
  guides(color = guide_legend("no. clusters"))

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_limma_supp_clustering_resolution_pAUC.pdf")
ggsave(fn, width = 4.75, height = 3.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



