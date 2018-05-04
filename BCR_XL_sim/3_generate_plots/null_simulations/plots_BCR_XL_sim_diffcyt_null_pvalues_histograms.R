##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: histograms of p-value distributions for null simulations
# - method: diffcyt methods
# 
# - null simulations
# 
# Lukas Weber, May 2018
##########################################################################################


library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL_sim/null_simulations"

load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_limma_null.RData"))
load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_LMM_null.RData"))

load(file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DS_limma_null.RData"))
load(file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DS_LMM_null.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/null_simulations"




################
# Generate plots
################

# histograms of p-value distribution (raw p-values): at cluster level, marker pS6 only

seed_names <- c("seed 1", "seed 2", "seed 3")



# ----------------
# diffcyt-DS-limma
# ----------------

plots_limma <- vector("list", length(seed_names))
names(plots_limma) <- seed_names

for (s in 1:length(seed_names)) {
  
  d_plot <- out_clusters_diffcyt_DS_limma_null[[s]][out_clusters_diffcyt_DS_limma_null[[s]]$marker == "pS6", ]
  
  dim(d_plot)
  
  # replace any NAs with 1s
  d_plot[is.na(d_plot$P.Value), "P.Value"] <- 1
  
  d_plot$method <- as.factor("diffcyt-DS-limma")
  
  p <- 
    ggplot(d_plot, aes(x = P.Value, fill = method)) + 
    geom_histogram(bins = 20, color = "black") + 
    scale_fill_manual(values = "firebrick1") + 
    ggtitle("BCR-XL-sim, null simulations: diffcyt-DS-limma", subtitle = paste("p-value distribution, random seed", s)) + 
    xlim(c(0, 1)) + 
    ylim(c(0, 13)) + 
    xlab("p-value") + 
    theme_bw()
  
  plots_limma[[s]] <- p
  
  fn <- file.path(DIR_PLOTS, "panels", paste0("results_BCR_XL_sim_diffcyt_DS_limma_null_pvalues", "_seed", s, "_hist.pdf"))
  ggsave(fn, width = 5.5, height = 3.75)
  
}



# --------------
# diffcyt-DS-LMM
# --------------

plots_LMM <- vector("list", length(seed_names))
names(plots_LMM) <- seed_names

for (s in 1:length(seed_names)) {
  
  d_plot <- out_clusters_diffcyt_DS_LMM_null[[s]][out_clusters_diffcyt_DS_LMM_null[[s]]$marker == "pS6", ]
  
  dim(d_plot)
  
  # replace any NAs with 1s
  d_plot[is.na(d_plot$p_val), "p_val"] <- 1
  
  d_plot$method <- as.factor("diffcyt-DS-LMM")
  
  p <- 
    ggplot(d_plot, aes(x = p_val, fill = method)) + 
    geom_histogram(bins = 20, color = "black") + 
    scale_fill_manual(values = "darkviolet") + 
    ggtitle("BCR-XL-sim, null simulations: diffcyt-DS-LMM", subtitle = paste("p-value distribution, random seed", s)) + 
    xlim(c(0, 1)) + 
    ylim(c(0, 13)) + 
    xlab("p-value") + 
    theme_bw()
  
  plots_LMM[[s]] <- p
  
  fn <- file.path(DIR_PLOTS, "panels", paste0("results_BCR_XL_sim_diffcyt_DS_LMM_null_pvalues", "_seed", s, "_hist.pdf"))
  ggsave(fn, width = 5.5, height = 3.75)
  
}




###################
# Multi-panel plots
###################

# ----------------
# diffcyt-DS-limma
# ----------------

plots_list <- plots_limma

# modify plot elements
plots_list <- lapply(plots_list, function(p) {
  p +
    labs(title = gsub("^.* ", "", p$labels$title), subtitle = gsub("^.*, ", "", p$labels$subtitle)) + 
    theme(legend.position = "none")
})

plots_multi <- plot_grid(plotlist = plots_list,
                         nrow = 1, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_single <- gsub(":.*$", ": p-value distributions", plots_limma[[1]]$labels$title)
plots_title <- ggdraw() + draw_label(title_single)
plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 5.5))

# add combined legend
legend_single <- get_legend(plots_limma[[1]] + theme(legend.position = "right"))
plots_multi <- plot_grid(plots_multi, legend_single, nrow = 1, rel_widths = c(5, 1))

# save multi-panel plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_limma_null_pvalues_hist.pdf")
ggsave(fn, width = 9, height = 2.5)



# --------------
# diffcyt-DS-LMM
# --------------

plots_list <- plots_LMM

# modify plot elements
plots_list <- lapply(plots_list, function(p) {
  p +
    labs(title = gsub("^.* ", "", p$labels$title), subtitle = gsub("^.*, ", "", p$labels$subtitle)) + 
    theme(legend.position = "none")
})

plots_multi <- plot_grid(plotlist = plots_list,
                         nrow = 1, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_single <- gsub(":.*$", ": p-value distributions", plots_LMM[[1]]$labels$title)
plots_title <- ggdraw() + draw_label(title_single)
plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 5.5))

# add combined legend
legend_single <- get_legend(plots_LMM[[1]] + theme(legend.position = "right"))
plots_multi <- plot_grid(plots_multi, legend_single, nrow = 1, rel_widths = c(5, 1))

# save multi-panel plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_LMM_null_pvalues_hist.pdf")
ggsave(fn, width = 9, height = 2.5)




############################################
# Additional multi-panel plot for main paper
############################################

# showing one panel for each method; one random seed only


# create plot with combined legend (use to get legend object only)

d_plot_limma <- out_clusters_diffcyt_DS_limma_null[[1]][out_clusters_diffcyt_DS_limma_null[[1]]$marker == "pS6", ]
d_plot_LMM <- out_clusters_diffcyt_DS_LMM_null[[1]][out_clusters_diffcyt_DS_LMM_null[[1]]$marker == "pS6", ]

d_plot_limma <- d_plot_limma[, c(1, 2, 7, 8)]
colnames(d_plot_limma) <- colnames(d_plot_LMM)

d_plot_limma$method <- "diffcyt-DS-limma"
d_plot_LMM$method <- "diffcyt-DS-LMM"

d_plot <- rbind(d_plot_limma, d_plot_LMM)

d_plot$method <- factor(d_plot$method)

# replace any NAs with 1s
d_plot[is.na(d_plot$p_val), "p_val"] <- 1

p_legend <- 
  ggplot(d_plot, aes(x = p_val, fill = method)) + 
  geom_histogram(position = "dodge", bins = 20, color = "black") + 
  scale_fill_manual(values = c("firebrick1", "darkviolet")) + 
  ggtitle("BCR-XL-sim, null simulations", subtitle = paste("p-value distribution, random seed 1")) + 
  xlim(c(0, 1)) + 
  ylim(c(0, 13)) + 
  xlab("p-value") + 
  theme_bw()


# create multi-panel plot

plots_list <- list(plots_limma[[1]], plots_LMM[[1]])

# modify plot elements
plots_list <- lapply(plots_list, function(p) {
  p +
    labs(title = gsub("^.* ", "", p$labels$title), subtitle = gsub("^.*, ", "", p$labels$subtitle)) + 
    theme(legend.position = "none")
})

plots_multi <- plot_grid(plotlist = plots_list,
                         nrow = 1, ncol = 2, align = "hv", axis = "bl")

# add combined title
title_single <- gsub(":.*$", ": p-values", plots_limma[[1]]$labels$title)
plots_title <- ggdraw() + draw_label(title_single)
plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 5.5))

# add combined legend
legend_single <- get_legend(p_legend + theme(legend.position = "right"))
plots_multi <- plot_grid(plots_multi, legend_single, nrow = 1, rel_widths = c(3.25, 1))

# save multi-panel plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_null_pvalues_2_panels_hist.pdf")
ggsave(fn, width = 6, height = 2.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



