##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: histograms of p-value distributions for null simulations
# - method: diffcyt methods
# 
# - null simulations
# 
# Lukas Weber, January 2018
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


# diffcyt-DS-limma

plots_limma <- vector("list", length(seed_names))
names(plots_limma) <- seed_names

for (s in 1:length(seed_names)) {
  
  d_plot <- out_clusters_diffcyt_DS_limma_null[[s]][out_clusters_diffcyt_DS_limma_null[[s]]$marker == "pS6", ]
  
  dim(d_plot)
  
  # replace any NAs with 1s
  d_plot[is.na(d_plot$P.Value), "P.Value"] <- 1
  
  p <- 
    ggplot(d_plot, aes(x = P.Value)) + 
    geom_histogram(bins = 20, color = "black", fill = "darkturquoise") + 
    ggtitle("BCR-XL-sim, null simulations: diffcyt-DS-limma", subtitle = paste("p-value distribution, random seed", s)) + 
    scale_y_continuous(breaks = seq(0, 16, by = 2)) + 
    xlim(c(0, 1)) + 
    xlab("p-value") + 
    theme_bw()
  
  plots_limma[[s]] <- p
  
  fn <- file.path(DIR_PLOTS, "panels", paste0("results_BCR_XL_sim_diffcyt_DS_limma_null_pvalues", "_seed", s, ".pdf"))
  ggsave(fn, width = 4.5, height = 3.75)
  
}


# diffcyt-DS-LMM

plots_LMM <- vector("list", length(seed_names))
names(plots_LMM) <- seed_names

for (s in 1:length(seed_names)) {
  
  d_plot <- out_clusters_diffcyt_DS_LMM_null[[s]][out_clusters_diffcyt_DS_LMM_null[[s]]$marker == "pS6", ]
  
  dim(d_plot)
  
  # replace any NAs with 1s
  d_plot[is.na(d_plot$p_vals), "p_vals"] <- 1
  
  p <- 
    ggplot(d_plot, aes(x = p_vals)) + 
    geom_histogram(bins = 20, color = "black", fill = "darkslategray4") + 
    ggtitle("BCR-XL-sim, null simulations: diffcyt-DS-LMM", subtitle = paste("p-value distribution, random seed", s)) + 
    scale_y_continuous(breaks = seq(0, 16, by = 2)) + 
    xlim(c(0, 1)) + 
    xlab("p-value") + 
    theme_bw()
  
  plots_LMM[[s]] <- p
  
  fn <- file.path(DIR_PLOTS, "panels", paste0("results_BCR_XL_sim_diffcyt_DS_LMM_null_pvalues", "_seed", s, ".pdf"))
  ggsave(fn, width = 4.5, height = 3.75)
  
}




##################
# Multi-panel plot
##################

plots_list <- c(plots_limma, plots_LMM)

# modify plot elements
plots_list <- lapply(plots_list, function(p) {
  p +
    labs(title = gsub("^.* ", "", p$labels$title), subtitle = gsub("^.*, ", "", p$labels$subtitle))
})

plots_multi <- plot_grid(plotlist = plots_list,
                         nrow = 2, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_single <- gsub(":.*$", ": p-value distributions", plots_limma[[1]]$labels$title)
plots_title <- ggdraw() + draw_label(title_single)
plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 12))

# save multi-panel plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_null_pvalues.pdf")
ggsave(fn, width = 8, height = 5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



