##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: density plots of p-value distributions for null simulations
# - method: diffcyt methods
# 
# - null simulations
# 
# Lukas Weber, April 2018
##########################################################################################


library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/null_simulations"

load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_edgeR_null.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_voom_null.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_GLMM_null.RData"))

load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_edgeR_null.RData"))
load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_voom_null.RData"))
load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_GLMM_null.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/null_simulations"




################
# Generate plots
################

# density plots of p-value distribution (raw p-values): at cluster level

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# color scheme
colors <- c("darkblue", "dodgerblue", "darkslategray3")



# ----------------
# diffcyt-DA-edgeR
# ----------------

# store plots in list
plots_edgeR <- vector("list", length(thresholds))

for (th in 1:length(thresholds)) {
  
  d_plot <- rbind(
    cbind(out_clusters_diffcyt_DA_edgeR_null[[th]][[1]], seed = "seed 1"), 
    cbind(out_clusters_diffcyt_DA_edgeR_null[[th]][[2]], seed = "seed 2"), 
    cbind(out_clusters_diffcyt_DA_edgeR_null[[th]][[3]], seed = "seed 3")
  )
  
  # replace any NAs with 1s
  d_plot[is.na(d_plot[, "PValue"]), "PValue"] <- 1
  
  d_plot$method <- as.factor("diffcyt-DA-edgeR")
  
  p <- 
    ggplot(d_plot, aes(x = PValue, linetype = seed, fill = method)) + 
    geom_density(adjust = 0.75, alpha = 0.5) + 
    ggtitle("AML-sim, null simulations: diffcyt-DA-edgeR", 
            subtitle = paste0("p-value distributions, threshold ", gsub("pc", "\\%", thresholds[th]))) + 
    scale_linetype_discrete(name = "random seed") + 
    scale_fill_manual(values = colors[1]) + 
    xlim(c(0, 1)) + 
    ylim(0, 2) + 
    xlab("p-value") + 
    theme_bw() + 
    guides(fill = guide_legend(order = 1), 
           linetype = guide_legend(order = 2))
  
  plots_edgeR[[th]] <- p
  
  fn <- file.path(DIR_PLOTS, "panels", paste0("results_BCR_XL_sim_diffcyt_DA_edgeR_null_pvalues_densities_", thresholds[th], ".pdf"))
  ggsave(fn, width = 5.5, height = 3.75)
  
}



# ---------------
# diffcyt-DA-voom
# ---------------

# store plots in list
plots_voom <- vector("list", length(thresholds))

for (th in 1:length(thresholds)) {
  
  d_plot <- rbind(
    cbind(out_clusters_diffcyt_DA_voom_null[[th]][[1]], seed = "seed 1"), 
    cbind(out_clusters_diffcyt_DA_voom_null[[th]][[2]], seed = "seed 2"), 
    cbind(out_clusters_diffcyt_DA_voom_null[[th]][[3]], seed = "seed 3")
  )
  
  # replace any NAs with 1s
  d_plot[is.na(d_plot[, "P.Value"]), "P.Value"] <- 1
  
  d_plot$method <- as.factor("diffcyt-DA-voom")
  
  p <- 
    ggplot(d_plot, aes(x = P.Value, linetype = seed, fill = method)) + 
    geom_density(adjust = 0.75, alpha = 0.5) + 
    ggtitle("AML-sim, null simulations: diffcyt-DA-voom", 
            subtitle = paste0("p-value distributions, threshold ", gsub("pc", "\\%", thresholds[th]))) + 
    scale_linetype_discrete(name = "random seed") + 
    scale_fill_manual(values = colors[2]) + 
    xlim(c(0, 1)) + 
    ylim(0, 2) + 
    xlab("p-value") + 
    theme_bw() + 
    guides(fill = guide_legend(order = 1), 
           linetype = guide_legend(order = 2))
  
  plots_voom[[th]] <- p
  
  fn <- file.path(DIR_PLOTS, "panels", paste0("results_BCR_XL_sim_diffcyt_DA_voom_null_pvalues_densities_", thresholds[th], ".pdf"))
  ggsave(fn, width = 5.5, height = 3.75)
  
}



# ---------------
# diffcyt-DA-GLMM
# ---------------

# store plots in list
plots_GLMM <- vector("list", length(thresholds))

for (th in 1:length(thresholds)) {
  
  d_plot <- rbind(
    cbind(out_clusters_diffcyt_DA_GLMM_null[[th]][[1]], seed = "seed 1"), 
    cbind(out_clusters_diffcyt_DA_GLMM_null[[th]][[2]], seed = "seed 2"), 
    cbind(out_clusters_diffcyt_DA_GLMM_null[[th]][[3]], seed = "seed 3")
  )
  
  # replace any NAs with 1s
  d_plot[is.na(d_plot[, "p_vals"]), "p_vals"] <- 1
  
  d_plot$method <- as.factor("diffcyt-DA-GLMM")
  
  p <- 
    ggplot(d_plot, aes(x = p_vals, linetype = seed, fill = method)) + 
    geom_density(adjust = 0.75, alpha = 0.5) + 
    ggtitle("AML-sim, null simulations: diffcyt-DA-GLMM", 
            subtitle = paste0("p-value distributions, threshold ", gsub("pc", "\\%", thresholds[th]))) + 
    scale_linetype_discrete(name = "random seed") + 
    scale_fill_manual(values = colors[3]) + 
    xlim(c(0, 1)) + 
    ylim(0, 2) + 
    xlab("p-value") + 
    theme_bw() + 
    guides(fill = guide_legend(order = 1), 
           linetype = guide_legend(order = 2))
  
  plots_GLMM[[th]] <- p
  
  fn <- file.path(DIR_PLOTS, "panels", paste0("results_BCR_XL_sim_diffcyt_DA_GLMM_null_pvalues_densities_", thresholds[th], ".pdf"))
  ggsave(fn, width = 5.5, height = 3.75)
  
}




###################
# Multi-panel plots
###################

# ----------------
# diffcyt-DA-edgeR
# ----------------

plots_list <- plots_edgeR

# modify plot elements
plots_list <- lapply(plots_list, function(p) {
  p +
    labs(title = gsub("^.*distributions, ", "", p$labels$subtitle), subtitle = element_blank()) + 
    theme(legend.position = "none")
})

plots_multi <- plot_grid(plotlist = plots_list,
                         nrow = 1, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_single <- gsub(":.*$", ": p-value distributions", plots_edgeR[[1]]$labels$title)
plots_title <- ggdraw() + draw_label(title_single)
plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 5.5))

# add combined legend
legend_single <- get_legend(plots_edgeR[[1]] + theme(legend.position = "right"))
plots_multi <- plot_grid(plots_multi, legend_single, nrow = 1, rel_widths = c(5, 1))

# save multi-panel plot
fn <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_DA_edgeR_null_pvalues_densities.pdf")
ggsave(fn, width = 9, height = 2.5)



# ---------------
# diffcyt-DA-voom
# ---------------

plots_list <- plots_voom

# modify plot elements
plots_list <- lapply(plots_list, function(p) {
  p +
    labs(title = gsub("^.*distributions, ", "", p$labels$subtitle), subtitle = element_blank()) + 
    theme(legend.position = "none")
})

plots_multi <- plot_grid(plotlist = plots_list,
                         nrow = 1, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_single <- gsub(":.*$", ": p-value distributions", plots_voom[[1]]$labels$title)
plots_title <- ggdraw() + draw_label(title_single)
plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 5.5))

# add combined legend
legend_single <- get_legend(plots_voom[[1]] + theme(legend.position = "right"))
plots_multi <- plot_grid(plots_multi, legend_single, nrow = 1, rel_widths = c(5, 1))

# save multi-panel plot
fn <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_DA_voom_null_pvalues_densities.pdf")
ggsave(fn, width = 9, height = 2.5)



# ---------------
# diffcyt-DA-GLMM
# ---------------

plots_list <- plots_GLMM

# modify plot elements
plots_list <- lapply(plots_list, function(p) {
  p +
    labs(title = gsub("^.*distributions, ", "", p$labels$subtitle), subtitle = element_blank()) + 
    theme(legend.position = "none")
})

plots_multi <- plot_grid(plotlist = plots_list,
                         nrow = 1, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_single <- gsub(":.*$", ": p-value distributions", plots_GLMM[[1]]$labels$title)
plots_title <- ggdraw() + draw_label(title_single)
plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 5.5))

# add combined legend
legend_single <- get_legend(plots_GLMM[[1]] + theme(legend.position = "right"))
plots_multi <- plot_grid(plots_multi, legend_single, nrow = 1, rel_widths = c(5, 1))

# save multi-panel plot
fn <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_DA_GLMM_null_pvalues_densities.pdf")
ggsave(fn, width = 9, height = 2.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



