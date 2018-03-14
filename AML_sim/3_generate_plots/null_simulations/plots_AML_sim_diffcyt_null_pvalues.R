##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: histograms of p-value distributions for null simulations
# - method: diffcyt methods
# 
# - null simulations
# 
# Lukas Weber, March 2018
##########################################################################################


library(ggplot2)
library(cowplot)  # note: cowplot masks 'ggsave' from ggplot2


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/null_simulations"

load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_edgeR_null.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_limma_null.RData"))
load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_GLMM_null.RData"))

load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_edgeR_null.RData"))
load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_limma_null.RData"))
load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_GLMM_null.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/null_simulations"




################
# Generate plots
################

# histograms of p-value distribution (raw p-values): at cluster level

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# random seeds
seed_names <- c("seed 1", "seed 2", "seed 3")

# color scheme
colors <- c("darkblue", "dodgerblue", "darkslategray3")



# ----------------
# diffcyt-DA-edgeR
# ----------------

# store plots in list
plots_edgeR <- vector("list", length(thresholds) * length(seed_names))

for (th in 1:length(thresholds)) {
  
  for (s in 1:length(seed_names)) {
    
    # index to store plots in list
    ix <- (th * length(seed_names)) - (length(seed_names) - s)
    
    d_plot <- out_clusters_diffcyt_DA_edgeR_null[[th]][[s]]
    
    # replace any NAs with 1s
    d_plot[is.na(d_plot[, "PValue"]), "PValue"] <- 1
    
    d_plot$method <- as.factor("diffcyt-DA-edgeR")
    
    p <- 
      ggplot(d_plot, aes(x = PValue, fill = method)) + 
      geom_histogram(bins = 20, color = "black") + 
      scale_fill_manual(values = colors[1]) + 
      ggtitle("AML-sim, null simulations: diffcyt-DA-edgeR", 
              subtitle = paste0("p-value distribution, threshold ", gsub("pc", "\\%", thresholds[th]), ", random seed ", s)) + 
      xlim(c(0, 1)) + 
      ylim(0, 35) + 
      xlab("p-value") + 
      theme_bw()
    
    plots_edgeR[[ix]] <- p
    
  }
}



# ----------------
# diffcyt-DA-limma
# ----------------

# store plots in list
plots_limma <- vector("list", length(thresholds) * length(seed_names))

for (th in 1:length(thresholds)) {
  
  for (s in 1:length(seed_names)) {
    
    # index to store plots in list
    ix <- (th * length(seed_names)) - (length(seed_names) - s)
    
    d_plot <- out_clusters_diffcyt_DA_limma_null[[th]][[s]]
    
    # replace any NAs with 1s
    d_plot[is.na(d_plot[, "P.Value"]), "P.Value"] <- 1
    
    d_plot$method <- as.factor("diffcyt-DA-limma")
    
    p <- 
      ggplot(d_plot, aes(x = P.Value, fill = method)) + 
      geom_histogram(bins = 20, color = "black") + 
      scale_fill_manual(values = colors[2]) + 
      ggtitle("AML-sim, null simulations: diffcyt-DA-limma", 
              subtitle = paste0("p-value distribution, threshold ", gsub("pc", "\\%", thresholds[th]), ", random seed ", s)) + 
      xlim(c(0, 1)) + 
      ylim(0, 35) + 
      xlab("p-value") + 
      theme_bw()
    
    plots_limma[[ix]] <- p
    
  }
}



# ---------------
# diffcyt-DA-GLMM
# ---------------

# store plots in list
plots_GLMM <- vector("list", length(thresholds) * length(seed_names))

for (th in 1:length(thresholds)) {
  
  for (s in 1:length(seed_names)) {
    
    # index to store plots in list
    ix <- (th * length(seed_names)) - (length(seed_names) - s)
    
    d_plot <- out_clusters_diffcyt_DA_GLMM_null[[th]][[s]]
    
    # replace any NAs with 1s
    d_plot[is.na(d_plot[, "p_vals"]), "p_vals"] <- 1
    
    d_plot$method <- as.factor("diffcyt-DA-GLMM")
    
    p <- 
      ggplot(d_plot, aes(x = p_vals, fill = method)) + 
      geom_histogram(bins = 20, color = "black") + 
      scale_fill_manual(values = colors[3]) + 
      ggtitle("AML-sim, null simulations: diffcyt-DA-GLMM", 
              subtitle = paste0("p-value distribution, threshold ", gsub("pc", "\\%", thresholds[th]), ", random seed ", s)) + 
      xlim(c(0, 1)) + 
      ylim(0, 35) + 
      xlab("p-value") + 
      theme_bw()
    
    plots_GLMM[[ix]] <- p
    
  }
}




##################################
# Multi-panel plots: one seed only
##################################

# ----------------
# diffcyt-DA-edgeR
# ----------------

plots_list <- plots_edgeR[c(1, 4, 7)]

# modify plot elements
plots_list <- lapply(plots_list, function(p) {
  p +
    labs(title = gsub("^.* ", "", p$labels$title), subtitle = gsub("^.*distribution, ", "", p$labels$subtitle)) + 
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
fn <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_DA_edgeR_null_pvalues_one_seed.pdf")
ggsave(fn, width = 9, height = 2.5)



# ----------------
# diffcyt-DA-limma
# ----------------

plots_list <- plots_limma[c(1, 4, 7)]

# modify plot elements
plots_list <- lapply(plots_list, function(p) {
  p +
    labs(title = gsub("^.* ", "", p$labels$title), subtitle = gsub("^.*distribution, ", "", p$labels$subtitle)) + 
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
fn <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_DA_limma_null_pvalues_one_seed.pdf")
ggsave(fn, width = 9, height = 2.5)



# ----------------
# diffcyt-DA-GLMM
# ----------------

plots_list <- plots_GLMM[c(1, 4, 7)]

# modify plot elements
plots_list <- lapply(plots_list, function(p) {
  p +
    labs(title = gsub("^.* ", "", p$labels$title), subtitle = gsub("^.*distribution, ", "", p$labels$subtitle)) + 
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
fn <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_DA_GLMM_null_pvalues_one_seed.pdf")
ggsave(fn, width = 9, height = 2.5)




##############################
# Multi-panel plots: all seeds
##############################

# ----------------
# diffcyt-DA-edgeR
# ----------------

plots_list <- plots_edgeR

# modify plot elements
plots_list <- lapply(plots_list, function(p) {
  p +
    labs(title = gsub("^.* ", "", p$labels$title), subtitle = gsub("^.*distribution, ", "", p$labels$subtitle)) + 
    theme(legend.position = "none")
})

plots_multi <- plot_grid(plotlist = plots_list,
                         nrow = 3, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_single <- gsub(":.*$", ": p-value distributions", plots_edgeR[[1]]$labels$title)
plots_title <- ggdraw() + draw_label(title_single)
plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 15))

# add combined legend
legend_single <- get_legend(plots_edgeR[[1]] + theme(legend.position = "right"))
plots_multi <- plot_grid(plots_multi, legend_single, nrow = 1, rel_widths = c(5, 1))

# save multi-panel plot
fn <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_DA_edgeR_null_pvalues_all_seeds.pdf")
ggsave(fn, width = 9, height = 6)



# ----------------
# diffcyt-DA-limma
# ----------------

plots_list <- plots_limma

# modify plot elements
plots_list <- lapply(plots_list, function(p) {
  p +
    labs(title = gsub("^.* ", "", p$labels$title), subtitle = gsub("^.*distribution, ", "", p$labels$subtitle)) + 
    theme(legend.position = "none")
})

plots_multi <- plot_grid(plotlist = plots_list,
                         nrow = 3, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_single <- gsub(":.*$", ": p-value distributions", plots_limma[[1]]$labels$title)
plots_title <- ggdraw() + draw_label(title_single)
plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 15))

# add combined legend
legend_single <- get_legend(plots_limma[[1]] + theme(legend.position = "right"))
plots_multi <- plot_grid(plots_multi, legend_single, nrow = 1, rel_widths = c(5, 1))

# save multi-panel plot
fn <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_DA_limma_null_pvalues_all_seeds.pdf")
ggsave(fn, width = 9, height = 6)



# ----------------
# diffcyt-DA-GLMM
# ----------------

plots_list <- plots_GLMM

# modify plot elements
plots_list <- lapply(plots_list, function(p) {
  p +
    labs(title = gsub("^.* ", "", p$labels$title), subtitle = gsub("^.*distribution, ", "", p$labels$subtitle)) + 
    theme(legend.position = "none")
})

plots_multi <- plot_grid(plotlist = plots_list,
                         nrow = 3, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_single <- gsub(":.*$", ": p-value distributions", plots_GLMM[[1]]$labels$title)
plots_title <- ggdraw() + draw_label(title_single)
plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 15))

# add combined legend
legend_single <- get_legend(plots_GLMM[[1]] + theme(legend.position = "right"))
plots_multi <- plot_grid(plots_multi, legend_single, nrow = 1, rel_widths = c(5, 1))

# save multi-panel plot
fn <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_DA_GLMM_null_pvalues_all_seeds.pdf")
ggsave(fn, width = 9, height = 6)




############################################
# Additional multi-panel plot for main paper
############################################

# showing one panel for each method


# create plot with combined legend (use to get legend object only)

d_plot_edgeR <- out_clusters_diffcyt_DA_edgeR_null[[1]][[1]][, c(1, 5, 6)]
d_plot_limma <- out_clusters_diffcyt_DA_limma_null[[1]][[1]][, c(1, 5, 6)]
d_plot_GLMM <- out_clusters_diffcyt_DA_GLMM_null[[1]][[1]]

colnames(d_plot_edgeR) <- colnames(d_plot_limma) <- colnames(d_plot_GLMM)

d_plot_edgeR$method <- "diffcyt-DA-edgeR"
d_plot_limma$method <- "diffcyt-DA-limma"
d_plot_GLMM$method <- "diffcyt-DA-GLMM"

d_plot <- rbind(d_plot_edgeR, d_plot_limma, d_plot_GLMM)

d_plot$method <- factor(d_plot$method)

# replace any NAs with 1s
d_plot[is.na(d_plot$p_vals), "p_vals"] <- 1
d_plot[is.na(d_plot$p_adj), "p_adj"] <- 1

p_legend <- 
  ggplot(d_plot, aes(x = p_vals, fill = method)) + 
  geom_histogram(position = "dodge", bins = 20, color = "black") + 
  scale_fill_manual(values = colors) + 
  ggtitle("AML-sim, null simulations", 
          subtitle = paste0("p-value distribution, threshold 5\\%, random seed 1")) + 
  xlim(c(0, 1)) + 
  ylim(c(0, 35)) + 
  xlab("p-value") + 
  theme_bw()


# create multi-panel plot

plots_list <- list(plots_edgeR[[1]], plots_limma[[1]], plots_GLMM[[1]])

# modify plot elements
plots_list <- lapply(plots_list, function(p) {
  p +
    labs(title = gsub("^.* ", "", p$labels$title), subtitle = gsub("^.*, ", "", p$labels$subtitle)) + 
    theme(legend.position = "none")
})

plots_multi <- plot_grid(plotlist = plots_list,
                         nrow = 1, ncol = 3, align = "hv", axis = "bl")

# add combined title
title_single <- gsub(":.*$", ": p-values", plots_edgeR[[1]]$labels$title)
plots_title <- ggdraw() + draw_label(title_single)
plots_multi <- plot_grid(plots_title, plots_multi, ncol = 1, rel_heights = c(1, 5.5))

# add combined legend
legend_single <- get_legend(p_legend + theme(legend.position = "right"))
plots_multi <- plot_grid(plots_multi, legend_single, nrow = 1, rel_widths = c(4.75, 1))

# save multi-panel plot
fn <- file.path(DIR_PLOTS, "results_AML_sim_diffcyt_null_pvalues_3_panels.pdf")
ggsave(fn, width = 8.5, height = 2.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



