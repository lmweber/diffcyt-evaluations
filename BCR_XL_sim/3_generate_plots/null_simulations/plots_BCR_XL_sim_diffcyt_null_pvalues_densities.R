##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: density plots of p-value distributions for null simulations
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

# density plots of p-value distribution (raw p-values): at cluster level, marker pS6 only


# ----------------
# diffcyt-DS-limma
# ----------------

d_plot <- rbind(
  cbind(out_clusters_diffcyt_DS_limma_null[[1]][out_clusters_diffcyt_DS_limma_null[[1]]$marker == "pS6", ], seed = "seed 1"), 
  cbind(out_clusters_diffcyt_DS_limma_null[[2]][out_clusters_diffcyt_DS_limma_null[[2]]$marker == "pS6", ], seed = "seed 2"), 
  cbind(out_clusters_diffcyt_DS_limma_null[[3]][out_clusters_diffcyt_DS_limma_null[[3]]$marker == "pS6", ], seed = "seed 3")
)

# replace any NAs with 1s
d_plot[is.na(d_plot$p_val), "p_val"] <- 1

d_plot$method <- as.factor("diffcyt-DS-limma")

p_limma <- 
  ggplot(d_plot, aes(x = p_val, linetype = seed, fill = method)) + 
  geom_density(adjust = 0.75, alpha = 0.5) + 
  ggtitle("BCR-XL-sim, null simulations: diffcyt-DS-limma", subtitle = paste("p-value distributions")) + 
  scale_linetype_discrete(name = "random seed") + 
  scale_fill_manual(values = "firebrick1") + 
  xlim(c(0, 1)) + 
  ylim(c(0, 2)) + 
  xlab("p-value") + 
  theme_bw() + 
  guides(fill = guide_legend(order = 1), 
         linetype = guide_legend(order = 2))

fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_limma_null_pvalues_densities.pdf")
ggsave(fn, width = 5.5, height = 3.75)



# --------------
# diffcyt-DS-LMM
# --------------

d_plot <- rbind(
  cbind(out_clusters_diffcyt_DS_LMM_null[[1]][out_clusters_diffcyt_DS_LMM_null[[1]]$marker == "pS6", ], seed = "seed 1"), 
  cbind(out_clusters_diffcyt_DS_LMM_null[[2]][out_clusters_diffcyt_DS_LMM_null[[2]]$marker == "pS6", ], seed = "seed 2"), 
  cbind(out_clusters_diffcyt_DS_LMM_null[[3]][out_clusters_diffcyt_DS_LMM_null[[3]]$marker == "pS6", ], seed = "seed 3")
)

# replace any NAs with 1s
d_plot[is.na(d_plot$p_val), "p_val"] <- 1

d_plot$method <- as.factor("diffcyt-DS-LMM")

p_LMM <- 
  ggplot(d_plot, aes(x = p_val, linetype = seed, fill = method)) + 
  geom_density(adjust = 0.75, alpha = 0.5) + 
  ggtitle("BCR-XL-sim, null simulations: diffcyt-DS-LMM", subtitle = paste("p-value distributions")) + 
  scale_linetype_discrete(name = "random seed") + 
  scale_fill_manual(values = "darkviolet") + 
  xlim(c(0, 1)) + 
  ylim(c(0, 2)) + 
  xlab("p-value") + 
  theme_bw() + 
  guides(fill = guide_legend(order = 1), 
         linetype = guide_legend(order = 2))

fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_LMM_null_pvalues_densities.pdf")
ggsave(fn, width = 5.5, height = 3.75)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



