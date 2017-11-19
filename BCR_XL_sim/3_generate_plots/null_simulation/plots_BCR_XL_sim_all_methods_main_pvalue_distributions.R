##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: histograms of p-value distributions for null simulation
# - method: diffcyt methods
# 
# - main results
# 
# Lukas Weber, November 2017
##########################################################################################


# note: showing 'diffcyt' methods only


library(ggplot2)


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL_sim/null_simulation"

load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_med_null.RData"))
load(file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DS_med_null.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/null_simulation"




################
# Generate plots
################

# --------------------------------------------------------------
# Histogram of p-value distribution: cell level, marker pS6 only
# --------------------------------------------------------------

d_plot <- out_diffcyt_DS_med_null

dim(d_plot)

ggplot(d_plot, aes(x = p_vals)) + 
  geom_histogram(bins = 25, color = "black", fill = "darkblue") + 
  ggtitle("BCR-XL-sim, null simulation, p-value distribution (cell level): diffcyt-DS-med") + 
  xlab("p-values") + 
  theme_bw()

fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_med_null_pvalues_cell.pdf")
ggsave(fn, width = 7.5, height = 6.5)



# -----------------------------------------------------------------
# Histogram of p-value distribution: cluster level, marker pS6 only
# -----------------------------------------------------------------

d_plot <- out_clusters_diffcyt_DS_med_null[out_clusters_diffcyt_DS_med_null$marker == "pS6(Yb172)Dd", ]

dim(d_plot)

# replace any NAs with 1s (same as at cell level)
d_plot[is.na(d_plot$P.Value), "P.Value"] <- 1

ggplot(d_plot, aes(x = P.Value)) + 
  geom_histogram(bins = 25, color = "black", fill = "darkblue") + 
  ggtitle("BCR-XL-sim, null simulation, p-value distribution (cluster level): diffcyt-DS-med") + 
  scale_y_continuous(breaks = seq(0, 16, by = 2)) + 
  xlab("p-values") + 
  theme_bw()

fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_DS_med_null_pvalues_cluster.pdf")
ggsave(fn, width = 7.5, height = 6.5)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



