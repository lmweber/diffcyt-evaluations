##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: clustering performance
# - method: diffcyt methods
# 
# - main results
# 
# Lukas Weber, February 2018
##########################################################################################


# note: clustering step is the same for all 'diffcyt' methods (diffcyt-DS-limma, diffcyt-DS-LMM)


library(SummarizedExperiment)
library(reshape2)
library(ggplot2)


# load saved results
DIR_RDATA <- "../../../../RData/BCR_XL_sim/main"

# note: only need to load one set of results, since clustering step is the same for all 'diffcyt' methods
load(file.path(DIR_RDATA, "outputs_BCR_XL_sim_diffcyt_DS_limma_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_BCR_XL_sim_diffcyt_DS_limma_main.RData"))
load(file.path(DIR_RDATA, "out_objects_BCR_XL_sim_diffcyt_DS_limma_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/BCR_XL_sim/main_clustering_performance"




##################################
# Calculate clustering performance
##################################

# ----------------------------------------
# load data objects and true B-cell status
# ----------------------------------------

d_se <- out_objects_diffcyt_DS_limma_main$d_se

B_cell <- out_diffcyt_DS_limma_main$B_cell

stopifnot(length(B_cell) == nrow(d_se))

# add true B cell status to data object
rowData(d_se)$B_cell <- B_cell


# -------------------------------------------------------------------------
# calculate clustering performance for all clusters containing true B cells
# -------------------------------------------------------------------------

# find matching clusters (clusters containing true B cells)

# check no missing clusters
stopifnot(all(names(table(rowData(d_se)[rowData(d_se)$B_cell == 1, ]$cluster_id)) == levels(rowData(d_se)$cluster_id)))

labels_matched <- unname(which(table(rowData(d_se)[rowData(d_se)$B_cell == 1, ]$cluster_id) > 0))
labels_matched

# total number of cells in each matching cluster
n_matched <- sapply(labels_matched, function(l) sum(rowData(d_se)$cluster_id == l))
n_matched

# number of true B cells in each matching cluster
n_correct <- sapply(labels_matched, function(l) sum(rowData(d_se)$cluster_id == l & rowData(d_se)$B_cell == 1))
n_correct

# total number of true B cells
n_B_cells <- sum(rowData(d_se)$B_cell == 1)
n_B_cells

# calculate precision, recall, F1 score for each matching cluster

stopifnot(length(n_matched) == length(n_correct), 
          length(n_B_cells) == 1)

pr <- n_correct / n_matched
re <- n_correct / n_B_cells
F1 <- 2 * (pr * re) / (pr + re)




################
# Generate plots
################

# create data frame for plotting

d_plot <- data.frame(
  cluster = labels_matched, 
  precision = pr, 
  recall = re, 
  F1_score = F1
)

# sort by F1 score
d_plot <- d_plot[rev(order(d_plot$F1_score)), ]

d_plot$cluster <- factor(d_plot$cluster, levels = as.character(d_plot$cluster))
d_plot <- melt(d_plot, id.vars = "cluster", variable.name = "measure")
d_plot$measure <- factor(d_plot$measure, levels = c("F1_score", "precision", "recall"))


# create plot

colors <- c("firebrick1", "forestgreen", "deepskyblue")

p <- 
  ggplot(d_plot, aes(x = cluster, y = value, color = measure)) + 
  geom_point(shape = 1, stroke = 1) + 
  scale_color_manual(values = colors) + 
  ylim(c(-0.025, 1.025)) + 
  ggtitle("BCR-XL-sim, diffcyt methods: clustering performance") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        axis.title.y = element_blank())

# save plot
fn <- file.path(DIR_PLOTS, "results_BCR_XL_sim_diffcyt_main_clustering_performance.pdf")
ggsave(fn, width = 6, height = 2.75)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



