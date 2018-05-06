##########################################################################################
# Generate plots
# 
# - data set: BCR-XL-sim
# - plot type: density plots
# 
# - data distributions for main benchmark data
# 
# Lukas Weber, May 2018
##########################################################################################


library(ggplot2)
library(ggridges)
library(flowCore)
library(magrittr)
library(reshape2)
library(SummarizedExperiment)


DIR_BENCHMARK_MAIN <- "../../../../../benchmark_data/BCR_XL_sim/data/main"
DIR_BENCHMARK_LESS_DISTINCT <- "../../../../../benchmark_data/BCR_XL_sim/data/less_distinct"

DIR_RDATA <- "../../../../RData/BCR_XL_sim/main"
load(file.path(DIR_RDATA, "out_objects_BCR_XL_sim_diffcyt_DS_limma_main.RData"))

DIR_PLOTS <- "../../../../plots/BCR_XL_sim/data_distributions"




#################
# Load data: main
#################

# filenames
files <- list.files(DIR_BENCHMARK_MAIN, pattern = "\\.fcs$", full.names = TRUE)
files_base <- files[grep("base\\.fcs$", files)]
files_spike <- files[grep("spike\\.fcs$", files)]

# indices of all marker columns, lineage markers, and functional markers
# (10 surface markers / 14 functional markers; see Bruggner et al. 2014, Table 1)
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)

# load objects (to identify cell type and cell state markers)
d_medians_by_cluster_marker <- out_objects_diffcyt_DS_limma_main$d_medians_by_cluster_marker


# -----------------------
# create data frame: base
# -----------------------

cofactor <- 5

d_base <- 
  # load expression values
  files_base %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select all markers
  lapply(function(d) d[, cols_markers]) %>% 
  # arcsinh transform
  lapply(function(d) asinh(d / cofactor)) %>%
  do.call(rbind, .)

# arrange cell type and cell state markers in two groups
d_base <- cbind(d_base[, metadata(d_medians_by_cluster_marker)$id_type_markers], 
                d_base[, metadata(d_medians_by_cluster_marker)$id_state_markers])

# arrange each group (cell type and cell state markers) alphabetically
n_type <- sum(metadata(d_medians_by_cluster_marker)$id_type_markers)
n_state <- sum(metadata(d_medians_by_cluster_marker)$id_state_markers)
d_base <- cbind(d_base[, seq_len(n_type)][, order(colnames(d_base)[seq_len(n_type)])], 
                d_base[, (seq_len(n_state)) + n_type][, order(colnames(d_base)[(seq_len(n_state)) + n_type])])

colnames(d_base) <- gsub("\\(.*$", "", colnames(d_base))

# labels for B cells
d_base_labels <- 
  # load expression values
  files_base %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select labels
  lapply(function(d) d[, colnames(d) == "B_cell", drop = FALSE]) %>% 
  do.call(rbind, .)

stopifnot(nrow(d_base) == nrow(d_base_labels))

d_base <- cbind(d_base, d_base_labels)


# ------------------------
# create data frame: spike
# ------------------------

cofactor <- 5

d_spike <- 
  # load expression values
  files_spike %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select all markers
  lapply(function(d) d[, cols_markers]) %>% 
  # arcsinh transform
  lapply(function(d) asinh(d / cofactor)) %>%
  do.call(rbind, .)

# arrange cell type and cell state markers in two groups
d_spike <- cbind(d_spike[, metadata(d_medians_by_cluster_marker)$id_type_markers], 
                 d_spike[, metadata(d_medians_by_cluster_marker)$id_state_markers])

# arrange each group (cell type and cell state markers) alphabetically
n_type <- sum(metadata(d_medians_by_cluster_marker)$id_type_markers)
n_state <- sum(metadata(d_medians_by_cluster_marker)$id_state_markers)
d_spike <- cbind(d_spike[, seq_len(n_type)][, order(colnames(d_spike)[seq_len(n_type)])], 
                 d_spike[, (seq_len(n_state)) + n_type][, order(colnames(d_spike)[(seq_len(n_state)) + n_type])])

colnames(d_spike) <- gsub("\\(.*$", "", colnames(d_spike))

# labels for B cells
d_spike_labels <- 
  # load expression values
  files_spike %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select labels
  lapply(function(d) d[, colnames(d) == "B_cell", drop = FALSE]) %>% 
  do.call(rbind, .)

stopifnot(nrow(d_spike) == nrow(d_spike_labels))

d_spike <- cbind(d_spike, d_spike_labels)


# check column names
stopifnot(all(colnames(d_base) == colnames(d_spike)))

names_type <- colnames(d_base)[1:n_type]
names_state <- colnames(d_base)[(n_type + 1):(n_type + n_state)]




###########################################
# Load data: 'less distinct' spike-in cells
###########################################

# -----------------
# 50% less distinct
# -----------------

files_50pc <- list.files(file.path(DIR_BENCHMARK_LESS_DISTINCT, "less_50pc"), pattern = "\\.fcs$", full.names = TRUE)
files_spike_50pc <- files_50pc[grep("spike", files_50pc)]

# create data frame
cofactor <- 5

d_spike_50pc <- 
  # load expression values
  files_spike_50pc %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select all markers
  lapply(function(d) d[, cols_markers]) %>% 
  # arcsinh transform
  lapply(function(d) asinh(d / cofactor)) %>%
  do.call(rbind, .)

# arrange cell type and cell state markers in two groups
d_spike_50pc <- cbind(d_spike_50pc[, metadata(d_medians_by_cluster_marker)$id_type_markers], 
                      d_spike_50pc[, metadata(d_medians_by_cluster_marker)$id_state_markers])

# arrange each group (cell type and cell state markers) alphabetically
n_type <- sum(metadata(d_medians_by_cluster_marker)$id_type_markers)
n_state <- sum(metadata(d_medians_by_cluster_marker)$id_state_markers)
d_spike_50pc <- cbind(d_spike_50pc[, seq_len(n_type)][, order(colnames(d_spike_50pc)[seq_len(n_type)])], 
                      d_spike_50pc[, (seq_len(n_state)) + n_type][, order(colnames(d_spike_50pc)[(seq_len(n_state)) + n_type])])

colnames(d_spike_50pc) <- gsub("\\(.*$", "", colnames(d_spike_50pc))

# labels for B cells
d_spike_50pc_labels <- 
  # load expression values
  files_spike_50pc %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select labels
  lapply(function(d) d[, colnames(d) == "B_cell", drop = FALSE]) %>% 
  do.call(rbind, .)

stopifnot(nrow(d_spike_50pc) == nrow(d_spike_50pc_labels))

d_spike_50pc <- cbind(d_spike_50pc, d_spike_50pc_labels)


# -----------------
# 75% less distinct
# -----------------

files_75pc <- list.files(file.path(DIR_BENCHMARK_LESS_DISTINCT, "less_75pc"), pattern = "\\.fcs$", full.names = TRUE)
files_spike_75pc <- files_75pc[grep("spike", files_75pc)]

# create data frame
cofactor <- 5

d_spike_75pc <- 
  # load expression values
  files_spike_75pc %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select all markers
  lapply(function(d) d[, cols_markers]) %>% 
  # arcsinh transform
  lapply(function(d) asinh(d / cofactor)) %>%
  do.call(rbind, .)

# arrange cell type and cell state markers in two groups
d_spike_75pc <- cbind(d_spike_75pc[, metadata(d_medians_by_cluster_marker)$id_type_markers], 
                      d_spike_75pc[, metadata(d_medians_by_cluster_marker)$id_state_markers])

# arrange each group (cell type and cell state markers) alphabetically
n_type <- sum(metadata(d_medians_by_cluster_marker)$id_type_markers)
n_state <- sum(metadata(d_medians_by_cluster_marker)$id_state_markers)
d_spike_75pc <- cbind(d_spike_75pc[, seq_len(n_type)][, order(colnames(d_spike_75pc)[seq_len(n_type)])], 
                      d_spike_75pc[, (seq_len(n_state)) + n_type][, order(colnames(d_spike_75pc)[(seq_len(n_state)) + n_type])])

colnames(d_spike_75pc) <- gsub("\\(.*$", "", colnames(d_spike_75pc))

# labels for B cells
d_spike_75pc_labels <- 
  # load expression values
  files_spike_75pc %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select labels
  lapply(function(d) d[, colnames(d) == "B_cell", drop = FALSE]) %>% 
  do.call(rbind, .)

stopifnot(nrow(d_spike_75pc) == nrow(d_spike_75pc_labels))

d_spike_75pc <- cbind(d_spike_75pc, d_spike_75pc_labels)




################
# Generate plots
################

# -----------------------------------------------------
# Ridges plots to compare B cells in 'base' vs. 'spike'
# -----------------------------------------------------

# prepare data frame for plotting
d_base <- cbind(as.data.frame(d_base), group = "base")
d_spike <- cbind(as.data.frame(d_spike), group = "spike")

d_plot <- rbind(d_base, d_spike)

# keep B cells only
d_plot <- d_plot[d_plot$B_cell == 1, -which(colnames(d_plot) == "B_cell")]

d_plot <- melt(d_plot, id.vars = "group", variable.name = "marker", value.name = "expression")

# marker types
d_plot$marker_class <- factor(as.numeric(d_plot$marker %in% names_state), labels = c("cell type", "cell state"))

d_plot$marker <- factor(d_plot$marker, levels = rev(levels(d_plot$marker)))

# color scheme
colors <- c("royalblue", "tomato")


# create plot
p <- 
  ggplot(d_plot, aes(x = expression, y = marker, fill = group)) + 
  geom_density_ridges(alpha = 0.5, lwd = 0.25) + 
  geom_point(aes(shape = NA, color = marker_class)) + 
  scale_fill_cyclical(values = colors, guide = "legend") + 
  xlab("arcsinh transformed expression") + 
  ylab("density") + 
  theme_bw() + 
  ggtitle("BCR-XL-sim: B cells") + 
  guides(fill = guide_legend(order = 1, 
                             override.aes = list(shape = NA)), 
         color = guide_legend(order = 2, title = "marker type", 
                              override.aes = list(shape = 15, size = 6, color = c("gold", "forestgreen"))))

p + 
  annotate("rect", xmin = -4.2, xmax = -3.6, ymin = 14.67, ymax = 24.33, fill = "gold") + 
  annotate("rect", xmin = -4.2, xmax = -3.6, ymin = 0.67, ymax = 14.33, fill = "forestgreen")

# save plot
filename <- file.path(DIR_PLOTS, "BCR_XL_sim_data_distributions_conditions_main.pdf")
ggsave(filename, width = 4, height = 7)



# ----------------------------------------------------------------------
# Ridges plots to compare B cells in 'base' vs. 'spike': 'less distinct'
# ----------------------------------------------------------------------

d_base_main <- cbind(d_base, less_distinct = "main")
d_base_50pc <- cbind(d_base, less_distinct = "less_50pc")
d_base_75pc <- cbind(d_base, less_distinct = "less_75pc")

d_spike_main <- cbind(d_spike, less_distinct = "main")
d_spike_50pc <- cbind(as.data.frame(d_spike_50pc), group = "spike", less_distinct = "less_50pc")
d_spike_75pc <- cbind(as.data.frame(d_spike_75pc), group = "spike", less_distinct = "less_75pc")

stopifnot(all(colnames(d_base_main) == colnames(d_base_50pc)), 
          all(colnames(d_base_main) == colnames(d_base_75pc)), 
          all(colnames(d_base_main) == colnames(d_spike_main)), 
          all(colnames(d_base_main) == colnames(d_spike_50pc)), 
          all(colnames(d_base_main) == colnames(d_spike_75pc)))

# data frame for plotting
d_plot <- rbind(d_base_main, d_base_50pc, d_base_75pc, d_spike_main, d_spike_50pc, d_spike_75pc)

d_plot$less_distinct <- factor(d_plot$less_distinct, levels = c("main", "less_50pc", "less_75pc"))

# check with contingency table
table(d_plot$group, d_plot$less_distinct)


# keep B cells only
d_plot <- d_plot[d_plot$B_cell == 1, -which(colnames(d_plot) == "B_cell")]

table(d_plot$group, d_plot$less_distinct)

d_plot <- melt(d_plot, id.vars = c("group", "less_distinct"), variable.name = "marker", value.name = "expression")

# marker types
d_plot$marker_class <- factor(as.numeric(d_plot$marker %in% names_state), labels = c("cell type", "cell state"))

d_plot$marker <- factor(d_plot$marker, levels = rev(levels(d_plot$marker)))

# color scheme
colors <- c("royalblue", "tomato")


# create plot
p <- 
  ggplot(d_plot, aes(x = expression, y = marker, fill = group)) + 
  geom_density_ridges(alpha = 0.5, lwd = 0.25) + 
  geom_point(aes(shape = NA, color = marker_class)) + 
  scale_fill_cyclical(values = colors, guide = "legend") + 
  facet_grid(~ less_distinct) + 
  xlab("arcsinh transformed expression") + 
  ylab("density") + 
  theme_bw() + 
  ggtitle("BCR-XL-sim: main simulation and 'less distinct' populations") + 
  guides(fill = guide_legend(order = 1, 
                             override.aes = list(shape = NA)), 
         color = guide_legend(order = 2, title = "marker type", 
                              override.aes = list(shape = 15, size = 6, color = c("gold", "forestgreen"))))

p + 
  annotate("rect", xmin = -4.2, xmax = -3.4, ymin = 14.67, ymax = 24.33, fill = "gold") + 
  annotate("rect", xmin = -4.2, xmax = -3.4, ymin = 0.67, ymax = 14.33, fill = "forestgreen")

# save plot
filename <- file.path(DIR_PLOTS, "BCR_XL_sim_data_distributions_conditions_less_distinct.pdf")
ggsave(filename, width = 7, height = 7)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



