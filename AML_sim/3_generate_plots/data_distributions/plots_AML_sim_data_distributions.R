##########################################################################################
# Generate plots
# 
# - data set: AML-sim
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


DIR_BENCHMARK_MAIN <- "../../../../../benchmark_data/AML_sim/data/main"
DIR_BENCHMARK_LESS_DISTINCT <- "../../../../../benchmark_data/AML_sim/data/less_distinct"

DIR_PLOTS <- "../../../../plots/AML_sim/data_distributions"




#################
# Load data: main
#################

# ----------------------
# Healthy blasts (H1-H5)
# ----------------------

# column indices of lineage markers used for clustering (see data preparation scripts)
cols_lineage <- c(35, 29, 14, 30, 12, 26, 17, 33, 41, 32, 22, 40, 27, 37, 23, 39)

# filenames
files_blasts_healthy <- list.files(file.path(DIR_BENCHMARK_MAIN, "blasts_all"), 
                                   pattern = "_H[0-9]+\\.fcs$", full.names = TRUE)

sample_names <- as.list(paste0("H", 1:5))

cofactor <- 5

data_H_samples <- 
  # load expression values
  files_blasts_healthy %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select lineage markers used for clustering
  lapply(function(d) d[, cols_lineage]) %>% 
  # arcsinh transform
  lapply(function(d) asinh(d / cofactor)) %>%
  # add sample labels
  mapply(function(d, s) data.frame(d, sample = s), ., sample_names, SIMPLIFY = FALSE) %>% 
  do.call(rbind, .)

colnames(data_H_samples) <- gsub("\\..*$", "", colnames(data_H_samples))
colnames(data_H_samples)[colnames(data_H_samples) == "HLA"] <- "HLA-DR"



# ---------
# CN blasts
# ---------

# filenames
files_blasts_CN <- list.files(file.path(DIR_BENCHMARK_MAIN, "blasts_all"), 
                              pattern = "_CN\\.fcs$", full.names = TRUE)

data_CN <- 
  # load expression values
  files_blasts_CN %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select lineage markers used for clustering
  lapply(function(d) d[, cols_lineage]) %>% 
  # arcsinh transform
  lapply(function(d) asinh(d / cofactor)) %>% 
  .[[1]]

colnames(data_CN) <- gsub("\\(.*$", "", colnames(data_CN))



# ----------
# CBF blasts
# ----------

# filenames
files_blasts_CBF <- list.files(file.path(DIR_BENCHMARK_MAIN, "blasts_all"), 
                               pattern = "_CBF\\.fcs$", full.names = TRUE)

data_CBF <- 
  # load expression values
  files_blasts_CBF %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select lineage markers used for clustering
  lapply(function(d) d[, cols_lineage]) %>% 
  # arcsinh transform
  lapply(function(d) asinh(d / cofactor)) %>% 
  .[[1]]

colnames(data_CBF) <- gsub("\\(.*$", "", colnames(data_CBF))




##############################################
# Load data: 'less distinct' blast populations
##############################################

# load data from 5% threshold (since this contains the largest number of cells), and
# combine cells from spike-in samples H1-H5 for each condition


# ---------
# CN blasts
# ---------

# 50% less distinct

files_blasts_50pc_CN <- list.files(file.path(DIR_BENCHMARK_LESS_DISTINCT, "less_50pc/CN"), 
                                   pattern = "^AML_sim_CN_H[0-9]+_5pc", full.names = TRUE)

data_50pc_CN_samples <- 
  # load expression values
  files_blasts_50pc_CN %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select spike-in cells only
  lapply(function(d) d[d[, "spikein"] == 1, ]) %>% 
  # select lineage markers used for clustering
  lapply(function(d) d[, cols_lineage]) %>% 
  # arcsinh transform
  lapply(function(d) asinh(d / cofactor)) %>%
  # add sample labels
  mapply(function(d, s) data.frame(d, sample = s), ., sample_names, SIMPLIFY = FALSE) %>% 
  do.call(rbind, .)

colnames(data_50pc_CN_samples) <- gsub("\\..*$", "", colnames(data_50pc_CN_samples))
colnames(data_50pc_CN_samples)[colnames(data_50pc_CN_samples) == "HLA"] <- "HLA-DR"


# 75% less distinct

files_blasts_75pc_CN <- list.files(file.path(DIR_BENCHMARK_LESS_DISTINCT, "less_75pc/CN"), 
                                   pattern = "^AML_sim_CN_H[0-9]+_5pc", full.names = TRUE)

data_75pc_CN_samples <- 
  # load expression values
  files_blasts_75pc_CN %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select spike-in cells only
  lapply(function(d) d[d[, "spikein"] == 1, ]) %>% 
  # select lineage markers used for clustering
  lapply(function(d) d[, cols_lineage]) %>% 
  # arcsinh transform
  lapply(function(d) asinh(d / cofactor)) %>%
  # add sample labels
  mapply(function(d, s) data.frame(d, sample = s), ., sample_names, SIMPLIFY = FALSE) %>% 
  do.call(rbind, .)

colnames(data_75pc_CN_samples) <- gsub("\\..*$", "", colnames(data_75pc_CN_samples))
colnames(data_75pc_CN_samples)[colnames(data_75pc_CN_samples) == "HLA"] <- "HLA-DR"



# ----------
# CBF blasts
# ----------

# 50% less distinct

files_blasts_50pc_CBF <- list.files(file.path(DIR_BENCHMARK_LESS_DISTINCT, "less_50pc/CBF"), 
                                    pattern = "^AML_sim_CBF_H[0-9]+_5pc", full.names = TRUE)

data_50pc_CBF_samples <- 
  # load expression values
  files_blasts_50pc_CBF %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select spike-in cells only
  lapply(function(d) d[d[, "spikein"] == 1, ]) %>% 
  # select lineage markers used for clustering
  lapply(function(d) d[, cols_lineage]) %>% 
  # arcsinh transform
  lapply(function(d) asinh(d / cofactor)) %>%
  # add sample labels
  mapply(function(d, s) data.frame(d, sample = s), ., sample_names, SIMPLIFY = FALSE) %>% 
  do.call(rbind, .)

colnames(data_50pc_CBF_samples) <- gsub("\\..*$", "", colnames(data_50pc_CBF_samples))
colnames(data_50pc_CBF_samples)[colnames(data_50pc_CBF_samples) == "HLA"] <- "HLA-DR"


# 75% less distinct

files_blasts_75pc_CBF <- list.files(file.path(DIR_BENCHMARK_LESS_DISTINCT, "less_75pc/CBF"), 
                                    pattern = "^AML_sim_CBF_H[0-9]+_5pc", full.names = TRUE)

data_75pc_CBF_samples <- 
  # load expression values
  files_blasts_75pc_CBF %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select spike-in cells only
  lapply(function(d) d[d[, "spikein"] == 1, ]) %>% 
  # select lineage markers used for clustering
  lapply(function(d) d[, cols_lineage]) %>% 
  # arcsinh transform
  lapply(function(d) asinh(d / cofactor)) %>%
  # add sample labels
  mapply(function(d, s) data.frame(d, sample = s), ., sample_names, SIMPLIFY = FALSE) %>% 
  do.call(rbind, .)

colnames(data_75pc_CBF_samples) <- gsub("\\..*$", "", colnames(data_75pc_CBF_samples))
colnames(data_75pc_CBF_samples)[colnames(data_75pc_CBF_samples) == "HLA"] <- "HLA-DR"




################
# Generate plots
################

# -------------------------------------------
# Ridges plots to compare samples H1-H5: main
# -------------------------------------------

d_plot <- melt(data_H_samples, id.vars = "sample", variable.name = "marker", value.name = "expression")

d_plot$marker <- factor(d_plot$marker, levels = levels(d_plot$marker)[rev(order(levels(d_plot$marker)))])


p <- ggplot(d_plot, aes(x = expression, y = marker)) + 
  geom_density_ridges(fill = "royalblue", alpha = 0.5, lwd = 0.25) + 
  facet_grid(~ sample) + 
  xlab("arcsinh transformed expression") + 
  ylab("density") + 
  theme_bw() + 
  ggtitle("AML-sim: marker distributions for blast cells, healthy samples")


# save plot
filename <- file.path(DIR_PLOTS, "AML_sim_data_distributions_healthy_samples.pdf")
ggsave(filename, width = 6, height = 7.25)



# --------------------------------------------------
# Ridges plots to compare healthy, CN, and CBF: main
# --------------------------------------------------

# combine all samples H1-H5 for healthy

data_H <- data_H_samples[, -which(colnames(data_H_samples) == "sample")]

d_plot_H <- cbind(data_H, condition = "healthy")
d_plot_CN <- cbind(as.data.frame(data_CN), condition = "CN")
d_plot_CBF <- cbind(as.data.frame(data_CBF), condition = "CBF")

stopifnot(all(colnames(d_plot_H) == colnames(d_plot_CN)), all(colnames(d_plot_H) == colnames(d_plot_CBF)))

d_plot_cnd <- rbind(d_plot_H, d_plot_CN, d_plot_CBF)

d_plot <- melt(d_plot_cnd, id.vars = "condition", variable.name = "marker", value.name = "expression")

d_plot$condition <- factor(d_plot$condition, levels = c("healthy", "CN", "CBF"))
d_plot$marker <- factor(d_plot$marker, levels = levels(d_plot$marker)[rev(order(levels(d_plot$marker)))])


# color scheme
colors <- c("royalblue", "gold", "tomato")


p <- ggplot(d_plot, aes(x = expression, y = marker, fill = condition)) + 
  geom_density_ridges(alpha = 0.5, lwd = 0.25) + 
  scale_fill_cyclical(values = colors, guide = "legend") + 
  xlab("arcsinh transformed expression") + 
  ylab("density") + 
  theme_bw() + 
  ggtitle("AML-sim: blast cells")


# save plot
filename <- file.path(DIR_PLOTS, "AML_sim_data_distributions_conditions_main.pdf")
ggsave(filename, width = 4, height = 7)



# -------------------------------------------------------------------------------
# Ridges plots to compare healthy, CN, and CBF: 'less distinct' blast populations
# -------------------------------------------------------------------------------

# combine all samples H1-H5

data_50pc_CN <- data_50pc_CN_samples[, -which(colnames(data_50pc_CN_samples) == "sample")]
data_50pc_CBF <- data_50pc_CBF_samples[, -which(colnames(data_50pc_CBF_samples) == "sample")]
data_75pc_CN <- data_75pc_CN_samples[, -which(colnames(data_75pc_CN_samples) == "sample")]
data_75pc_CBF <- data_75pc_CBF_samples[, -which(colnames(data_75pc_CBF_samples) == "sample")]

data_50pc_CN <- cbind(data_50pc_CN, condition = "CN", less_distinct = "less_50pc")
data_50pc_CBF <- cbind(data_50pc_CBF, condition = "CBF", less_distinct = "less_50pc")
data_75pc_CN <- cbind(data_75pc_CN, condition = "CN", less_distinct = "less_75pc")
data_75pc_CBF <- cbind(data_75pc_CBF, condition = "CBF", less_distinct = "less_75pc")

d_plot_cnd_main <- cbind(d_plot_cnd, less_distinct = "main")

# required for facet_grid()
d_plot_H_50pc <- cbind(d_plot_cnd[d_plot_cnd$condition == "healthy", ], less_distinct = "less_50pc")
d_plot_H_75pc <- cbind(d_plot_cnd[d_plot_cnd$condition == "healthy", ], less_distinct = "less_75pc")

stopifnot(all(colnames(data_50pc_CN) == colnames(d_plot_cnd_main)), 
          all(colnames(data_50pc_CBF) == colnames(d_plot_cnd_main)), 
          all(colnames(data_75pc_CN) == colnames(d_plot_cnd_main)), 
          all(colnames(data_75pc_CBF) == colnames(d_plot_cnd_main)), 
          all(colnames(d_plot_H_50pc) == colnames(d_plot_cnd_main)), 
          all(colnames(d_plot_H_75pc) == colnames(d_plot_cnd_main)))

d_plot <- rbind(d_plot_cnd_main, data_50pc_CN, data_50pc_CBF, data_75pc_CN, data_75pc_CBF, d_plot_H_50pc, d_plot_H_75pc)

# check with contingency table
table(d_plot$condition, d_plot$less_distinct)


d_plot$condition <- factor(d_plot$condition, levels = c("healthy", "CN", "CBF"))
d_plot$less_distinct <- factor(d_plot$less_distinct, levels = c("main", "less_50pc", "less_75pc"))

d_plot <- melt(d_plot, id.vars = c("condition", "less_distinct"), variable.name = "marker", value.name = "expression")

d_plot$marker <- factor(d_plot$marker, levels = levels(d_plot$marker)[rev(order(levels(d_plot$marker)))])


p <- ggplot(d_plot, aes(x = expression, y = marker, fill = condition)) + 
  geom_density_ridges(alpha = 0.5, lwd = 0.25) + 
  scale_fill_cyclical(values = colors, guide = "legend") + 
  facet_grid(~ less_distinct) + 
  xlim(c(-2, 8)) + 
  xlab("arcsinh transformed expression") + 
  ylab("density") + 
  theme_bw() + 
  ggtitle("AML-sim: main simulation and 'less distinct' populations")


# save plot
filename <- file.path(DIR_PLOTS, "AML_sim_data_distributions_conditions_less_distinct.pdf")
ggsave(filename, width = 7, height = 7)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



