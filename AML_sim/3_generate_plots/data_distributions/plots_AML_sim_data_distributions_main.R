##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: density plots
# 
# - data distributions for main benchmark data
# 
# Lukas Weber, October 2017
##########################################################################################


library(ggplot2)
library(ggridges)
library(flowCore)
library(magrittr)
library(reshape2)


DIR_BENCHMARK <- "../../../../../benchmark_data/AML_sim/data/main"

DIR_PLOTS <- "../../../../plots/AML_sim/data_distributions"




###########
# Load data
###########

# ----------------------
# Healthy blasts (H1-H5)
# ----------------------

# column indices of lineage markers used for clustering (see data preparation scripts)
cols_lineage <- c(35, 29, 14, 30, 12, 26, 17, 33, 41, 32, 22, 40, 27, 37, 23, 39)

# filenames
files_blasts_healthy <- list.files(file.path(DIR_BENCHMARK, "blasts_all"), pattern = "_H[0-9]+\\.fcs$", full.names = TRUE)

sample_names <- as.list(paste0("H", 1:5))

cofactor <- 5

data_H_samples <- 
  # load expression values
  files_blasts_healthy %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select lineage markers used for clustering
  lapply(function(d) d[, cols_lineage]) %>% 
  # asinh transform
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
files_blasts_CN <- list.files(file.path(DIR_BENCHMARK, "blasts_all"), pattern = "_CN\\.fcs$", full.names = TRUE)

data_CN <- 
  # load expression values
  files_blasts_CN %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select lineage markers used for clustering
  lapply(function(d) d[, cols_lineage]) %>% 
  # asinh transform
  lapply(function(d) asinh(d / cofactor)) %>% 
  .[[1]]

colnames(data_CN) <- gsub("\\(.*$", "", colnames(data_CN))



# ----------
# CBF blasts
# ----------

# filenames
files_blasts_CBF <- list.files(file.path(DIR_BENCHMARK, "blasts_all"), pattern = "_CBF\\.fcs$", full.names = TRUE)

data_CBF <- 
  # load expression values
  files_blasts_CBF %>% 
  lapply(read.FCS, transformation = FALSE, truncate_max_range = FALSE) %>% 
  lapply(exprs) %>% 
  # select lineage markers used for clustering
  lapply(function(d) d[, cols_lineage]) %>% 
  # asinh transform
  lapply(function(d) asinh(d / cofactor)) %>% 
  .[[1]]

colnames(data_CBF) <- gsub("\\(.*$", "", colnames(data_CBF))




################
# Generate plots
################

# -------------------------------------
# Ridges plots to compare samples H1-H5
# -------------------------------------

d_plot <- melt(data_H_samples, id.vars = "sample", variable.name = "marker", value.name = "expression")

d_plot$marker <- factor(d_plot$marker, levels = levels(d_plot$marker)[rev(order(levels(d_plot$marker)))])


p <- ggplot(d_plot, aes(x = expression, y = marker)) + 
  geom_density_ridges(fill = "royalblue", alpha = 0.5) + 
  facet_grid(~ sample) + 
  xlab("asinh transformed expression") + 
  ylab("density per marker") + 
  theme_bw() + 
  ggtitle("AML-sim, marker distributions, healthy samples")


# save plot
filename <- file.path(DIR_PLOTS, "AML_sim_data_distributions_healthy_samples.pdf")
ggsave(filename, width = 6, height = 6.5)



# --------------------------------------------
# Ridges plots to compare healthy, CN, and CBF
# --------------------------------------------

# combine all samples H1-H5 for healthy

data_H <- data_H_samples[, -which(colnames(data_H_samples) == "sample")]

d_plot_H <- cbind(data_H, condition = "healthy")
d_plot_CN <- cbind(as.data.frame(data_CN), condition = "CN")
d_plot_CBF <- cbind(as.data.frame(data_CBF), condition = "CBF")

stopifnot(all(colnames(d_plot_H) == colnames(d_plot_CN)), all(colnames(d_plot_H) == colnames(d_plot_CBF)))

d_plot <- rbind(d_plot_H, d_plot_CN, d_plot_CBF)

d_plot <- melt(d_plot, id.vars = "condition", variable.name = "marker", value.name = "expression")

d_plot$condition <- factor(d_plot$condition, levels = c("healthy", "CN", "CBF"))
d_plot$marker <- factor(d_plot$marker, levels = levels(d_plot$marker)[rev(order(levels(d_plot$marker)))])


# color scheme
colors <- c("royalblue", "gold", "tomato")


p <- ggplot(d_plot, aes(x = expression, y = marker, fill = condition)) + 
  geom_density_ridges(alpha = 0.5) + 
  scale_fill_cyclical(values = colors, guide = "legend") + 
  xlab("asinh transformed expression") + 
  ylab("density per marker") + 
  theme_bw() + 
  ggtitle("AML-sim, marker distributions")


# save plot
filename <- file.path(DIR_PLOTS, "AML_sim_data_distributions_by_condition.pdf")
ggsave(filename, width = 4.5, height = 7)




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



