##########################################################################################
# Script to prepare benchmark data set 'BCR-XL-sim'
# 
# The original 'BCR-XL' data set is sourced from Bodenmiller et al. (2012), and was
# previously used for benchmark evaluations by Bruggner et al. (2014) (Citrus paper).
# 
# Raw data downloaded from Cytobank (experiment 15713)
# - see Citrus wiki (section "PBMC Example 1"):
# https://github.com/nolanlab/citrus/wiki/PBMC-Example-1
# - direct link to Cytobank repository:
# https://community.cytobank.org/cytobank/experiments/15713/download_files
# 
# Cell population labels are reproduced from Nowicka et al. (2017), F1000Research, using a
# strategy of expert-guided manual merging of automatically generated clusters from the
# FlowSOM algorithm. Code to reproduce the cell population labels is available in the
# script "cell_population_labels_BCR_XL.R".
# 
# The simulations in this script are generated as follows:
# - select reference (unstimulated) samples from the main 'BCR-XL' data set
# - randomly split each sample into two halves
# - in one half, replace B cells with equivalent number of B cells from the corresponding
# stimulated sample
# - adjust "difficulty" of the simulation by scaling the average difference in pS6 signal
# 
# Methods are then evaluated by their ability to detect the known strong differential
# signal in pS6 expression.
# 
# Lukas Weber, August 2017
##########################################################################################


library(flowCore)




###########
# FILENAMES
###########

# .fcs files
DIR_RAW_DATA <- "../../../../benchmark_data/BCR_XL_sim/raw_data/experiment_15713_files"
files <- list.files(DIR_RAW_DATA, pattern = "\\.fcs$", full.names = TRUE)

files_BCRXL <- files[grep("patient[1-8]_BCR-XL\\.fcs$", files)]
files_ref <- files[grep("patient[1-8]_Reference\\.fcs$", files)]

files_all <- c(files_BCRXL, files_ref)

data <- lapply(files_all, function(f) exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)))


# vector of conditions
conditions <- as.factor(gsub("\\.fcs$", "", gsub("^.*_", "", files_all)))
conditions

# vector of patient IDs
patient_IDs <- as.factor(gsub("_.*$", "", gsub("^PBMC8_30min_", "", basename(files_all))))
patient_IDs


# cell population labels
DIR_LABELS <- "../../../../benchmark_data/BCR_XL_sim/population_IDs"
files_labels <- list.files(DIR_LABELS, pattern = "\\.csv$", full.names = TRUE)

files_labels_BCRXL <- files_labels[grep("patient[1-8]_BCR-XL\\.csv$", files_labels)]
files_labels_ref <- files_labels[grep("patient[1-8]_Reference\\.csv$", files_labels)]

files_labels_all <- c(files_labels_BCRXL, files_labels_ref)

labels <- lapply(files_labels_all, read.csv)




#############################
# EXPORT FILES: NO SIMULATION
#############################

# 'nosim': export data files without any modifications

DIR_DATA_NOSIM <- "../../../../benchmark_data/BCR_XL_sim/data/nosim"

cmds <- paste("cp", files_all, DIR_DATA_NOSIM)

for (cmd in cmds) {
  system(cmd)
}




##################
# SIMULATION: FULL
##################

# 'full' simulation contains the full differential expression signal for pS6


# add column of population labels for each sample
data <- mapply(function(d, l) {
  stopifnot(nrow(d) == nrow(l))
  d <- as.data.frame(d)
  d <- cbind(d, l)
  d
}, data, labels, SIMPLIFY = FALSE)

# add column indicating B cells for each sample
data <- lapply(data, function(d) {
  is_B_cell <- as.numeric(d$population %in% c("B-cells IgM-", "B-cells IgM+"))
  d <- cbind(d, is_B_cell)
  d
})


# ---------------------------------------------------------------
# split each Reference sample into two halves: 'base' and 'spike'
# ---------------------------------------------------------------

data_ref <- data[conditions == "Reference"]

n_cells_ref <- sapply(data_ref, nrow)

set.seed(123)

# generate random indices
inds <- lapply(n_cells_ref, function(n) {
  i_base <- sort(sample(seq_len(n), floor(n / 2)))
  i_spike <- setdiff(seq_len(n), i_base)
  list(base = i_base, spike = i_spike)
})

inds_base <- lapply(inds, function(l) l[[1]])
inds_spike <- lapply(inds, function(l) l[[2]])

# subset data
data_base <- mapply(function(d, i) d[i, ], data_ref, inds_base, SIMPLIFY = FALSE)
data_spike <- mapply(function(d, i) d[i, ], data_ref, inds_spike, SIMPLIFY = FALSE)


# ------------------------------------------------------------------------------------
# replace B cells in 'spike' samples with an equivalent number of B cells from 'BCRXL'
# (stimulated) condition
# ------------------------------------------------------------------------------------

# note: for some samples, there may not be enough B cells available; simply use all
# available B cells in this case

# B cells from 'BCRXL' (stimulated) condition
data_BCRXL <- data[conditions == "BCR-XL"]

B_cells_BCRXL <- lapply(data_BCRXL, function(d) {
  d[d$is_B_cell == 1, ]
})

# number of B cells available
sapply(B_cells_BCRXL, nrow)

# number of B cells needed
n_spike <- sapply(data_spike, function(d) {
  sum(d$is_B_cell == 1)
})
n_spike

# select correct number of B cells from 'BCRXL' (stimulated) condition for each sample

set.seed(123)

B_cells_spike <- mapply(function(b, n) {
  # reduce 'n' if not enough B cells available
  n <- min(n, nrow(b))
  inds <- sample(seq_len(nrow(b)), n)
  b[inds, ]
}, B_cells_BCRXL, n_spike, SIMPLIFY = FALSE)

sapply(B_cells_spike, nrow)

# replace B cells in 'spike' samples

data_spike <- mapply(function(d, b) {
  d <- d[d$is_B_cell == 0, ]
  d <- rbind(d, b)
  rownames(d) <- NULL
  d
}, data_spike, B_cells_spike, SIMPLIFY = FALSE)


# ------------
# export files
# ------------

DIR_DATA_FULL <- "../../../../benchmark_data/BCR_XL_sim/data/sim_full"

data_export <- c(data_spike, data_base)

conditions_spike <- c(rep("spike", length(data_spike)), rep("base", length(data_base)))

# add column indicating spike-in cells (all B cells in 'spike' samples)
data_export <- mapply(function(d, cnd) {
  is_spikein <- as.numeric((d$is_B_cell == 1) & (cnd == "spike"))
  d$is_spikein <- is_spikein
  d
}, data_export, conditions_spike, SIMPLIFY = FALSE)


# filenames
filenames <- file.path(DIR_DATA_FULL, paste0("BCR_XL_sim_", patient_IDs, "_", conditions_spike, ".fcs"))

for (i in seq_along(data_export)) {
  
  # convert to flowFrame
  d_export_i <- data_export[[i]]
  d_export_i$population <- as.numeric(d_export_i$population)
  d_export_i <- flowFrame(as.matrix(d_export_i))
  
  # save as .fcs file
  write.FCS(d_export_i, filename = filenames[i])
}




#################
# SIMULATION: 0.5
#################




##################
# SIMULATION: 0.25
##################



