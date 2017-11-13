##########################################################################################
# Script to prepare benchmark data set 'BCR-XL-sim'
# 
# The original 'BCR-XL' data set is sourced from Bodenmiller et al. (2012), and was
# previously used for benchmark evaluations by Bruggner et al. (2014) (Citrus paper).
# 
# Raw data downloaded from Cytobank (experiment 15713)
# - see Citrus wiki (section 'PBMC Example 1'):
# https://github.com/nolanlab/citrus/wiki/PBMC-Example-1
# - direct link to Cytobank repository:
# https://community.cytobank.org/cytobank/experiments/15713/download_files
# 
# Cell population labels are reproduced from Nowicka et al. (2017), F1000Research, using a
# strategy of expert-guided manual merging of automatically generated clusters from the
# FlowSOM algorithm. Code to reproduce the cell population labels is available in the
# script 'cell_population_labels_BCR_XL.R'.
# 
# The simulations in this script are generated as follows:
# - select reference (unstimulated) samples from the main 'BCR-XL' data set
# - randomly split each sample into two halves
# - in one half, replace B cells with equivalent number of B cells from the corresponding
# stimulated sample
# - adjust 'difficulty' of the simulation by scaling the average difference in pS6 signal
# 
# Methods are then evaluated by their ability to detect the known strong differential
# signal in pS6 expression.
# 
# Lukas Weber, November 2017
##########################################################################################


# 'main' simulation: contains the full differential expression signal for pS6


library(flowCore)



# ---------
# Filenames
# ---------

# .fcs files

DIR_BENCHMARK <- "../../../../../benchmark_data"

DIR_RAW_DATA <- file.path(DIR_BENCHMARK, "BCR_XL_sim/raw_data/experiment_15713_files")

files <- list.files(DIR_RAW_DATA, pattern = "\\.fcs$", full.names = TRUE)

files_BCRXL <- files[grep("patient[1-8]_BCR-XL\\.fcs$", files)]
files_ref <- files[grep("patient[1-8]_Reference\\.fcs$", files)]

files_all <- c(files_BCRXL, files_ref)


# cell population labels

DIR_LABELS <- file.path(DIR_BENCHMARK, "BCR_XL_sim/population_IDs")

files_labels <- list.files(DIR_LABELS, pattern = "\\.csv$", full.names = TRUE)

files_labels_BCRXL <- files_labels[grep("patient[1-8]_BCR-XL\\.csv$", files_labels)]
files_labels_ref <- files_labels[grep("patient[1-8]_Reference\\.csv$", files_labels)]

files_labels_all <- c(files_labels_BCRXL, files_labels_ref)



# ---------
# Load data
# ---------

data <- lapply(files_all, function(f) exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)))

# sample IDs
sample_IDs <- gsub("\\.fcs$", "", gsub("^PBMC8_30min_", "", basename(files_all)))
sample_IDs

names(data) <- sample_IDs

# conditions
conditions <- as.factor(gsub("^patient[0-9+]_", "", sample_IDs))
conditions

# patient IDs
patient_IDs <- as.factor(gsub("_.*$", "", sample_IDs))
patient_IDs

# cell population labels
labels <- lapply(files_labels_all, read.csv)



# ---------------------
# Add population labels
# ---------------------

data <- mapply(function(d, l) {
  stopifnot(nrow(d) == nrow(l))
  stopifnot(levels(l$population) == levels(labels[[1]]$population))
  
  B_cell <- l$population %in% c("B-cells IgM-", "B-cells IgM+")
  
  # add column of population labels and column indicating all B cells
  # (also convert population labels to numeric to allow writing .fcs files)
  cbind(d, population = as.numeric(l$population), B_cell = as.numeric(B_cell))
  
}, data, labels, SIMPLIFY = FALSE)



# ----------------------
# Export unmodified data
# ----------------------

# export key for population labels/names

DIR_DATA_OUT <- file.path(DIR_BENCHMARK, "BCR_XL_sim/data")

labels_key <- data.frame(label = 1:length(levels(labels[[1]]$population)), 
                         name = levels(labels[[1]]$population))

file_key <- file.path(DIR_DATA_OUT, "population_names.csv")
write.csv(labels_key, file = file_key, row.names = FALSE)


# export unmodified data

for (i in 1:length(data)) {
  data_i <- data[[i]]
  
  # include spike-in status column so all .fcs files have same shape
  data_out_i <- cbind(data_i, spikein = 0)
  
  file_i <- file.path(DIR_DATA_OUT, "no_sim", basename(files_all)[i])
  write.FCS(flowFrame(data_out_i), file_i)
}



# ---------------------------------------------------------------
# Split each reference sample into two halves: 'base' and 'spike'
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



# -------------------------------------------------------------------------------------
# Replace B cells in 'spike' samples with an equivalent number of B cells from 'BCR-XL'
# (stimulated) condition
# -------------------------------------------------------------------------------------

# note: for some samples, not enough B cells are available; use all available B cells in
# this case

# B cells from 'BCR-XL' (stimulated) condition
data_BCRXL <- data[conditions == "BCR-XL"]

B_cells_BCRXL <- lapply(data_BCRXL, function(d) {
  d[d[, "B_cell"] == 1, ]
})

# number of B cells available
sapply(B_cells_BCRXL, nrow)

# number of B cells needed
n_spike <- sapply(data_spike, function(d) {
  sum(d[, "B_cell"] == 1)
})
n_spike

# select correct number of B cells from 'BCR-XL' (stimulated) condition for each sample

set.seed(123)

B_cells_spike <- mapply(function(b, n) {
  # reduce 'n' if not enough B cells available
  n <- min(n, nrow(b))
  ixs <- sample(seq_len(nrow(b)), n)
  b[ixs, ]
}, B_cells_BCRXL, n_spike, SIMPLIFY = FALSE)

sapply(B_cells_spike, nrow)

# replace B cells in 'spike' samples

data_spike <- mapply(function(d, b) {
  d <- d[d[, "B_cell"] == 0, ]
  d <- rbind(d, b)
  rownames(d) <- NULL
  d
}, data_spike, B_cells_spike, SIMPLIFY = FALSE)



# -----------
# Export data
# -----------

data_export <- c(data_base, data_spike)

conditions_spike <- c(rep("base", length(data_base)), rep("spike", length(data_spike)))

# add column indicating spike-in cells (all B cells in 'spike' samples)
data_export <- mapply(function(d, cnd) {
  spikein <- as.numeric((d[, "B_cell"] == 1) & (cnd == "spike"))
  cbind(d, spikein = spikein)
}, data_export, conditions_spike, SIMPLIFY = FALSE)


filenames <- file.path(DIR_DATA_OUT, "main", paste0("BCR_XL_sim_", patient_IDs, "_", conditions_spike, ".fcs"))

for (i in 1:length(data_export)) {
  write.FCS(flowFrame(data_export[[i]]), filename = filenames[i])
}



# ---------------------------------
# Save timestamp file for Makefiles
# ---------------------------------

file_timestamp <- file.path(DIR_DATA_OUT, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



