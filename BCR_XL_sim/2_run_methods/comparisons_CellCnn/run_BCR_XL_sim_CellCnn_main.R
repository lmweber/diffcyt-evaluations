##########################################################################################
# Script to run methods
# 
# - method: CellCnn
# - data set: BCR-XL-sim
# 
# - main results
# 
# Lukas Weber, November 2017
##########################################################################################


library(flowCore)
library(SummarizedExperiment)


DIR_BENCHMARK <- "../../../../../benchmark_data/BCR_XL_sim/data/main"
DIR_CELLCNN <- "../../../../../CellCnn/CellCnn"
DIR_CELLCNN_FILES <- "../../../../CellCnn_files/BCR_XL_sim/main"
DIR_RDATA <- "../../../../RData/BCR_XL_sim/comparisons_CellCnn"
DIR_SESSION_INFO <- "../../../../session_info/BCR_XL_sim/comparisons_CellCnn"




##############################
# Delete previous output files
##############################

# delete output files from previous CellCnn runs (but leave directory structure intact)

cmd_clean <- paste("find", DIR_CELLCNN_FILES, "-type f -delete")

system(cmd_clean)




###########################
# Load data, pre-processing
###########################

# filenames
files <- list.files(DIR_BENCHMARK, pattern = "\\.fcs$", full.names = TRUE)
files_base <- files[grep("base\\.fcs$", files)]
files_spike <- files[grep("spike\\.fcs$", files)]

# load data
files_load <- c(files_base, files_spike)
files_load

d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)

# sample IDs, group IDs, patient IDs
sample_IDs <- gsub("^BCR_XL_sim_", "", 
                   gsub("\\.fcs$", "", basename(files_load)))
sample_IDs

group_IDs <- factor(gsub("^.*_", "", sample_IDs), levels = c("base", "spike"))
group_IDs

patient_IDs <- factor(gsub("_.*$", "", sample_IDs))
patient_IDs

# check
data.frame(sample_IDs, group_IDs, patient_IDs)

# indices of all marker columns, lineage markers, and functional markers
# (10 surface markers / 14 functional markers; see Bruggner et al. 2014, Table 1)
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)


# --------------
# markers to use
# --------------

cols_to_use <- cols_markers


# -------------------------------
# marker names for CellCnn inputs
# -------------------------------

check <- c()
for (i in 1:length(d_input)) {
  check[i] <- all(colnames(d_input[[i]]) == colnames(d_input[[1]]))
}
all(check)

marker_names <- colnames(d_input[[1]])[cols_to_use]


# --------------
# transform data
# --------------

# not required, since CellCnn automatically transforms data




##################
# CellCnn pipeline
##################

# -----------------------------
# Export .fcs files for CellCnn
# -----------------------------

for (i in 1:length(sample_IDs)) {
  filename <- file.path(DIR_CELLCNN_FILES, "data", basename(files_load[i]))
  write.FCS(d_input[[i]], filename)
}


# --------------------------------
# Generate input files for CellCnn
# --------------------------------

# generate .csv files with input arguments for CellCnn (in required format)

files_in <- basename(files_load)
files_in


# create data frame of sample names and conditions (for CellCnn input .csv file)

label <- group_IDs
label <- as.numeric(label) - 1
label

df_samples <- data.frame(fcs_filename = files_in, label = label)
df_samples

# re-arrange alphabetically (otherwise CellCnn reads input files in incorrect order)
df_samples <- df_samples[order(df_samples$fcs_filename), ]
df_samples


# create data frame of column names (markers) (for CellCnn input .csv file)

df_markers <- t(data.frame(marker_names))
df_markers


# save as .csv files

fn_samples <- file.path(DIR_CELLCNN_FILES, "inputs", "input_samples.csv")
write.csv(df_samples, fn_samples, quote = FALSE, row.names = FALSE)

# need to use 'write.table' to allow removing column names
fn_markers <- file.path(DIR_CELLCNN_FILES, "inputs", "input_markers.csv")
write.table(df_markers, fn_markers, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)


# -----------------------------
# Run CellCnn from command line
# -----------------------------

# for installation instructions and examples see: https://github.com/eiriniar/CellCnn

# note: for Linux systems, need to paste the following additional lines at the beginning of 
# the script 'plotting.py' (before the line 'import matplotlib.pyplot as plt'):
# import matplotlib
# matplotlib.use('agg')

# note: additional advice from authors:
# (1) '--no_arcsinh' argument to disable arcsinh transform
# (2) '--ncell 300' argument to increase size of each training set [this should be at least 
# 'n_cells_min = (avg. no. cells per sample) * (no. samples per condition) / 1000'; 
# if no memory constraints then increase to 5 to 10 times 'n_cells_min']
# (3) '--subset_selection outlier' for extremely rare populations

# command to run CellCnn analysis
cmd <- paste("python", paste0(DIR_CELLCNN, "/cellCnn/run_analysis.py"), 
             paste0("-f ", DIR_CELLCNN_FILES, "/inputs/input_samples.csv"), 
             paste0("-m ", DIR_CELLCNN_FILES, "/inputs/input_markers.csv"), 
             paste0("-i ", DIR_CELLCNN_FILES, "/data/"), 
             paste0("-o ", DIR_CELLCNN_FILES, "/out_CellCnn/"), 
             paste0("--ncell 300 --export_csv"), 
             #"--no_arcsinh",  ## currently not working correctly
             paste("--group_a", "base", "--group_b", "spike"))

# run from command line
runtime_analysis <- system.time(
  system(cmd)
)


# command to export selected cells
cmd <- paste("python", paste0(DIR_CELLCNN, "/cellCnn/run_analysis.py"), 
             paste0("-f ", DIR_CELLCNN_FILES, "/inputs/input_samples.csv"), 
             paste0("-m ", DIR_CELLCNN_FILES, "/inputs/input_markers.csv"), 
             paste0("-i ", DIR_CELLCNN_FILES, "/data/"), 
             paste0("-o ", DIR_CELLCNN_FILES, "/out_CellCnn/"), 
             paste("--group_a", "base", "--group_b", "spike"), 
             "--filter_response_thres 0.3 --load_results --export_selected_cells")

# run from command line
runtime_select <- system.time(
  system(cmd)
)


# runtime
runtime_total <- runtime_analysis[["elapsed"]] + runtime_select[["elapsed"]]
print(runtime_total)

runtime_CellCnn_main <- runtime_total




##############################
# Return results at cell level
##############################

# Note: CellCnn returns continuous 'scores' at the cell level, indicating the
# likelihood of each cell belonging to each detected 'filter' (population). If there
# are multiple detected filters, we sum the scores to give a single total score per
# cell.


# number of cells per sample (including spike-in cells)
n_cells <- sapply(d_input, nrow)

# identify B cells (these contain the true differential signal; from both 'base' and
# 'spike' conditions)
is_B_cell <- unlist(sapply(d_input, function(d) exprs(d)[, "B_cell"]))
stopifnot(length(is_B_cell) == sum(n_cells))


# CellCnn output files

path_out <- file.path(DIR_CELLCNN_FILES, "out_CellCnn/selected_cells")

# if no files exist, CellCnn did not run correctly; return all zeros in this case
if (length(list.files(path_out)) == 0) {
  filter_continuous <- list(rep(0, sum(n_cells)))
}

files_cells <- paste0(path_out, "/", gsub("\\.fcs$", "", basename(files_load)), "_selected_cells.csv")

# get cells in selected filters
filter_continuous <- vector("list", length(files_cells))

for (f in 1:length(files_cells)) {
  d <- try(read.csv(files_cells[f]), silent = TRUE)
  
  if (!(class(d) == "try-error")) {
    # note: if there are multiple filters for one sample, combine them using 'rowSums'
    # (filters are stored in odd-numbered columns of the .csv file)
    ix <- seq(1, ncol(d), by = 2)
    filt_sum <- rowSums(d[, ix, drop = FALSE])
    filter_continuous[[f]] <- filt_sum
    
  } else {
    # if .csv file is empty, fill with zeros instead
    filter_continuous[[f]] <- rep(0, n_cells[[f]])
  }
}

filter_continuous <- unlist(filter_continuous)


# set up data frame with results and true B-cell status at cell level

scores <- filter_continuous

# replace any NAs to ensure same cells are returned for all methods
scores[is.na(scores)] <- 0

stopifnot(length(scores) == sum(n_cells), 
          length(scores) == length(is_B_cell))

res <- data.frame(scores = scores, 
                  B_cell = is_B_cell)

# store results
out_CellCnn_main <- res




#####################
# Save output objects
#####################

save(out_CellCnn_main, runtime_CellCnn_main, 
     file = file.path(DIR_RDATA, "outputs_BCR_XL_sim_CellCnn_main.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_BCR_XL_sim_CellCnn_main.txt"))
sessionInfo()
sink()



