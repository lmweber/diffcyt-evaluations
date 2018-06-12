##########################################################################################
# Script to save lists of benchmark data files (including paths) used for 'diffcyt'
# evaluations
# 
# Lukas Weber, June 2018
##########################################################################################


DIR_CURRENT <- getwd()

DIR_BENCHMARK <- "../../../benchmark_data"
DIR_OUT <- "../../../benchmark_data/benchmark_data_paths"

# set working directory so paths are shown as relative paths within benchmark data directory
setwd(DIR_BENCHMARK)


# -------------------------
# List benchmark data files
# -------------------------

# list files separately for each benchmark dataset (note: includes .fcs, .csv, and .xlsx
# files, depending on the dataset)


# AML-sim

paths_AML_sim <- c(
  "AML_sim/data/less_distinct", 
  "AML_sim/data/main", 
  "AML_sim/data/null_simulations", 
  "AML_sim/data/random_seeds"
)

files_AML_sim <- list.files(paths_AML_sim, pattern = "\\.fcs$", recursive = TRUE, full.names = TRUE)


# BCR-XL-sim

paths_BCR_XL_sim_fcs <- c(
  "BCR_XL_sim/data/less_distinct", 
  "BCR_XL_sim/data/main", 
  "BCR_XL_sim/data/null_simulations", 
  "BCR_XL_sim/data/random_seeds"
)

paths_BCR_XL_sim_csv <- c(
  "BCR_XL_sim/data/population_names", 
  "BCR_XL_sim/population_IDs"
)

files_BCR_XL_sim_fcs <- list.files(paths_BCR_XL_sim_fcs, pattern = c("\\.fcs$"), recursive = TRUE, full.names = TRUE)
files_BCR_XL_sim_csv <- list.files(paths_BCR_XL_sim_csv, pattern = c("\\.csv$"), recursive = TRUE, full.names = TRUE)

files_BCR_XL_sim <- c(files_BCR_XL_sim_fcs, files_BCR_XL_sim_csv)


# Anti-PD-1

paths_Anti_PD_1 <- c(
  "Anti_PD_1/CK_2016-06-23_03all", 
  "Anti_PD_1/CK_2016-06-29_03all3", 
  "Anti_PD_1/CK_metadata", 
  "Anti_PD_1/CK_panels"
)

files_Anti_PD_1 <- list.files(paths_Anti_PD_1, pattern = c("\\.fcs$|\\.xlsx$"), recursive = TRUE, full.names = TRUE)


# BCR-XL

paths_BCR_XL <- c(
  "BCR_XL/Bodenmiller_BCR_XL_fcs_files", 
  "BCR_XL/Bodenmiller_BCR_XL_population_IDs"
)

files_BCR_XL <- list.files(paths_BCR_XL, pattern = c("\\.fcs$|\\.csv$"), recursive = TRUE, full.names = TRUE)


# ------------------
# Save in .csv files
# ------------------

# save lists of benchmark data files in .csv files

# return to original working directory
setwd(DIR_CURRENT)

write.table(files_AML_sim,    file.path(DIR_OUT, "benchmark_data_paths_AML_sim.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
write.table(files_BCR_XL_sim, file.path(DIR_OUT, "benchmark_data_paths_BCR_XL_sim.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
write.table(files_Anti_PD_1,  file.path(DIR_OUT, "benchmark_data_paths_Anti_PD_1.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
write.table(files_BCR_XL,     file.path(DIR_OUT, "benchmark_data_paths_BCR_XL.csv"), sep = ",", row.names = FALSE, col.names = FALSE)


