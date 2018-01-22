##########################################################################################
# Script to prepare benchmark data set 'AML-sim'
#
# The 'AML-sim' data set is constructed by computationally 'spiking in' small percentages
# of AML (acute myeloid leukemia) blast cells into samples of healthy BMMCs (bone marrow
# mononuclear cells). This simulates the phenotype of minimal residual disease (MRD) in
# AML patients. Raw data is sourced from Levine et al. (2015) (PhenoGraph paper). The data
# generation strategy is modified from Arvaniti et al. (2017) (CellCnn paper), who
# generated a similar benchmark data set for their evaluations.
#
# Raw data downloaded from Cytobank:
# - all cells (also contains gating scheme for CD34+ CD45 mid cells, i.e. blasts):
# https://community.cytobank.org/cytobank/experiments/46098/illustrations/121588
# - blasts (repository cloned from the one for 'all cells' above, using the gating scheme
# for CD34+ CD45 mid cells; this allows .fcs files for the subset to be exported):
# https://community.cytobank.org/cytobank/experiments/63534/illustrations/125318
#
# Notes:
# - Gating plots for blasts are also shown in Levine et al. (2015), Supplemental Data S3B.
# - Individuals SJ1, SJ2, and SJ3 each contain two replicates; these are in separate .fcs
# files. The original Cytobank repository combines the two replicates for each individual
# (see 'Individuals' dimension setup); so we should use combined cells from both .fcs
# files in downstream analysis. However, when I tried this for SJ1, the percentage of
# blasts did not match to the published numbers (shown in Levine et al. 2015, Supplemental
# Data S3B); so we have not used these samples in the final analysis.
# - Arvaniti et al. (2017) (CellCnn paper) classified patients SJ10, SJ12, SJ13 as CN
# (cytogenetically normal), and SJ1, SJ2, SJ3, SJ4, SJ5 as CBF (core-binding factor
# translocation); we re-use these classifications here.
# - Sample names and filenames in the raw data are shuffled (e.g. file H3 is actually
# sample H1). The matching scheme can be seen in the 'Individuals' setup in Cytobank, or
# in the downloaded .tsv files 'experiment_46098_annotations.tsv' and
# 'experiment_63534_annotations.tsv'.
#
# Lukas Weber, January 2018
##########################################################################################


# modified to create 'null' simulations: no true spike-in cells


library(flowCore)



# ---------
# Filenames
# ---------

DIR_BENCHMARK <- "../../../../../benchmark_data"

DIR_RAW_DATA_ALL <- file.path(DIR_BENCHMARK, "AML_sim/raw_data/all_cells/experiment_46098_files")
DIR_RAW_DATA_BLASTS <- file.path(DIR_BENCHMARK, "AML_sim/raw_data/CD34_CD45mid_cells/experiment_63534_files")

DIR_DATA_OUT <- file.path(DIR_BENCHMARK, "AML_sim/data/null_simulations")

files_all <- list.files(DIR_RAW_DATA_ALL, pattern = "\\.fcs$", full.names = TRUE)
files_healthy <- files_all[grep("H[0-9]+", files_all)]

files_blasts <- list.files(DIR_RAW_DATA_BLASTS, pattern = "\\.fcs$", full.names = TRUE)

# metadata spreadsheets to match shuffled sample names
file_match_samples_all <- file.path(DIR_RAW_DATA_ALL, "experiment_46098_annotations.tsv")
file_match_samples_blasts <- file.path(DIR_RAW_DATA_BLASTS, "experiment_63534_annotations.tsv")



# -----------------------------------------------
# Load data for healthy samples H1-H5 (all cells)
# -----------------------------------------------

data_healthy <- lapply(files_healthy, function(f) exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)))

# show shuffled sample names
tbl_match_healthy <- read.delim(file_match_samples_all)
tbl_match_healthy_sub <- tbl_match_healthy[grep("H[0-9]+", tbl_match_healthy[, "FCS.Filename"]), ]
tbl_match_healthy_sub[, c("FCS.Filename", "Individuals")]

# match correct sample names; store as names of list items
names(data_healthy) <- tbl_match_healthy_sub[, "Individuals"]

length(data_healthy)
sapply(data_healthy, dim)



# ---------------------
# Randomized replicates
# ---------------------

# use different random seed for each replicate

n_replicates <- 3

for (r in 1:n_replicates) {
  
  
  
  # ---------------------
  # Split healthy samples
  # ---------------------
  
  # Null simulation: split each healthy sample (H1-H5) into 2 equal parts.
  
  data_healthy_null_1 <- data_healthy_null_2 <- vector("list", length(data_healthy))
  names(data_healthy_null_1) <- names(data_healthy_null_2) <- names(data_healthy)
  
  # modified random seed for each replicate
  seed <- 10000 + 100 * r
  
  for (i in 1:length(data_healthy)) {
    data_i <- data_healthy[[i]]
    
    # modified random seed for each replicate
    set.seed(seed + i)
    
    # null simulation: 2 equal parts
    n <- round(nrow(data_i) / 2)
    
    ix_null_1 <- sample(1:nrow(data_i), n)
    ix_null_2 <- setdiff(1:nrow(data_i), ix_null_1)
    
    data_healthy_null_1[[i]] <- data_i[ix_null_1, ]
    data_healthy_null_2[[i]] <- data_i[ix_null_2, ]
  }
  
  sapply(data_healthy_null_1, dim)
  sapply(data_healthy_null_2, dim)
  
  
  
  # Export 'null_1' and 'null_2' for all samples (H1-H5)
  
  # save .fcs files
  for (i in 1:length(data_healthy_null_1)) {
    data_null_1_i <- data_healthy_null_1[[i]]
    data_null_2_i <- data_healthy_null_2[[i]]
    
    nm_i <- names(data_healthy_null_1)[i]
    
    # include spike-in status column so all .fcs files have same shape
    data_null_1_out_i <- cbind(data_null_1_i, spikein = 0)
    data_null_2_out_i <- cbind(data_null_2_i, spikein = 0)
    
    filename_null_1 <- file.path(DIR_DATA_OUT, paste0("seed", r), "null1", paste0("AML_sim_", nm_i, "_seed", r, "_null1.fcs"))
    filename_null_2 <- file.path(DIR_DATA_OUT, paste0("seed", r), "null2", paste0("AML_sim_", nm_i, "_seed", r, "_null2.fcs"))
    
    write.FCS(flowFrame(data_null_1_out_i), filename_null_1)
    write.FCS(flowFrame(data_null_2_out_i), filename_null_2)
  }
  
}



# ---------------------------------
# Save timestamp file for Makefiles
# ---------------------------------

file_timestamp <- file.path(DIR_DATA_OUT, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



