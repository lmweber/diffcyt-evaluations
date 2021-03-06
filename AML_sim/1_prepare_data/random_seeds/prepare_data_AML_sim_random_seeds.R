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


# modified to generate randomized benchmark data sets using different random seeds


library(flowCore)



# ---------
# Filenames
# ---------

DIR_BENCHMARK <- "../../../../../benchmark_data"

DIR_RAW_DATA_ALL <- file.path(DIR_BENCHMARK, "AML_sim/raw_data/all_cells/experiment_46098_files")
DIR_RAW_DATA_BLASTS <- file.path(DIR_BENCHMARK, "AML_sim/raw_data/CD34_CD45mid_cells/experiment_63534_files")

DIR_DATA_OUT <- file.path(DIR_BENCHMARK, "AML_sim/data/random_seeds")

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



# -------------------------------------------------
# Load data for healthy samples H1-H5 (blast cells)
# -------------------------------------------------

# note sample names and filenames are shuffled
tbl_match_blasts <- read.delim(file_match_samples_blasts)
tbl_match_blasts[grep("H[0-9]+", tbl_match_blasts[, "FCS.Filename"]), c("FCS.Filename", "Individuals")]

files_healthy_blasts <- files_blasts[1:5]

data_healthy_blasts <- lapply(files_healthy_blasts, function(f) exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)))

names(data_healthy_blasts) <- names(data_healthy)

# check numbers of cells
sapply(data_healthy_blasts, dim)



# ----------------------------------------------------------------------------
# Load data for AML patients: blast cells (CN: patient SJ10; CBF: patient SJ4)
# ----------------------------------------------------------------------------

# note sample names and filenames are shuffled
tbl_match_blasts <- read.delim(file_match_samples_blasts)
tbl_match_blasts[-grep("H[0-9]+", tbl_match_blasts[, "FCS.Filename"]), c("FCS.Filename", "Individuals")]

file_SJ10 <- files_blasts[6]
file_SJ10  # note: sample 'SJ10' has filename 'SJ11'

file_SJ4 <- files_blasts[20]
file_SJ4  # note: sample 'SJ4' has filename 'SJ5'

# load data for SJ10 (CN)
data_SJ10 <- exprs(read.FCS(file_SJ10, transformation = FALSE, truncate_max_range = FALSE))

# load data for SJ4 (CBF)
data_SJ4 <- exprs(read.FCS(file_SJ4, transformation = FALSE, truncate_max_range = FALSE))


# check column names match for all samples (healthy and blasts)
for (i in 1:5) {
  print(all.equal(colnames(data_healthy[[i]]), colnames(data_SJ4)))
}
all.equal(colnames(data_SJ10), colnames(data_SJ4))



# ---------------------
# Randomized replicates
# ---------------------

# use different random seed for each replicate

n_replicates <- 3

for (r in 1:n_replicates) {
  
  
  
  # ---------------------
  # Split healthy samples
  # ---------------------
  
  # Split each healthy sample (H1-H5) into 3 equal parts. One part will be used as the
  # healthy sample, and the other parts will each have spike-in cells added (for conditions
  # CN and CBF).
  
  data_healthy_base <- data_healthy_CN <- data_healthy_CBF <- 
    vector("list", length(data_healthy))
  names(data_healthy_base) <- names(data_healthy_CN) <- names(data_healthy_CBF) <- 
    names(data_healthy)
  
  # modified random seed for each replicate
  seed <- 1000 + 100 * r
  
  for (i in 1:length(data_healthy)) {
    data_i <- data_healthy[[i]]
    
    # modified random seed for each replicate
    set.seed(seed + i)
    
    n <- round(nrow(data_i) / 3)
    
    ix_base <- sample(1:nrow(data_i), n)
    ix_CN <- sample((1:nrow(data_i))[-ix_base], n)
    ix_CBF <- setdiff(1:nrow(data_i), c(ix_base, ix_CN))
    
    data_healthy_base[[i]] <- data_i[ix_base, ]
    data_healthy_CN[[i]] <- data_i[ix_CN, ]
    data_healthy_CBF[[i]] <- data_i[ix_CBF, ]
  }
  
  sapply(data_healthy_base, dim)
  sapply(data_healthy_CN, dim)
  sapply(data_healthy_CBF, dim)
  
  
  
  # Export healthy samples (H1-H5)
  
  # save .fcs files
  for (i in 1:length(data_healthy_base)) {
    data_i <- data_healthy_base[[i]]
    nm_i <- names(data_healthy_base)[i]
    
    # include spike-in status column so all .fcs files have same shape
    data_out_i <- cbind(data_i, spikein = 0)
    
    filename <- file.path(DIR_DATA_OUT, paste0("seed", r), "healthy", 
                          paste0("AML_sim_healthy_", nm_i, "_randomseed", r, ".fcs"))
    write.FCS(flowFrame(data_out_i), filename)
  }
  
  
  
  # -------------------------
  # Create spike-in data sets
  # -------------------------
  
  # AML blast cells are subsampled at various thresholds (5%, 1%, 0.1%) of the total number
  # of healthy cells for each sample, and combined with the healthy cells to create the
  # spike-in data sets.
  
  thresholds <- c(0.05, 0.01, 0.001)  # 5%, 1%, 0.1%
  
  
  # condition CN (patient SJ10)
  
  data_blasts_AML <- data_SJ10
  cnd <- "CN"
  
  # modified random seed for each replicate
  seed <- 2000 + 100 * r
  
  for (i in 1:length(data_healthy_CN)) {
    data_i <- data_healthy_CN[[i]]
    nm_i <- names(data_healthy_CN)[i]
    
    for (th in thresholds) {
      
      # modified random seed for each replicate
      set.seed(seed + 10 * th + i)
      
      n_spikein <- ceiling(th * nrow(data_i))
      is_spikein <- c(rep(0, nrow(data_i)), rep(1, n_spikein))
      
      cat("n =", n_spikein, "\n")
      
      # subsample blasts
      spikein_i <- data_blasts_AML[sample(1:nrow(data_blasts_AML), n_spikein), , drop = FALSE]
      
      data_out_i <- rbind(data_i, spikein_i)
      data_out_i <- cbind(data_out_i, spikein = is_spikein)
      
      filename <- file.path(DIR_DATA_OUT, paste0("seed", r), cnd, 
                            paste0("AML_sim_", cnd, "_", nm_i, "_", th * 100, "pc", 
                                   "_randomseed", r, ".fcs"))
      write.FCS(flowFrame(data_out_i), filename)
    }
  }
  
  
  
  # condition CBF (patient SJ4)
  
  data_blasts_AML <- data_SJ4
  cnd <- "CBF"
  
  # modified random seed for each replicate
  seed <- 3000 + 100 * r
  
  for (i in 1:length(data_healthy_CBF)) {
    data_i <- data_healthy_CBF[[i]]
    nm_i <- names(data_healthy_CBF)[i]
    
    for (th in thresholds) {
      
      # modified random seed for each replicate
      set.seed(seed + 10 * th + i)
      
      n_spikein <- ceiling(th * nrow(data_i))
      is_spikein <- c(rep(0, nrow(data_i)), rep(1, n_spikein))
      
      cat("n =", n_spikein, "\n")
      
      # subsample blasts
      spikein_i <- data_blasts_AML[sample(1:nrow(data_blasts_AML), n_spikein), , drop = FALSE]
      
      data_out_i <- rbind(data_i, spikein_i)
      data_out_i <- cbind(data_out_i, spikein = is_spikein)
      
      filename <- file.path(DIR_DATA_OUT, paste0("seed", r), cnd, 
                            paste0("AML_sim_", cnd, "_", nm_i, "_", th * 100, "pc", 
                                   "_randomseed", r, ".fcs"))
      write.FCS(flowFrame(data_out_i), filename)
    }
  }
  
}



# ---------------------------------
# Save timestamp file for Makefiles
# ---------------------------------

file_timestamp <- file.path(DIR_DATA_OUT, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



