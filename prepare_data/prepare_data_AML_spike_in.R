##########################################################################################
# Script to prepare benchmark data set 'AML-spike-in'
#
# The 'AML-spike-in' data set is constructed by 'spiking in' (computationally, i.e. in 
# silico) small percentages of AML (acute myeloid leukemia) blast cells into samples of 
# healthy BMMCs (bone marrow mononuclear cells). This simulates the phenotype of minimal 
# residual disease (MRD) in AML patients. Raw data is sourced from Levine et al. (2015) 
# (PhenoGraph paper). A similar strategy was used by Arvaniti et al. (2017) (CellCnn 
# paper) to generate one of their benchmark data sets, although using different settings.
#
# Raw data downloaded from Cytobank:
# - all cells (also contains gating scheme for CD34+ CD45 mid cells, i.e. blasts):
# https://community.cytobank.org/cytobank/experiments/46098/illustrations/121588
# - blasts (this repository was cloned from the one for all cells above, using the gating 
# scheme for CD34+ CD45mid cells; this allows FCS files for the subset to be exported): 
# https://community.cytobank.org/cytobank/experiments/63534/illustrations/125318
#
# Notes:
# - Gating plots for blasts are also shown in Levine et al. (2015), Supplemental Data S3B.
# - Individuals SJ1, SJ2, and SJ3 each contain two replicates; these are in separate FCS 
# files. The original Cytobank repository combines the two replicates for each individual 
# (see 'Individuals' dimension setup). So should use combined cells from both FCS files in
# downstream analysis. However, when I tried this for SJ1, the percentage of blasts did 
# not match to the published numbers (shown in Levine et al. 2015, Supplemental Data S3B);
# so we have not used these samples in the final analysis.
# - Arvaniti et al. (2017) (CellCnn paper) classified patients SJ10, SJ12, SJ13 as CN 
# (cytogenetically normal), and SJ1, SJ2, SJ3, SJ4, SJ5 as CBF (core-binding factor 
# translocation); we re-use these classifications here.
# - Sample names and file names in the raw data are shuffled! (e.g. file H3 is actually
# sample H1). The matching scheme can be seen in the 'Individuals' setup in Cytobank, or
# in the downloaded .tsv files 'experiment_46098_annotations.tsv' and
# 'experiment_63534_annotations.tsv'.
#
# Lukas Weber, May 2017
##########################################################################################


library(flowCore)



# ---------
# Filenames
# ---------

DIR_RAW_DATA_ALL <- "../../../benchmark_data/AML_spike_in/raw_data/all_cells/experiment_46098_files"
DIR_RAW_DATA_BLASTS <- "../../../benchmark_data/AML_spike_in/raw_data/CD34_CD45mid_cells/experiment_63534_files"

DIR_DATA <- "../../../benchmark_data/AML_spike_in/data"

files_all <- list.files(DIR_RAW_DATA_ALL, pattern = "\\.fcs$", full.names = TRUE)
files_healthy <- files_all[grep("H[0-9]+", files_all)]

files_blasts <- list.files(DIR_RAW_DATA_BLASTS, pattern = "\\.fcs$", full.names = TRUE)

# to match shuffled sample names
file_match_samples_all <- file.path(DIR_RAW_DATA_ALL, "experiment_46098_annotations.tsv")
file_match_samples_blasts <- file.path(DIR_RAW_DATA_BLASTS, "experiment_63534_annotations.tsv")



# ------------------------------------
# Load data from healthy samples H1-H5
# ------------------------------------

data_healthy <- lapply(files_healthy, function(f) exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)))

# show shuffled sample names
tbl_match_healthy <- read.delim(file_match_samples_all)
tbl_match_healthy_sub <- tbl_match_healthy[grep("H[0-9]+", tbl_match_healthy[, "FCS.Filename"]), ]
tbl_match_healthy_sub[, c("FCS.Filename", "Individuals")]

# match correct sample names; store as names of list items
names(data_healthy) <- tbl_match_healthy_sub[, "Individuals"]

length(data_healthy)



# ----------------------------------------------------------------
# Load data from AML patients (CN: patient SJ10; CBF: patient SJ4)
# ----------------------------------------------------------------

# note sample names and filenames are shuffled
tbl_match_blasts <- read.delim(file_match_samples_blasts)
tbl_match_blasts[-grep("H[0-9]+", tbl_match_blasts[, "FCS.Filename"]), c("FCS.Filename", "Individuals")]

file_SJ10 <- files_blasts[6]
file_SJ10  # note: has filename SJ11

file_SJ4 <- files_blasts[20]
file_SJ4  # note: has filename SJ5

# load data for SJ10 (CN)
data_SJ10 <- exprs(read.FCS(file_SJ10, transformation = FALSE, truncate_max_range = FALSE))

# load data for SJ4 (CBF)
data_SJ4 <- exprs(read.FCS(file_SJ4, transformation = FALSE, truncate_max_range = FALSE))


# check column names match for all samples (healthy and blasts)
for (i in 1:5) {
  print(all.equal(colnames(data_healthy[[i]]), colnames(data_SJ4)))
}
all.equal(colnames(data_SJ10), colnames(data_SJ4))



# ----------------------
# Check numbers of cells
# ----------------------

# check numbers of cells

# healthy
sapply(data_healthy, dim)

# SJ10: should be 80.7% of total (Levine et al. 2015, Supplemental Data S3B)
dim(data_SJ10)
dim(exprs(read.FCS(paste0("../../../benchmark_data/AML_spike_in/raw_data/all_cells/experiment_46098_files/", 
                          "SJ11d_min3_s0.10_m10_debar1_NoDrug_Basal1_Viable_NoDrug_Basal1_SJ11d.fcs"))))
# check
37615 / 46601  # 80.7%

# SJ4: should be 55.2% of total (Levine et al. 2015, Supplemental Data S3B)
dim(data_SJ4)
dim(exprs(read.FCS(paste0("../../../benchmark_data/AML_spike_in/raw_data/all_cells/experiment_46098_files/", 
                          "SJ5d_min5_s0.15_m10_debar1_NoDrug_Basal1_Viable_NoDrug_Basal1_SJ5d.fcs"))))
# check
14520 / 26321  # 55.2%



# -------------------------------------------
# Up/downsample and create spike-in data sets
# -------------------------------------------

# Healthy samples (H1-H5) are upsampled 10x to create simulated data sets with enough
# cells to allow rare populations to be spiked in at various thresholds. For the smallest
# data set (H1) with the lowest threshold (0.01%), this gives ~16 spike-in cells.

n <- 10

set.seed(100)

data_healthy_up <- lapply(data_healthy, function(d) {
  d[sample(1:nrow(d), n * nrow(d), replace = TRUE), ]
})

# save as new .fcs files
for (i in 1:length(data_healthy_up)) {
  data_i <- data_healthy_up[[i]]
  nm_i <- names(data_healthy_up)[i]
  filename <- file.path(DIR_DATA, "healthy", paste0("AML_spike_in_healthy_", nm_i, ".fcs"))
  write.FCS(flowFrame(data_i), filename)
}


# Blast cells are subsampled at various thresholds (1%, 0.1%, 0.01%) of the number of 
# healthy cells for each sample, and combined with the healthy cells to create the
# spike-in data sets.

thresholds <- c(0.01, 0.001, 0.0001)  # 1%, 0.1%, 0.01%


# condition CN (patient SJ10)

data_blasts <- data_SJ10
cnd <- "CN"

set.seed(101)

for (i in 1:length(data_healthy_up)) {
  data_i <- data_healthy_up[[i]]
  nm_i <- names(data_healthy_up)[i]
  
  for (th in thresholds) {
    n_spikein <- ceiling(th * nrow(data_i))
    cat("n =", n_spikein, "\n")
    
    # subsample blasts
    spikein_i <- data_blasts[sample(1:nrow(data_blasts), n_spikein), ]
    data_out_i <- rbind(data_i, spikein_i)
    
    filename <- file.path(DIR_DATA, cnd, paste0("AML_spike_in_", cnd, "_", nm_i, "_", th * 100, "pc.fcs"))
    write.FCS(flowFrame(data_out_i), filename)
  }
}


# condition CBF (patient SJ4)

data_blasts <- data_SJ4
cnd <- "CBF"

set.seed(102)

for (i in 1:length(data_healthy_up)) {
  data_i <- data_healthy_up[[i]]
  nm_i <- names(data_healthy_up)[i]
  
  for (th in thresholds) {
    n_spikein <- ceiling(th * nrow(data_i))
    cat("n =", n_spikein, "\n")
    
    # subsample blasts
    spikein_i <- data_blasts[sample(1:nrow(data_blasts), n_spikein), ]
    data_out_i <- rbind(data_i, spikein_i)
    
    filename <- file.path(DIR_DATA, cnd, paste0("AML_spike_in_", cnd, "_", nm_i, "_", th * 100, "pc.fcs"))
    write.FCS(flowFrame(data_out_i), filename)
  }
}


