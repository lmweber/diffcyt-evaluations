##########################################################################################
# Script to run methods
# 
# - method: Citrus
# - data set: BCR-XL-sim
# 
# - main results
# 
# Lukas Weber, May 2018
##########################################################################################


library(flowCore)
library(citrus)


DIR_BENCHMARK <- "../../../../../benchmark_data/BCR_XL_sim/data/main"
DIR_CITRUS_FILES <- "../../../../Citrus_files/BCR_XL_sim/main"
DIR_RDATA <- "../../../../RData/BCR_XL_sim/comparisons_Citrus"
DIR_SESSION_INFO <- "../../../../session_info/BCR_XL_sim/comparisons_Citrus"




##############################
# Delete previous output files
##############################

# delete output files from previous Citrus runs (but leave directory structure intact)

cmd_clean <- paste("find", DIR_CITRUS_FILES, "-type f -delete")

system(cmd_clean)




###########################
# Load data, pre-processing
###########################

# filenames

files <- list.files(DIR_BENCHMARK, pattern = "\\.fcs$", full.names = TRUE)
files_base <- files[grep("base\\.fcs$", files)]
files_spike <- files[grep("spike\\.fcs$", files)]

files_load <- c(files_base, files_spike)
files_load

# load data

d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)

# sample information

sample_id <- gsub("^BCR_XL_sim_", "", 
                  gsub("\\.fcs$", "", basename(files_load)))
sample_id

group_id <- factor(gsub("^.*_", "", sample_id), levels = c("base", "spike"))
group_id

patient_id <- factor(gsub("_.*$", "", sample_id))
patient_id

experiment_info <- data.frame(group_id, patient_id, sample_id)
experiment_info

# marker information

# indices of all marker columns, lineage markers, and functional markers
# (10 surface markers / 14 functional markers; see Bruggner et al. 2014, Table 1)
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)

marker_name <- colnames(d_input[[1]])
marker_name <- gsub("\\(.*$", "", marker_name)

marker_class <- rep("none", length(marker_name))
marker_class[cols_lineage] <- "type"
marker_class[cols_func] <- "state"
marker_class <- factor(marker_class, levels = c("type", "state", "none"))

# exclude CD45 from clustering
marker_class[marker_name == "CD45"] <- "none"

marker_info <- data.frame(marker_name, marker_class)
marker_info




######################################
# Additional pre-processing for Citrus
######################################

# --------------
# markers to use
# --------------

# excluding CD45
cols_to_use_clustering <- colnames(d_input[[1]])[cols_lineage]
cols_to_use_clustering <- cols_to_use_clustering[-grep("CD45", cols_to_use_clustering)]

cols_to_use_medians <- colnames(d_input[[1]])[cols_func]


# --------------
# transform data
# --------------

# 'arcsinh' transform with 'cofactor' = 5 (see Bendall et al. 2011, Supp. Fig. S2)

cofactor <- 5

d_input <- lapply(d_input, function(d) {
  e <- exprs(d)
  e[, cols_markers] <- asinh(e[, cols_markers] / cofactor)
  flowFrame(e)
})




#################
# Citrus pipeline
#################

# using modified code from auto-generated file 'runCitrus.R'


# ----------------------------------------
# Export transformed .fcs files for Citrus
# ----------------------------------------

for (i in 1:length(sample_id)) {
  filename <- file.path(DIR_CITRUS_FILES, "data_transformed", 
                        gsub("\\.fcs$", "_transf.fcs", basename(files_load[i])))
  write.FCS(d_input[[i]], filename)
}


# ------------------------
# Define inputs for Citrus
# ------------------------

# model type
family <- "classification"
modelTypes <- "glmnet"
nFolds <- 1

# feature type (use 'medians' for functional markers)
featureType <- "medians"

# define clustering and functional markers
clusteringColumns <- cols_to_use_clustering
medianColumns <- cols_to_use_medians

# number of cells and minimum cluster size
fileSampleSize <- 5000
minimumClusterSizePercent <- 0.01  # 1%

# transformation: not required since already done above
transformColumns <- NULL
transformCofactor <- NULL
scaleColumns <- NULL

# experimental design
labels <- group_id

# directories
dataDirectory <- file.path(DIR_CITRUS_FILES, "data_transformed")
outputDirectory <- file.path(DIR_CITRUS_FILES, "citrusOutput")

# files
fileList <- data.frame(defaultCondition = gsub("\\.fcs$", "_transf.fcs", basename(files_load)))

# number of processor threads
n_cores <- 1


# run Citrus

runtime_Citrus <- system.time({
  
  Rclusterpp.setThreads(n_cores)
  
  set.seed(12345)
  
  results <- citrus.full(
    fileList = fileList, 
    labels = labels, 
    clusteringColumns = clusteringColumns, 
    medianColumns = medianColumns, 
    dataDirectory = dataDirectory, 
    outputDirectory = outputDirectory, 
    family = family, 
    modelTypes = modelTypes, 
    nFolds = nFolds, 
    fileSampleSize = fileSampleSize, 
    featureType = featureType, 
    minimumClusterSizePercent = minimumClusterSizePercent, 
    transformColumns = transformColumns, 
    transformCofactor = transformCofactor, 
    scaleColumns = scaleColumns
  )
  
})

# Citrus plots
plot(results, outputDirectory)

# runtime
runtime_total <- runtime_Citrus[["elapsed"]]
print(runtime_total)

runtime_Citrus_main <- runtime_total




##############################
# Return results at cell level
##############################

# get differential clusters, match cells to clusters, and save results at cell level

# note: Citrus does not give any continuous-valued scores, e.g. p-values or q-values,
# so it is not possible to rank the selected clusters


# number of cells per sample (including spike-in cells)
n_cells <- sapply(d_input, nrow)

# identify B cells (these contain the true differential signal; from both 'base' and
# 'spike' conditions)
is_B_cell <- unlist(sapply(d_input, function(d) exprs(d)[, "B_cell"]))
stopifnot(length(is_B_cell) == sum(n_cells))


# differential clusters

clusters <- as.numeric(results$conditionRegressionResults$defaultCondition$glmnet$differentialFeatures[["cv.min"]][["clusters"]])


# match clusters to cells

res_cells <- rep(NA, sum(n_cells))
files_rep <- rep(seq_along(sample_id), n_cells)

for (i in seq_along(clusters)) {
  
  clus <- clusters[i]
  
  # identify cells in cluster
  ix_clus <- results$citrus.foldClustering$allClustering$clusterMembership[[clus]]
  
  data <- results$citrus.combinedFCSSet$data
  ix_events <- data[, "fileEventNumber"]
  ix_files <- data[, "fileId"]
  
  cells_subsampled <- rep(NA, nrow(data))
  cells_subsampled[ix_clus] <- 1  # 1 = cell is in differential cluster
  
  # match to indices in original data; taking into account subsampling
  for (z in seq_along(group_id)) {
    ix_events_z <- ix_events[ix_files == z]
    res_cells[files_rep == z][ix_events_z] <- cells_subsampled[ix_files == z]
  }
  
  # cluster marker expression values (if required)
  #results$citrus.combinedFCSSet$data[results$citrus.foldClustering$allClustering$clusterMembership[[clus]], clusteringColumns]
}

# set all NA values to 0 to allow evaluation
# (note: Citrus results are binary; 1 = cell selected, 0 = cell not selected)
res_cells[is.na(res_cells)] <- 0


# set up data frame with results and true B-cell status at cell level

scores <- res_cells

# replace any NAs to ensure same cells are returned for all methods
scores[is.na(scores)] <- 0

stopifnot(length(scores) == sum(n_cells), 
          length(scores) == length(is_B_cell))

res <- data.frame(scores = scores, 
                  B_cell = is_B_cell)

# store results
out_Citrus_main <- res




#####################
# Save output objects
#####################

save(out_Citrus_main, runtime_Citrus_main, 
     file = file.path(DIR_RDATA, "outputs_BCR_XL_sim_Citrus_main.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_BCR_XL_sim_Citrus_main.txt"))
sessionInfo()
sink()



