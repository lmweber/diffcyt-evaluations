##########################################################################################
# Script to run methods
# 
# - method: Citrus
# - data set: AML-sim
# 
# - main results
# 
# Lukas Weber, May 2018
##########################################################################################


library(flowCore)
library(citrus)


DIR_BENCHMARK <- "../../../../../benchmark_data/AML_sim/data/main"
DIR_CITRUS_FILES <- "../../../../Citrus_files/AML_sim/main"
DIR_RDATA <- "../../../../RData/AML_sim/comparisons_Citrus"
DIR_SESSION_INFO <- "../../../../session_info/AML_sim/comparisons_Citrus"




##############################
# Delete previous output files
##############################

# delete output files from previous Citrus runs (but leave directory structure intact)

cmd_clean <- paste("find", DIR_CITRUS_FILES, "-type f -delete")

system(cmd_clean)




################################
# Loop to run for each threshold
################################

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# condition names
cond_names <- c("CN", "CBF")

# lists to store objects and runtime
out_Citrus_main <- runtime_Citrus_main <- vector("list", length(thresholds))
names(out_Citrus_main) <- names(runtime_Citrus_main) <- thresholds




for (th in 1:length(thresholds)) {
  
  ###########################
  # Load data, pre-processing
  ###########################
  
  # filenames
  
  files_healthy <- list.files(file.path(DIR_BENCHMARK, "healthy"), 
                              pattern = "\\.fcs$", full.names = TRUE)
  files_CN <- list.files(file.path(DIR_BENCHMARK, "CN"), 
                         pattern = paste0("_", thresholds[th], "\\.fcs$"), full.names = TRUE)
  files_CBF <- list.files(file.path(DIR_BENCHMARK, "CBF"), 
                          pattern = paste0("_", thresholds[th], "\\.fcs$"), full.names = TRUE)
  
  files_load <- c(files_healthy, files_CN, files_CBF)
  files_load
  
  # load data
  
  d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)
  
  # sample IDs, group IDs, patient IDs
  sample_id <- gsub("(_[0-9]+pc$)|(_0\\.[0-9]+pc$)", "", 
                    gsub("^AML_sim_", "", 
                         gsub("\\.fcs$", "", basename(files_load))))
  sample_id
  
  group_id <- factor(gsub("_.*$", "", sample_id), levels = c("healthy", "CN", "CBF"))
  group_id
  
  patient_id <- factor(gsub("^.*_", "", sample_id))
  patient_id
  
  experiment_info <- data.frame(group_id, patient_id, sample_id)
  experiment_info
  
  # marker information
  
  # indices of all marker columns, lineage markers, and functional markers
  # (16 surface markers / 15 functional markers; see Levine et al. 2015, Supplemental 
  # Information, p. 4)
  cols_markers <- 11:41
  cols_lineage <- c(35, 29, 14, 30, 12, 26, 17, 33, 41, 32, 22, 40, 27, 37, 23, 39)
  cols_func <- setdiff(cols_markers, cols_lineage)
  
  stopifnot(all(sapply(seq_along(d_input), function(i) all(colnames(d_input[[i]]) == colnames(d_input[[1]])))))
  
  marker_name <- colnames(d_input[[1]])
  marker_name <- gsub("\\(.*$", "", marker_name)
  
  marker_class <- rep("none", length(marker_name))
  marker_class[cols_lineage] <- "type"
  marker_class[cols_func] <- "state"
  marker_class <- factor(marker_class, levels = c("type", "state", "none"))
  
  marker_info <- data.frame(marker_name, marker_class)
  marker_info
  
  
  
  
  ######################################
  # Additional pre-processing for Citrus
  ######################################
  
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
  
  out_Citrus_main[[th]] <- runtime_Citrus_main[[th]] <- vector("list", length(cond_names))
  names(out_Citrus_main[[th]]) <- names(runtime_Citrus_main[[th]]) <- cond_names
  
  
  for (j in 1:length(cond_names)) {
    
    # ----------------------------------------
    # Export transformed .fcs files for Citrus
    # ----------------------------------------
    
    ix_keep <- group_id %in% c("healthy", cond_names[j])
    
    sample_id_keep <- sample_id[ix_keep]
    group_id_keep <- droplevels(group_id[ix_keep])
    files_load_keep <- files_load[ix_keep]
    
    d_input_keep <- d_input[ix_keep]
    
    for (i in 1:length(sample_id_keep)) {
      path <- file.path(DIR_CITRUS_FILES, thresholds[th], cond_names[j], "data_transformed")
      filename <- file.path(path, gsub("\\.fcs$", "_transf.fcs", basename(files_load_keep[i])))
      write.FCS(d_input_keep[[i]], filename)
    }
    
    
    # ------------------------
    # Define inputs for Citrus
    # ------------------------
    
    # model type
    family <- "classification"
    modelTypes <- "glmnet"
    nFolds <- 1
    
    # feature type
    featureType <- "abundances"
    
    # define clustering markers
    clusteringColumns <- colnames(d_input[[1]])[cols_lineage]
    medianColumns <- NULL
    
    # number of cells and minimum cluster size
    fileSampleSize <- 5000
    minimumClusterSizePercent <- 0.001  # 0.1%
    
    # transformation: not required since already done above
    transformColumns <- NULL
    transformCofactor <- NULL
    scaleColumns <- NULL
    
    # experimental design
    labels <- group_id_keep
    
    # directories
    dataDirectory <- file.path(DIR_CITRUS_FILES, thresholds[th], cond_names[j], "data_transformed")
    outputDirectory <- file.path(DIR_CITRUS_FILES, thresholds[th], cond_names[j], "citrusOutput")
    
    # files
    fileList <- data.frame(defaultCondition = gsub("\\.fcs$", "_transf.fcs", basename(files_load_keep)))
    
    # number of processor threads
    n_cores <- 2
    
    
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
    
    runtime_Citrus_main[[th]][[j]] <- runtime_total
    
    
    
    
    ##############################
    # Return results at cell level
    ##############################
    
    # get differential clusters, match cells to clusters, and save results at cell level
    
    # note: Citrus does not give any continuous-valued scores, e.g. p-values or q-values,
    # so it is not possible to rank the selected clusters
    
    
    # number of cells per sample (including spike-in cells)
    n_cells <- sapply(d_input, nrow)
    
    # identify true spike-in cells (from 'spike' condition)
    is_spikein <- unlist(sapply(d_input, function(d) exprs(d)[, "spikein"]))
    stopifnot(length(is_spikein) == sum(n_cells))
    
    # select samples for this condition and healthy
    ix_keep_cnd <- group_id %in% c("healthy", cond_names[j])
    
    
    # differentially abundant clusters
    
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
      
      # match to indices in original data; taking into account subsampling and subset of
      # files in this condition
      which_files <- (1:length(ix_keep))[ix_keep]
      ix_files_keep <- as.numeric(as.character(factor(ix_files, labels = which_files)))
      
      for (z in seq_along(which_files)) {
        ix_events_z <- ix_events[ix_files_keep == which_files[z]]
        res_cells[files_rep == which_files[z]][ix_events_z] <- cells_subsampled[ix_files_keep == which_files[z]]
      }
      
      # cluster marker expression values (if required)
      #results$citrus.combinedFCSSet$data[results$citrus.foldClustering$allClustering$clusterMembership[[clus]], clusteringColumns]
    }
    
    # set all NA values to 0 to allow evaluation
    # (note: Citrus results are binary; 1 = cell selected, 0 = cell not selected)
    res_cells[is.na(res_cells)] <- 0
    
    
    # set up data frame with results and true spike-in status at cell level
    
    which_cnd <- rep(ix_keep_cnd, n_cells)
    is_spikein_cnd <- is_spikein[which_cnd]
    stopifnot(length(res_cells) == length(which_cnd), 
              length(res_cells[which_cnd]) == length(is_spikein_cnd))
    
    scores <- res_cells[which_cnd]
    
    # replace any NAs to ensure same set of cells is returned for all methods
    scores[is.na(scores)] <- 0
    
    stopifnot(length(scores) == length(is_spikein_cnd))
    
    # return values for this condition and healthy
    res <- data.frame(scores = scores, 
                      spikein = is_spikein_cnd)
    
    # store results
    out_Citrus_main[[th]][[j]] <- res
    
  }
}




#####################
# Save output objects
#####################

save(out_Citrus_main, runtime_Citrus_main, 
     file = file.path(DIR_RDATA, "outputs_AML_sim_Citrus_main.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_AML_sim_Citrus_main.txt"))
sessionInfo()
sink()



