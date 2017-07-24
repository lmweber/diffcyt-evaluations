##########################################################################################
# Script to run methods
# 
# - method: Citrus-all-markers
# - data set: AML-spike-in
# 
# Lukas Weber, July 2017
##########################################################################################


library(flowCore)
library(citrus)




################################
# Loop to run for each threshold
################################

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")

# lists to store objects
out_Citrus_all_markers <- vector("list", length(thresholds))
names(out_Citrus_all_markers) <- thresholds




for (th in 1:length(thresholds)) {
  
  ######################################
  # Load data, pre-processing, transform
  ######################################
  
  # ---------
  # load data
  # ---------
  
  # filenames
  files_healthy <- list.files("../../../../benchmark_data/AML_spike_in/data/healthy", 
                              pattern = "\\.fcs$", full.names = TRUE)
  files_CN <- list.files("../../../../benchmark_data/AML_spike_in/data/CN", 
                         pattern = paste0("_", thresholds[th], "\\.fcs$"), full.names = TRUE)
  files_CBF <- list.files("../../../../benchmark_data/AML_spike_in/data/CBF", 
                          pattern = paste0("_", thresholds[th], "\\.fcs$"), full.names = TRUE)
  
  # load data
  files_load <- c(files_healthy, files_CN, files_CBF)
  files_load
  
  d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)
  
  # sample IDs, group IDs, block IDs
  sample_IDs <- gsub("(_[0-9]+pc$)|(_0\\.[0-9]+pc$)", "", 
                     gsub("^AML_spike_in_", "", 
                          gsub("\\.fcs$", "", basename(files_load))))
  sample_IDs
  
  group_IDs <- gsub("_.*$", "", sample_IDs)
  group_IDs
  
  block_IDs <- gsub("^.*_", "", sample_IDs)
  block_IDs
  
  # set group_IDs reference level (for differential tests)
  group_IDs <- factor(group_IDs, levels = c("healthy", "CN", "CBF"))
  group_IDs
  
  # check all match correctly
  data.frame(sample_IDs, group_IDs, block_IDs)
  
  # indices of all marker columns, lineage markers, and functional markers
  # (16 surface markers / 15 functional markers; see Levine et al. 2015, Supplemental 
  # Information, p. 4)
  cols_markers <- 11:41
  cols_lineage <- c(35, 29, 14, 30, 12, 26, 17, 33, 41, 32, 22, 40, 27, 37, 23, 39)
  cols_func <- setdiff(cols_markers, cols_lineage)
  
  
  # ---------------------------
  # choose which markers to use
  # ---------------------------
  
  cols_to_use <- cols_markers
  
  
  # --------------
  # transform data
  # --------------
  
  # 'asinh' transform with 'cofactor' = 5 (see Bendall et al. 2011, Supp. Fig. S2)
  
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
  
  out_Citrus_all_markers[[th]] <- vector("list", length(cond_names))
  names(out_Citrus_all_markers[[th]]) <- cond_names
  
  
  for (j in 1:length(cond_names)) {
    
    # ----------------------------------------
    # Export transformed .fcs files for Citrus
    # ----------------------------------------
    
    ix_keep <- group_IDs %in% c("healthy", cond_names[j])
    
    sample_IDs_keep <- sample_IDs[ix_keep]
    group_IDs_keep <- droplevels(group_IDs[ix_keep])
    files_load_keep <- files_load[ix_keep]
    
    d_input_keep <- d_input[ix_keep]
    
    for (i in 1:length(sample_IDs_keep)) {
      path <- paste0("../../../Citrus_files/data_transformed/AML_spike_in/all_markers/", thresholds[th], "/", cond_names[j])
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
    
    # number of cells and minimum cluster size
    fileSampleSize <- 5000
    minimumClusterSizePercent <- 0.01
    
    # columns for clustering
    clusteringColumns <- colnames(d_input[[1]])[cols_to_use]
    medianColumns <- NULL
    
    # experimental design
    labels <- group_IDs_keep
    
    # transformation: not required since already done above
    transformColumns <- NULL
    transformCofactor <- NULL
    scaleColumns <- NULL
    
    # directories
    dataDirectory <- paste0("../../../Citrus_files/data_transformed/AML_spike_in/all_markers/", thresholds[th], "/", cond_names[j])
    outputDirectory <- file.path(dataDirectory, "citrusOutput")
    
    # files
    fileList <- data.frame(defaultCondition = gsub("\\.fcs$", "_transf.fcs", basename(files_load_keep)))
    
    
    # run Citrus
    
    set.seed(123)
    
    runtime_Citrus <- system.time(
      results <- citrus.full(
        fileList = fileList, 
        labels = labels, 
        clusteringColumns = clusteringColumns, 
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
        scaleColumns = scaleColumns, 
        medianColumns = medianColumns
      )
    )
    
    print(runtime_Citrus)  ## ~2 minutes on laptop
    
    # Citrus plots
    plot(results, outputDirectory)
    
    
    
    
    ##############################
    # Return results at cell level
    ##############################
    
    # get differential clusters, match cells to clusters, and save results at cell level
    
    # note: Citrus does not give any continuous-valued scores, e.g. p-values or q-values,
    # so it is not possible to rank the selected clusters
    
    
    # number of cells per sample (including spike-in cells)
    n_cells <- sapply(d_input, nrow)
    
    # spike-in status for each cell
    is_spikein <- unlist(sapply(d_input, function(d) exprs(d)[, "spikein"]))
    stopifnot(length(is_spikein) == sum(n_cells))
    
    # select samples for this condition
    ix_keep_cnd <- group_IDs == cond_names[j]
    
    
    # differentially abundant clusters
    
    clusters <- as.numeric(results$conditionRegressionResults$defaultCondition$glmnet$differentialFeatures[["cv.min"]][["clusters"]])
    
    
    # match clusters to cells
    
    res_cells <- rep(NA, sum(n_cells))
    files_rep <- rep(seq_along(sample_IDs), n_cells)
    
    for (i in seq_along(clusters)) {
      
      clus <- clusters[i]
      
      # identify cells in cluster
      ix_clus <- results$citrus.foldClustering$allClustering$clusterMembership[[clus]]
      
      data <- results$citrus.combinedFCSSet$data
      ix_events <- data[, "fileEventNumber"]
      ix_files <- data[, "fileId"]
      
      cells_subsampled <- rep(NA, nrow(data))
      cells_subsampled[ix_clus] <- 1  ## 1 = cell is in differential cluster
      
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
    stopifnot(length(res_cells) == length(which_cnd), length(res_cells[which_cnd]) == length(is_spikein_cnd))
    
    scores <- res_cells[which_cnd]
    
    # replace any NAs to ensure same set of cells is returned for all methods
    scores[is.na(scores)] <- 0
    
    # return values for this condition only
    res <- data.frame(scores = scores, 
                      spikein = is_spikein_cnd)
    
    # store results
    out_Citrus_all_markers[[th]][[j]] <- res
    
  }
}




#####################
# Save output objects
#####################

save(out_Citrus_all_markers, file = "../../../RData/AML_spike_in/outputs_AML_spike_in_Citrus_all_markers.RData")




#####################
# Session information
#####################

sink("../../../session_info/AML_spike_in/session_info_AML_spike_in_Citrus_all_markers.txt")
sessionInfo()
sink()



