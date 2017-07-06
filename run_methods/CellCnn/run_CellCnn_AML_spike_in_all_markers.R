##########################################################################################
# Script to run methods
# 
# - method: CellCnn
# - data set: AML-spike-in
# 
# Lukas Weber, July 2017
##########################################################################################


library(flowCore)
library(SummarizedExperiment)
library(diffcyt)




################################
# Loop to run for each threshold
################################

# spike-in thresholds (must match filenames)
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")

# lists to store objects
is_spikein <- vector("list", length(thresholds))
names(is_spikein) <- thresholds




for (th in 1:length(thresholds)) {
  
  ##############################
  # Load data and pre-processing
  ##############################
  
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
  
  
  # ------------------------------------------------------
  # remove spike-in indicator columns and store separately
  # ------------------------------------------------------
  
  # 'AML-spike-in' data set includes indicator columns for spike-in cells in conditions
  # 'CN' and 'CBF'; these need to be removed and stored separately
  
  is_spikein[[th]] <- vector("list", length(sample_IDs))
  names(is_spikein[[th]]) <- sample_IDs
  
  for (i in 1:length(sample_IDs)) {
    exprs_i <- exprs(d_input[[i]])
    if (group_IDs[i] == "healthy") {
      is_spikein[[th]][[i]] <- rep(0, nrow(exprs_i))
    } else {
      is_spikein[[th]][[i]] <- exprs_i[, "spikein"]
      # remove column from expression matrix
      exprs(d_input[[i]]) <- exprs_i[, -match("spikein", colnames(exprs_i))]
    }
  }
  
  
  # -------------------------------
  # marker names for CellCnn inputs
  # -------------------------------
  
  markers_to_use <- cols_markers
  
  check <- c()
  for (i in 1:length(d_input)) {
    check[i] <- all(colnames(d_input[[i]]) == colnames(d_input[[1]]))
  }
  all(check)
  
  marker_names <- colnames(d_input[[1]])[markers_to_use]
  
  
  
  
  #################################
  # Transform and export .fcs files
  #################################
  
  # --------------
  # transform data
  # --------------
  
  cofactor <- 5
  
  d_input <- lapply(d_input, function(d) {
    e <- exprs(d)
    e[, cols_markers] <- asinh(e[, cols_markers] / cofactor)
    e
  })
  
  
  # -----------------------------------------
  # export transformed .fcs files for CellCnn
  # -----------------------------------------
  
  for (i in 1:length(sample_IDs)) {
    path <- paste0("../../../CellCnn_files/data_transformed/AML_spike_in/", thresholds[th])
    filename <- file.path(path, gsub("\\.fcs$", "_transf.fcs", basename(files_load[i])))
    write.FCS(flowFrame(d_input[[i]]), filename)
  }
  
  
  
  
  ######################################
  # Generate input files and run CellCnn
  ######################################
  
  # note: generate inputs and run CellCnn separately for each condition: CN vs. healthy,
  # CBF vs. healthy
  
  
  for (j in 1:length(cond_names)) {
    
    ##################################
    # Generate input files for CellCnn
    ##################################
    
    # generate .csv files with input arguments for CellCnn (in required format)
    
    files_transf <- gsub("\\.fcs$", "_transf.fcs", basename(files_load))
    files_transf <- files_transf[group_IDs %in% c("healthy", cond_names[j])]
    files_transf
    
    
    # create data frame of sample names and conditions (for CellCnn input .csv file)
    
    label <- group_IDs[group_IDs %in% c("healthy", cond_names[j])]
    label <- as.numeric(droplevels(label)) - 1
    label
    
    df_samples <- data.frame(fcs_filename = files_transf, label = label)
    df_samples
    
    # re-arrange alphabetically (otherwise CellCnn reads input files in incorrect order)
    df_samples <- df_samples[order(df_samples$fcs_filename), ]
    df_samples
    
    
    # create data frame of column names (markers) (for CellCnn input .csv file)
    
    df_markers <- t(data.frame(marker_names))
    df_markers
    
    
    # save as .csv files
    
    fn_samples <- paste0("../../../CellCnn_files/inputs/AML_spike_in/", thresholds[th], "/input_samples.csv")
    write.csv(df_samples, fn_samples, quote = FALSE, row.names = FALSE)
    
    # need to use 'write.table' to allow removing column names
    fn_markers <- paste0("../../../CellCnn_files/inputs/AML_spike_in/", thresholds[th], "/input_markers.csv")
    write.table(df_markers, fn_markers, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    
    
    
    ###############################
    # Run CellCnn from command line
    ###############################
    
    # for installation instructions and examples see: https://github.com/eiriniar/CellCnn
    
    DIR_CellCnn <- "../../../../CellCnn/CellCnn/"
    
    
    # command to run main analysis
    cmd <- paste("python", paste0(DIR_CellCnn, "cellCnn/run_analysis.py"), 
                 paste0("-f ../../../CellCnn_files/inputs/AML_spike_in/", thresholds[th], "/input_samples.csv"), 
                 paste0("-m ../../../CellCnn_files/inputs/AML_spike_in/", thresholds[th], "/input_markers.csv"), 
                 paste0("-i ../../../CellCnn_files/data_transformed/AML_spike_in/", thresholds[th], "/"), 
                 paste0("-o ../../../CellCnn_files/out_CellCnn/AML_spike_in/", thresholds[th], "/"), 
                 "--export_csv", 
                 paste("--group_a", "Healthy", "--group_b", cond_names[j]))
    
    # run from command line
    runtime_main <- system.time(
      system(cmd)
    )
    
    runtime_main
    
    sink(paste0("../../../CellCnn_files/runtime/AML_spike_in/", thresholds[th], "/runtime_main.txt"))
    runtime_main
    sink()
    
    
    # command to export selected cells
    cmd <- paste("python", paste0(DIR_CellCnn, "cellCnn/run_analysis.py"), 
                 paste0("-f ../../../CellCnn_files/inputs/AML_spike_in/", thresholds[th], "/input_samples.csv"), 
                 paste0("-m ../../../CellCnn_files/inputs/AML_spike_in/", thresholds[th], "/input_markers.csv"), 
                 paste0("-i ../../../CellCnn_files/data_transformed/AML_spike_in/", thresholds[th], "/"), 
                 paste0("-o ../../../CellCnn_files/out_CellCnn/AML_spike_in/", thresholds[th], "/"), 
                 "--plot", 
                 paste("--group_a", "Healthy", "--group_b", cond_names[j]), 
                 "--filter_response_thres 0.3 --load_results --export_selected_cells")
    
    # run from command line
    runtime_select <- system.time(
      system(cmd)
    )
    
    runtime_select
    
    sink(paste0("../../../CellCnn_files/runtime/AML_spike_in/", thresholds[th], "/runtime_select.txt"))
    runtime_select
    sink()
  }
  
}




#####################
# Save output objects
#####################

save.image("../../../RData/output_CellCnn_AML_spike_in_all_markers.RData")




#####################
# Session information
#####################

sink("../../../session_info/session_info_CellCnn_AML_spike_in_all_markers.txt")
sessionInfo()
sink()



