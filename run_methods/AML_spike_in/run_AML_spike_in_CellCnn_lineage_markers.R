##########################################################################################
# Script to run methods
# 
# - method: CellCnn-lineage-markers
# - data set: AML-spike-in
# 
# Lukas Weber, July 2017
##########################################################################################


# note: run from command line with 'Rscript <filename.R>'


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
files_load_thresholds <- is_spikein_thresholds <- vector("list", length(thresholds))
names(files_load_thresholds) <- names(is_spikein_thresholds) <- thresholds




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
  
  
  # ---------------------------
  # choose which markers to use
  # ---------------------------
  
  cols_to_use <- cols_lineage
  
  
  # -----------------------------------
  # store filenames and spike-in status
  # -----------------------------------
  
  # filenames
  files_load_thresholds[[th]] <- files_load
  
  # spike-in status for each cell
  is_spikein_thresholds[[th]] <- vector("list", length(sample_IDs))
  names(is_spikein_thresholds[[th]]) <- sample_IDs
  
  for (i in 1:length(sample_IDs)) {
    exprs_i <- exprs(d_input[[i]])
    is_spikein_thresholds[[th]][[i]] <- exprs_i[, "spikein"]
  }
  
  
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
  
  cofactor <- 5
  
  d_input <- lapply(d_input, function(d) {
    e <- exprs(d)
    e[, cols_markers] <- asinh(e[, cols_markers] / cofactor)
    e
  })
  
  
  
  
  ##################################################################
  # Export .fcs files, generate CellCnn input files, and run CellCnn
  ##################################################################
  
  # note: run CellCnn separately for each condition: CN vs. healthy, CBF vs. healthy
  
  
  
  for (j in 1:length(cond_names)) {
    
    ###########################################
    # Export transformed .fcs files for CellCnn
    ###########################################
    
    ix_keep <- group_IDs %in% c("healthy", cond_names[j])
    
    sample_IDs_keep <- sample_IDs[ix_keep]
    files_load_keep <- files_load[ix_keep]
    d_input_keep <- d_input[ix_keep]
    
    for (i in 1:length(sample_IDs_keep)) {
      path <- paste0("../../../CellCnn_files/AML_spike_in/lineage_markers/data_transformed/", thresholds[th], "/", cond_names[j])
      filename <- file.path(path, gsub("\\.fcs$", "_transf.fcs", basename(files_load_keep[i])))
      write.FCS(flowFrame(d_input_keep[[i]]), filename)
    }
    
    
    
    
    ##################################
    # Generate input files for CellCnn
    ##################################
    
    # generate .csv files with input arguments for CellCnn (in required format)
    
    files_transf <- gsub("\\.fcs$", "_transf.fcs", basename(files_load))
    files_transf <- files_transf[ix_keep]
    files_transf
    
    
    # create data frame of sample names and conditions (for CellCnn input .csv file)
    
    label <- group_IDs[ix_keep]
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
    
    fn_samples <- paste0("../../../CellCnn_files/AML_spike_in/lineage_markers/inputs/", thresholds[th], "/", cond_names[j], "/input_samples.csv")
    write.csv(df_samples, fn_samples, quote = FALSE, row.names = FALSE)
    
    # need to use 'write.table' to allow removing column names
    fn_markers <- paste0("../../../CellCnn_files/AML_spike_in/lineage_markers/inputs/", thresholds[th], "/", cond_names[j], "/input_markers.csv")
    write.table(df_markers, fn_markers, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    
    
    
    ###############################
    # Run CellCnn from command line
    ###############################
    
    # for installation instructions and examples see: https://github.com/eiriniar/CellCnn
    
    DIR_CellCnn <- "../../../../CellCnn/CellCnn/"
    
    
    # command to run main analysis
    cmd <- paste("python", paste0(DIR_CellCnn, "cellCnn/run_analysis.py"), 
                 paste0("-f ../../../CellCnn_files/AML_spike_in/lineage_markers/inputs/", thresholds[th], "/", cond_names[j], "/input_samples.csv"), 
                 paste0("-m ../../../CellCnn_files/AML_spike_in/lineage_markers/inputs/", thresholds[th], "/", cond_names[j], "/input_markers.csv"), 
                 paste0("-i ../../../CellCnn_files/AML_spike_in/lineage_markers/data_transformed/", thresholds[th], "/", cond_names[j], "/"), 
                 paste0("-o ../../../CellCnn_files/AML_spike_in/lineage_markers/out_CellCnn/", thresholds[th], "/", cond_names[j], "/"), 
                 "--export_csv", 
                 paste("--group_a", "Healthy", "--group_b", cond_names[j]))
    
    # run from command line
    runtime_main <- system.time(
      system(cmd)
    )
    
    runtime_main
    
    sink(paste0("../../../CellCnn_files/AML_spike_in/lineage_markers/runtime/", thresholds[th], "/", cond_names[j], "/runtime_main.txt"))
    runtime_main
    sink()
    
    
    # command to export selected cells
    cmd <- paste("python", paste0(DIR_CellCnn, "cellCnn/run_analysis.py"), 
                 paste0("-f ../../../CellCnn_files/AML_spike_in/lineage_markers/inputs/", thresholds[th], "/", cond_names[j], "/input_samples.csv"), 
                 paste0("-m ../../../CellCnn_files/AML_spike_in/lineage_markers/inputs/", thresholds[th], "/", cond_names[j], "/input_markers.csv"), 
                 paste0("-i ../../../CellCnn_files/AML_spike_in/lineage_markers/data_transformed/", thresholds[th], "/", cond_names[j], "/"), 
                 paste0("-o ../../../CellCnn_files/AML_spike_in/lineage_markers/out_CellCnn/", thresholds[th], "/", cond_names[j], "/"), 
                 "--plot", 
                 paste("--group_a", "Healthy", "--group_b", cond_names[j]), 
                 "--filter_response_thres 0.3 --load_results --export_selected_cells")
    
    # run from command line
    runtime_select <- system.time(
      system(cmd)
    )
    
    runtime_select
    
    sink(paste0("../../../CellCnn_files/AML_spike_in/lineage_markers/runtime/", thresholds[th], "/", cond_names[j], "/runtime_select.txt"))
    runtime_select
    sink()
  }
  
}




#####################
# Save output objects
#####################

save.image("../../../RData/AML_spike_in/outputs_AML_spike_in_CellCnn_lineage_markers.RData")




#####################
# Session information
#####################

sink("../../../session_info/AML_spike_in/session_info_AML_spike_in_CellCnn_lineage_markers.txt")
sessionInfo()
sink()



