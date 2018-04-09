##########################################################################################
# Script to run methods
# 
# - method: CellCnn
# - data set: AML-sim
# 
# - main results
# 
# Lukas Weber, April 2018
##########################################################################################


# note: run from command line with 'Rscript <filename.R>'

# note: CellCnn returns errors and does not complete when using 'cell type' markers only;
# use all markers instead


library(flowCore)
library(SummarizedExperiment)


DIR_BENCHMARK <- "../../../../../benchmark_data/AML_sim/data/main"
DIR_CELLCNN <- "../../../../../CellCnn/CellCnn"
DIR_MINICONDA2 <- "~/miniconda2/bin"
DIR_CELLCNN_FILES <- "../../../../CellCnn_files/AML_sim/main"
DIR_RDATA <- "../../../../RData/AML_sim/comparisons_CellCnn"
DIR_SESSION_INFO <- "../../../../session_info/AML_sim/comparisons_CellCnn"




##############################
# Delete previous output files
##############################

# delete output files from previous CellCnn runs (but leave directory structure intact)

cmd_clean <- paste("find", DIR_CELLCNN_FILES, "-type f -delete")

system(cmd_clean)




################################
# Loop to run for each threshold
################################

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# condition names
cond_names <- c("CN", "CBF")

# lists to store objects and runtime
out_CellCnn_main <- runtime_CellCnn_main <- vector("list", length(thresholds))
names(out_CellCnn_main) <- names(runtime_CellCnn_main) <- thresholds




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
  marker_class[cols_lineage] <- "cell_type"
  marker_class[cols_func] <- "cell_state"
  marker_class <- factor(marker_class, levels = c("cell_type", "cell_state", "none"))
  
  marker_info <- data.frame(marker_name, marker_class)
  marker_info
  
  
  
  
  #######################################
  # Additional pre-processing for CellCnn
  #######################################
  
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
  
  markers_to_use <- colnames(d_input[[1]])[cols_to_use]
  
  
  # --------------
  # transform data
  # --------------
  
  # not required, since CellCnn automatically transforms data (cannot be disabled)
  
  
  
  
  ##################
  # CellCnn pipeline
  ##################
  
  # note: run CellCnn separately for each condition: CN vs. healthy, CBF vs. healthy
  
  
  out_CellCnn_main[[th]] <- runtime_CellCnn_main[[th]] <- vector("list", length(cond_names))
  names(out_CellCnn_main[[th]]) <- names(runtime_CellCnn_main[[th]]) <- cond_names
  
  
  for (j in 1:length(cond_names)) {
    
    # -----------------------------
    # Export .fcs files for CellCnn
    # -----------------------------
    
    ix_keep <- group_id %in% c("healthy", cond_names[j])
    
    sample_id_keep <- sample_id[ix_keep]
    sample_id_keep <- droplevels(sample_id[ix_keep])
    files_load_keep <- files_load[ix_keep]
    
    d_input_keep <- d_input[ix_keep]
    
    for (i in 1:length(sample_id_keep)) {
      filename <- file.path(DIR_CELLCNN_FILES, thresholds[th], cond_names[j], "data", basename(files_load_keep[i]))
      write.FCS(d_input_keep[[i]], filename)
    }
    
    
    # --------------------------------
    # Generate input files for CellCnn
    # --------------------------------
    
    # generate .csv files with input arguments for CellCnn (in required format)
    
    files_in <- basename(files_load_keep)
    files_in
    
    
    # create data frame of sample names and conditions (for CellCnn input .csv file)
    
    label <- sample_id_keep
    label <- as.numeric(label) - 1
    label
    
    df_samples <- data.frame(fcs_filename = files_in, label = label)
    df_samples
    
    # re-arrange alphabetically (otherwise CellCnn reads input files in incorrect order)
    df_samples <- df_samples[order(df_samples$fcs_filename), ]
    df_samples
    
    
    # create data frame of column names (markers) (for CellCnn input .csv file)
    
    df_markers <- t(data.frame(markers_to_use))
    df_markers
    
    
    # save as .csv files
    
    fn_samples <- file.path(DIR_CELLCNN_FILES, thresholds[th], cond_names[j], "inputs", "input_samples.csv")
    write.csv(df_samples, fn_samples, quote = FALSE, row.names = FALSE)
    
    # need to use 'write.table' to allow removing column names
    fn_markers <- file.path(DIR_CELLCNN_FILES, thresholds[th], cond_names[j], "inputs", "input_markers.csv")
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
    # (1) '--no_arcsinh' argument to disable arcsinh transform (but does not seem to work)
    # (2) '--ncell 300' argument to increase size of each training set [this should be at
    # least 'n_cells_min = (avg. no. cells per sample) * (no. samples per condition) / 1000';
    # if no memory constraints then increase to 5 to 10 times 'n_cells_min']
    # (3) '--subset_selection outlier' for extremely rare populations
    
    # note: use option '--subset_selection outlier' for extremely rare populations (threshold 0.1%)
    cmd_outlier <- ifelse(th == 3, "--subset_selection outlier ", "")
    
    # command to activate virtual environment
    # note: see CellCnn installation instructions for how to set up virtual environment
    # note: 'source activate cellcnn_env' works from command line, but also need to provide
    # path to 'activate' when using 'system' command from R
    cmd_env <- paste("source", file.path(DIR_MINICONDA2, "activate"), "cellcnn_env")
    
    # command to run CellCnn analysis
    cmd_run <- paste("python", paste0(DIR_CELLCNN, "/cellCnn/run_analysis.py"), 
                     paste0("-f ", DIR_CELLCNN_FILES, "/", thresholds[th], "/", cond_names[j], "/inputs/input_samples.csv"), 
                     paste0("-m ", DIR_CELLCNN_FILES, "/", thresholds[th], "/", cond_names[j], "/inputs/input_markers.csv"), 
                     paste0("-i ", DIR_CELLCNN_FILES, "/", thresholds[th], "/", cond_names[j], "/data/"), 
                     paste0("-o ", DIR_CELLCNN_FILES, "/", thresholds[th], "/", cond_names[j], "/out_CellCnn/"), 
                     paste0(cmd_outlier, "--ncell 300 --export_csv"), 
                     #"--no_arcsinh",  ## currently not working correctly
                     paste("--group_a", "Healthy", "--group_b", cond_names[j]))
    
    
    # command to export selected cells
    cmd_export <- paste("python", paste0(DIR_CELLCNN, "/cellCnn/run_analysis.py"), 
                        paste0("-f ", DIR_CELLCNN_FILES, "/", thresholds[th], "/", cond_names[j], "/inputs/input_samples.csv"), 
                        paste0("-m ", DIR_CELLCNN_FILES, "/", thresholds[th], "/", cond_names[j], "/inputs/input_markers.csv"), 
                        paste0("-i ", DIR_CELLCNN_FILES, "/", thresholds[th], "/", cond_names[j], "/data/"), 
                        paste0("-o ", DIR_CELLCNN_FILES, "/", thresholds[th], "/", cond_names[j], "/out_CellCnn/"), 
                        paste("--group_a", "Healthy", "--group_b", cond_names[j]), 
                        "--filter_response_thres 0.3 --load_results --export_selected_cells")
    
    # combine commands (in a single environment)
    cmd_combined <- paste(c(cmd_env, cmd_run, cmd_export), collapse = "; ")
    
    # run from command line
    runtime_combined <- system.time(
      system(cmd_combined)
    )
    
    
    # runtime
    runtime_total <- runtime_combined[["elapsed"]]
    print(runtime_total)
    
    runtime_CellCnn_main[[th]][[j]] <- runtime_total
    
    
    
    
    ##############################
    # Return results at cell level
    ##############################
    
    # Note: CellCnn returns continuous 'scores' at the cell level, representing the likelihood
    # of each cell belonging to each detected 'filter' (population). If there are multiple
    # detected filters, we sum the scores to give a single total score per cell.
    
    
    # number of cells per sample (including spike-in cells)
    n_cells <- sapply(d_input, nrow)
    
    # identify true spike-in cells (from 'spike' condition)
    is_spikein <- unlist(sapply(d_input, function(d) exprs(d)[, "spikein"]))
    stopifnot(length(is_spikein) == sum(n_cells))
    
    # select samples for this condition and healthy
    ix_keep_cnd <- sample_id %in% c("healthy", cond_names[j])
    
    
    # CellCnn output files
    
    path_out <- file.path(DIR_CELLCNN_FILES, thresholds[th], cond_names[j], "out_CellCnn", "selected_cells")
    
    # if no files exist, CellCnn did not run correctly; return all zeros in this case
    if (length(list.files(path_out)) == 0) {
      filter_continuous <- list(rep(0, sum(n_cells[ix_keep_cnd])))
    }
    
    files_cells <- paste0(path_out, "/", gsub("\\.fcs$", "", basename(files_load_keep)), "_selected_cells.csv")
    
    # get cells in selected filters for this condition
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
        filter_continuous[[f]] <- rep(0, n_cells[ix_keep_cnd][[f]])
      }
    }
    
    filter_continuous <- unlist(filter_continuous)
    
    
    # set up data frame with results and true spike-in status at cell level
    
    which_cnd <- rep(ix_keep_cnd, n_cells)
    is_spikein_cnd <- is_spikein[which_cnd]
    stopifnot(length(filter_continuous) == sum(which_cnd), 
              length(filter_continuous) == length(is_spikein_cnd))
    
    scores <- filter_continuous
    
    # replace any NAs to ensure same set of cells is returned for all methods
    scores[is.na(scores)] <- 0
    
    stopifnot(length(scores) == length(is_spikein_cnd))
    
    # return values for this condition and healthy
    res <- data.frame(scores = scores, 
                      spikein = is_spikein_cnd)
    
    # store results
    out_CellCnn_main[[th]][[j]] <- res
    
  }
}




#####################
# Save output objects
#####################

save(out_CellCnn_main, runtime_CellCnn_main, 
     file = file.path(DIR_RDATA, "outputs_AML_sim_CellCnn_main.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_AML_sim_CellCnn_main.txt"))
sessionInfo()
sink()



