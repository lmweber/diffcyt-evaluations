##########################################################################################
# Script to run methods
# 
# - method: CellCnn
# - data set: AML-sim
# 
# - main results
# 
# Lukas Weber, August 2017
##########################################################################################


# note: run from command line with 'Rscript <filename.R>'

# note: use all markers for CellCnn (using lineage markers only gives errors)


library(flowCore)
library(SummarizedExperiment)


DIR_BENCHMARK <- "../../../../../benchmark_data/AML_sim/data/main"
DIR_CELLCNN <- "../../../../../CellCnn/CellCnn"
DIR_CELLCNN_FILES <- "../../../../CellCnn_files"
DIR_RDATA <- "../../../../RData/AML_sim/main"
DIR_SESSION_INFO <- "../../../../session_info/AML_sim/main"




################################
# Loop to run for each threshold
################################

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")

# lists to store objects
out_CellCnn_main <- vector("list", length(thresholds))
names(out_CellCnn_main) <- thresholds




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
  
  # load data
  files_load <- c(files_healthy, files_CN, files_CBF)
  files_load
  
  d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)
  
  # sample IDs, group IDs, patient IDs
  sample_IDs <- gsub("(_[0-9]+pc$)|(_0\\.[0-9]+pc$)", "", 
                     gsub("^AML_sim_", "", 
                          gsub("\\.fcs$", "", basename(files_load))))
  sample_IDs
  
  group_IDs <- factor(gsub("_.*$", "", sample_IDs), levels = c("healthy", "CN", "CBF"))
  group_IDs
  
  patient_IDs <- factor(gsub("^.*_", "", sample_IDs))
  patient_IDs
  
  # check
  data.frame(sample_IDs, group_IDs, patient_IDs)
  
  # indices of all marker columns, lineage markers, and functional markers
  # (16 surface markers / 15 functional markers; see Levine et al. 2015, Supplemental 
  # Information, p. 4)
  cols_markers <- 11:41
  cols_lineage <- c(35, 29, 14, 30, 12, 26, 17, 33, 41, 32, 22, 40, 27, 37, 23, 39)
  cols_func <- setdiff(cols_markers, cols_lineage)
  
  
  # ------------------------------------
  # choose markers to use for clustering
  # ------------------------------------
  
  cols_clustering <- cols_markers
  
  
  # -------------------------------
  # marker names for CellCnn inputs
  # -------------------------------
  
  check <- c()
  for (i in 1:length(d_input)) {
    check[i] <- all(colnames(d_input[[i]]) == colnames(d_input[[1]]))
  }
  all(check)
  
  marker_names <- colnames(d_input[[1]])[cols_clustering]
  
  
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
  
  
  
  
  ##################
  # CellCnn pipeline
  ##################
  
  # note: run CellCnn separately for each condition: CN vs. healthy, CBF vs. healthy
  
  
  out_CellCnn_main[[th]] <- vector("list", length(cond_names))
  names(out_CellCnn_main[[th]]) <- cond_names
  
  
  for (j in 1:length(cond_names)) {
    
    # -----------------------------------------
    # Export transformed .fcs files for CellCnn
    # -----------------------------------------
    
    ix_keep <- group_IDs %in% c("healthy", cond_names[j])
    
    sample_IDs_keep <- sample_IDs[ix_keep]
    files_load_keep <- files_load[ix_keep]
    d_input_keep <- d_input[ix_keep]
    
    for (i in 1:length(sample_IDs_keep)) {
      path <- paste0(DIR_CELLCNN_FILES, "/data_transformed/AML_sim/main/", thresholds[th], "/", cond_names[j])
      filename <- file.path(path, gsub("\\.fcs$", "_transf.fcs", basename(files_load_keep[i])))
      write.FCS(d_input_keep[[i]], filename)
    }
    
    
    # --------------------------------
    # Generate input files for CellCnn
    # --------------------------------
    
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
    
    fn_samples <- paste0(DIR_CELLCNN_FILES, "/inputs/AML_sim/main/", thresholds[th], "/", cond_names[j], "/input_samples.csv")
    write.csv(df_samples, fn_samples, quote = FALSE, row.names = FALSE)
    
    # need to use 'write.table' to allow removing column names
    fn_markers <- paste0(DIR_CELLCNN_FILES, "/inputs/AML_sim/main/", thresholds[th], "/", cond_names[j], "/input_markers.csv")
    write.table(df_markers, fn_markers, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    
    # -----------------------------
    # Run CellCnn from command line
    # -----------------------------
    
    # for installation instructions and examples see: https://github.com/eiriniar/CellCnn
    
    
    # command to run CellCnn analysis
    cmd <- paste("python", paste0(DIR_CELLCNN, "/cellCnn/run_analysis.py"), 
                 paste0("-f ", DIR_CELLCNN_FILES, "/inputs/AML_sim/main/", thresholds[th], "/", cond_names[j], "/input_samples.csv"), 
                 paste0("-m ", DIR_CELLCNN_FILES, "/inputs/AML_sim/main/", thresholds[th], "/", cond_names[j], "/input_markers.csv"), 
                 paste0("-i ", DIR_CELLCNN_FILES, "/data_transformed/AML_sim/main/", thresholds[th], "/", cond_names[j], "/"), 
                 paste0("-o ", DIR_CELLCNN_FILES, "/out_CellCnn/AML_sim/main/", thresholds[th], "/", cond_names[j], "/"), 
                 "--export_csv", 
                 paste("--group_a", "Healthy", "--group_b", cond_names[j]))
    
    # run from command line
    runtime_analysis <- system.time(
      system(cmd)
    )
    
    runtime_analysis
    
    sink(paste0(DIR_CELLCNN_FILES, "/runtime/AML_sim/main/", thresholds[th], "/", cond_names[j], "/runtime_analysis.txt"))
    runtime_analysis
    sink()
    
    
    # command to export selected cells
    cmd <- paste("python", paste0(DIR_CELLCNN, "/cellCnn/run_analysis.py"), 
                 paste0("-f ", DIR_CELLCNN_FILES, "/inputs/AML_sim/main/", thresholds[th], "/", cond_names[j], "/input_samples.csv"), 
                 paste0("-m ", DIR_CELLCNN_FILES, "/inputs/AML_sim/main/", thresholds[th], "/", cond_names[j], "/input_markers.csv"), 
                 paste0("-i ", DIR_CELLCNN_FILES, "/data_transformed/AML_sim/main/", thresholds[th], "/", cond_names[j], "/"), 
                 paste0("-o ", DIR_CELLCNN_FILES, "/out_CellCnn/AML_sim/main/", thresholds[th], "/", cond_names[j], "/"), 
                 "--plot", 
                 paste("--group_a", "Healthy", "--group_b", cond_names[j]), 
                 "--filter_response_thres 0.3 --load_results --export_selected_cells")
    
    # run from command line
    runtime_select <- system.time(
      system(cmd)
    )
    
    runtime_select
    
    sink(paste0(DIR_CELLCNN_FILES, "/runtime/AML_sim/main/", thresholds[th], "/", cond_names[j], "/runtime_select.txt"))
    runtime_select
    sink()
    
    
    
    
    ##############################
    # Return results at cell level
    ##############################
    
    # Note: CellCnn returns continuous 'scores' at the cell level, indicating the
    # likelihood of each cell belonging to each detected 'filter' (population). If there
    # are multiple detected filters, we sum the scores to give a single total score per
    # cell.
    
    
    # number of cells per sample (including spike-in cells)
    n_cells <- sapply(d_input, nrow)
    
    # spike-in status for each cell
    is_spikein <- unlist(sapply(d_input, function(d) exprs(d)[, "spikein"]))
    stopifnot(length(is_spikein) == sum(n_cells))
    
    # select samples for this condition
    ix_keep_cnd <- group_IDs == cond_names[j]
    
    
    # CellCnn output files
    
    path_out <- paste0(DIR_CELLCNN_FILES, "/out_CellCnn/AML_sim/main/", thresholds[th], "/", cond_names[j], "/selected_cells")
    
    # if no files exist, CellCnn did not run correctly; return all zeros in this case
    if (length(list.files(path_out)) == 0) {
      filter_continuous_cnd <- list(rep(0, sum(n_cells[ix_keep_cnd])))
    }
    
    files_cnd <- paste0(path_out, "/", gsub("\\.fcs$", "", basename(files_load[ix_keep_cnd])), "_transf_selected_cells.csv")
    
    # get cells in selected filters for this condition
    filter_continuous_cnd <- vector("list", length(files_cnd))
    
    for (f in 1:length(files_cnd)) {
      d <- try(read.csv(files_cnd[f]), silent = TRUE)
      
      if (!(class(d) == "try-error")) {
        # note: if there are multiple filters for one sample, combine them using 'rowSums'
        # (filters are stored in odd-numbered columns of the .csv file)
        ix <- seq(1, ncol(d), by = 2)
        filt_sum <- rowSums(d[, ix, drop = FALSE])
        filter_continuous_cnd[[f]] <- filt_sum
        
      } else {
        # if .csv file is empty, fill with zeros instead
        filter_continuous_cnd[[f]] <- rep(0, n_cells[ix_keep_cnd][[f]])
      }
    }
    
    filter_continuous_cnd <- unlist(filter_continuous_cnd)
    
    
    # set up data frame with results and true spike-in status at cell level
    
    which_cnd <- rep(ix_keep_cnd, n_cells)
    is_spikein_cnd <- is_spikein[which_cnd]
    stopifnot(length(filter_continuous_cnd) == length(is_spikein_cnd))
    
    scores <- filter_continuous_cnd
    
    # replace any NAs to ensure same set of cells is returned for all methods
    scores[is.na(scores)] <- 0
    
    # return values for this condition only
    res <- data.frame(scores = scores, 
                      spikein = is_spikein_cnd)
    
    # store results
    out_CellCnn_main[[th]][[j]] <- res
    
  }
}




#####################
# Save output objects
#####################

save(out_CellCnn_main, file = file.path(DIR_RDATA, "/outputs_AML_sim_CellCnn_main.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "/session_info_AML_sim_CellCnn_main.txt"))
sessionInfo()
sink()



