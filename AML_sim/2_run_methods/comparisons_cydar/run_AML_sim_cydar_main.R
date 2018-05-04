##########################################################################################
# Script to run methods
# 
# - method: cydar
# - data set: AML-sim
# 
# - main results
# 
# Lukas Weber, May 2018
##########################################################################################


library(cydar)
library(flowCore)
library(ncdfFlow)
library(edgeR)
library(SummarizedExperiment)
library(Rtsne)


DIR_BENCHMARK <- "../../../../../benchmark_data/AML_sim/data/main"
DIR_CYDAR_FILES <- "../../../../cydar_files/AML_sim/main"
DIR_RDATA <- "../../../../RData/AML_sim/comparisons_cydar"
DIR_SESSION_INFO <- "../../../../session_info/AML_sim/comparisons_cydar"




##############################
# Delete previous output files
##############################

# delete output files from previous cydar runs (but leave directory structure intact)

cmd_clean <- paste("find", DIR_CYDAR_FILES, "-type f -delete")

system(cmd_clean)




################################
# Loop to run for each threshold
################################

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")

# condition names
cond_names <- c("CN", "CBF")

# lists to store objects and runtime
out_cydar_main <- runtime_cydar_main <- vector("list", length(thresholds))
names(out_cydar_main) <- names(runtime_cydar_main) <- thresholds




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
  
  # load data as 'ncdfFlowSet' object (for cydar)
  
  d_input <- read.ncdfFlowSet(files_load, transformation = FALSE, truncate_max_range = FALSE)
  
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
  
  
  
  
  #####################################
  # Additional pre-processing for cydar
  #####################################
  
  # --------------
  # markers to use
  # --------------
  
  cols_to_use <- cols_lineage
  
  
  
  
  ################
  # cydar pipeline
  ################
  
  # following steps in Bioconductor vignette
  
  
  # -----------------------------
  # Pre-processing of intensities
  # -----------------------------
  
  runtime_pre <- system.time({
    
    # Pooling cells together
    # note: keep all cells
    
    pool.ff <- poolCells(d_input, equalize = FALSE)
    
    
    # Estimating transformation parameters
    # note: we use 'asinh' transform with 'cofactor' = 5 (see Bendall et al. 2011, Supp. Fig. S2), 
    # instead of 'logicle' as shown in Bioconductor vignette
    
    cofactor <- 5
    
    transf_exprs <- asinh(exprs(pool.ff)[, cols_markers] / cofactor)
    
    proc.ff <- pool.ff
    exprs(proc.ff)[, cols_markers] <- transf_exprs
    
    
    # Gating out uninteresting cells
    # note: not required for this data set
    
    # keep only markers of interest
    processed.exprs <- proc.ff[, cols_to_use]
    
    
    # Normalizing intensities across batches
    # note: not required for this data set
    
  })
  
  
  # --------------------------------
  # Counting cells into hyperspheres
  # --------------------------------
  
  runtime_count <- system.time({
    
    # note: data object must be list (or ncdfFlowSet object), with one list item per sample
    
    # convert processed data to named list (this step is not shown in Bioconductor vignette)
    
    # check sample order
    sampleNames(d_input)
    sample_id
    
    n_cells <- as.vector(fsApply(d_input, nrow))
    sample_id_rep <- rep(sample_id, n_cells)
    sample_id_rep <- factor(sample_id_rep, levels = unique(sample_id_rep))
    
    d_list <- split.data.frame(exprs(processed.exprs), sample_id_rep)
    
    
    set.seed(1234)
    
    # prepare 'CyData' object
    cd <- prepareCellData(d_list)
    
    # assign cells to hyperspheres
    cd <- countCells(cd, tol = 0.5)
    
    # check hypersphere counts
    head(assay(cd))
    dim(assay(cd))
    
    # check hypersphere positions (median intensities)
    head(intensities(cd))
    
  })
  
  
  # --------------------------------------------------------------------------
  # Test for differentially abundant (DA) hyperspheres and control spatial FDR
  # --------------------------------------------------------------------------
  
  # using 'edgeR' for testing
  
  runtime_test <- system.time({
    
    # create object for edgeR
    y <- DGEList(assay(cd), lib.size = cd$totals)
    
    # filtering
    keep <- aveLogCPM(y) >= aveLogCPM(1, mean(cd$totals))
    cd <- cd[keep, ]
    y <- y[keep, ]
    
    # design matrix
    # note: including 'patient_id' for paired design
    # note: design matrix specification with intercept term
    design <- model.matrix(~ group_id + patient_id)
    
    # estimate dispersions and fit models
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust = TRUE)
    
  })
  
  
  # set up contrasts and test separately for each condition
  
  out_cydar_main[[th]] <- runtime_cydar_main[[th]] <- vector("list", length(cond_names))
  names(out_cydar_main[[th]]) <- names(runtime_cydar_main[[th]]) <- cond_names
  
  for (j in 1:length(cond_names)) {
    
    runtime_j <- system.time({
      
      # set up contrast
      contr_string <- paste0("group_id", cond_names[j])
      contrast <- makeContrasts(contr_string, levels = design)
      
      # calculate differential tests
      res_cydar <- glmQLFTest(fit, contrast = contrast)
      
      # raw p-values
      p_vals <- res_cydar$table$PValue
      
      # controlling the spatial false discovery rate (FDR)
      q_vals <- spatialFDR(intensities(cd), res_cydar$table$PValue)
      
    })
    
    # significant hyperspheres
    is.sig <- q_vals <= 0.1
    print(table(is.sig))
    
    
    # ---------------------------------------------------------------------
    # Visualizing and interpreting the results (from Bioconductor vignette)
    # ---------------------------------------------------------------------
    
    # 2D PCA plots
    if (sum(is.sig) > 0) {
      pdf(file.path(DIR_CYDAR_FILES, thresholds[th], cond_names[j], "cydar_populations_PCA_main.pdf"))
      
      sig.coords <- intensities(cd)[is.sig, ]
      sig.res <- res_cydar$table[is.sig, ]
      coords <- prcomp(sig.coords)
      
      plotCellLogFC(coords$x[, 1], coords$x[, 2], sig.res$logFC)
      
      dev.off()
    }
    
    
    # 2D tSNE plots (does not work for all thresholds/conditions)
    # if (sum(is.sig) > 0) {
    #   pdf(file.path(DIR_CYDAR_FILES, thresholds[th], cond_names[j], "cydar_populations_tSNE_main.pdf"))
    # 
    #   sig.coords <- intensities(cd)[is.sig, ]
    #   sig.res <- res_cydar$table[is.sig, ]
    # 
    #   set.seed(123)
    #   out_tsne <- Rtsne(sig.coords, pca = FALSE, verbose = TRUE)
    #   tsne_coords <- as.data.frame(out_tsne$Y)
    #   colnames(tsne_coords) <- c("tSNE_1", "tSNE_2")
    # 
    #   plotCellLogFC(tsne_coords[, 1], tsne_coords[, 2], sig.res$logFC)
    # 
    #   dev.off()
    # }
    
    
    # Median marker intensities of each hypersphere
    if (sum(is.sig) > 0) {
      pdf(file.path(DIR_CYDAR_FILES, thresholds[th], cond_names[j], "cydar_medians_main.pdf"), width = 9, height = 5.5)
      
      par(mfrow = c(3, 6), mar = c(2.1, 1.1, 3.1, 1.1))
      limits <- intensityRanges(cd, p = 0.05)
      all.markers <- rownames(markerData(cd))
      for (i in order(all.markers)) { 
        plotCellIntensity(coords$x[, 1], coords$x[, 2], sig.coords[, i], 
                          irange = limits[, i], main = all.markers[i])
      }
      par(mfrow = c(1, 1))
      
      dev.off()
    }
    
    
    # -------
    # Runtime
    # -------
    
    runtime_total <- runtime_pre[["elapsed"]] + runtime_count[["elapsed"]] + runtime_test[["elapsed"]] + runtime_j[["elapsed"]]
    print(runtime_total)
    
    runtime_cydar_main[[th]][[j]] <- runtime_total
    
    
    
    
    ##############################
    # Return results at cell level
    ##############################
    
    # Note: cydar returns q-values at the hypersphere level. Since hyperspheres overlap,
    # the q-values are not unique for each cell. To evaluate performance at the cell
    # level, we assign a unique q-value to each cell, by selecting the smallest q-value
    # for any hypersphere containing that cell.
    
    
    # number of cells per sample (including spike-in cells)
    n_cells <- as.vector(fsApply(d_input, nrow))
    
    # identify true spike-in cells (from 'spike' condition)
    is_spikein <- unlist(fsApply(d_input, function(d) exprs(d)[, "spikein"], simplify = FALSE))
    stopifnot(length(is_spikein) == sum(n_cells))
    
    # select samples for this condition and healthy
    ix_keep_cnd <- group_id %in% c("healthy", cond_names[j])
    
    
    # get smallest q-value for each cell, across all hyperspheres
    
    stopifnot(length(p_vals) == length(q_vals))
    
    # cell assignments
    cells <- cellAssignments(cd)
    # 'unpack' indices (e.g. '142183 -142188' means all values from 142183 to 142188)
    cells <- unpackIndices(cells)
    
    length(cells)
    stopifnot(length(cells) == length(p_vals))
    
    # repeat p/q-values for each hypersphere
    cells_rep <- unlist(cells)
    p_vals_rep <- rep(p_vals, sapply(cells, length))
    q_vals_rep <- rep(q_vals, sapply(cells, length))
    stopifnot(length(p_vals_rep) == length(cells_rep), 
              length(q_vals_rep) == length(cells_rep))
    # split by cell indices
    p_vals_rep_split <- split(p_vals_rep, cells_rep)
    q_vals_rep_split <- split(q_vals_rep, cells_rep)
    # get minimum p/q-value for each unique cell
    cells_p_vals <- sapply(p_vals_rep_split, function(v) min(v))
    cells_q_vals <- sapply(q_vals_rep_split, function(v) min(v))
    
    # fill in NAs for missing cells
    p_vals_all <- q_vals_all <- rep(NA, sum(n_cells))
    names(p_vals_all) <- names(q_vals_all) <- 1:sum(n_cells)
    p_vals_all[names(cells_p_vals)] <- cells_p_vals
    q_vals_all[names(cells_q_vals)] <- cells_q_vals
    
    stopifnot(length(p_vals_all) == length(q_vals_all))
    
    
    # set up data frame with results and true spike-in status at cell level
    
    which_cnd <- rep(ix_keep_cnd, n_cells)
    is_spikein_cnd <- is_spikein[which_cnd]
    stopifnot(length(p_vals_all) == length(which_cnd), 
              length(p_vals_all[which_cnd]) == length(is_spikein_cnd))
    
    res_p_vals <- p_vals_all[which_cnd]
    res_p_adj <- q_vals_all[which_cnd]
    
    # replace any NAs to ensure same set of cells is returned for all methods
    res_p_vals[is.na(res_p_vals)] <- 1
    res_p_adj[is.na(res_p_adj)] <- 1
    
    stopifnot(length(res_p_vals) == length(is_spikein_cnd), 
              length(res_p_adj) == length(is_spikein_cnd))
    
    # return values for this condition and healthy
    res <- data.frame(p_vals = res_p_vals, 
                      q_vals = res_p_adj, 
                      spikein = is_spikein_cnd)
    
    # store results
    out_cydar_main[[th]][[j]] <- res
    
  }
}




#####################
# Save output objects
#####################

save(out_cydar_main, runtime_cydar_main, 
     file = file.path(DIR_RDATA, "outputs_AML_sim_cydar_main.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_AML_sim_cydar_main.txt"))
sessionInfo()
sink()



