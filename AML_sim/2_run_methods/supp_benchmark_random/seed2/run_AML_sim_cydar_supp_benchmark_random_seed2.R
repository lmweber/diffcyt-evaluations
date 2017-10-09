##########################################################################################
# Script to run methods
# 
# - method: cydar
# - data set: AML-sim
# 
# - supplementary results: randomized benchmark data sets (random seed 'seed2')
# 
# Lukas Weber, October 2017
##########################################################################################


library(flowCore)
library(ncdfFlow)
library(cydar)
library(edgeR)


DIR_BENCHMARK <- "../../../../../../benchmark_data/AML_sim/data/random_seeds/seed2"
DIR_RDATA <- "../../../../../RData/AML_sim/supp_benchmark_random/seed2"
DIR_SESSION_INFO <- "../../../../../session_info/AML_sim/supp_benchmark_random/seed2"




################################
# Loop to run for each threshold
################################

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")

# lists to store objects
out_cydar_supp_benchmark_random_seed2 <- runtime_cydar_supp_benchmark_random_seed2 <- vector("list", length(thresholds))
names(out_cydar_supp_benchmark_random_seed2) <- names(runtime_cydar_supp_benchmark_random_seed2) <- thresholds




for (th in 1:length(thresholds)) {
  
  ###########################
  # Load data, pre-processing
  ###########################
  
  # filenames
  files_healthy <- list.files(file.path(DIR_BENCHMARK, "healthy"), 
                              pattern = "_randomseed[0-9]+\\.fcs$", full.names = TRUE)
  files_CN <- list.files(file.path(DIR_BENCHMARK, "CN"), 
                         pattern = paste0("_", thresholds[th], "_randomseed[0-9]+\\.fcs$"), full.names = TRUE)
  files_CBF <- list.files(file.path(DIR_BENCHMARK, "CBF"), 
                          pattern = paste0("_", thresholds[th], "_randomseed[0-9]+\\.fcs$"), full.names = TRUE)
  
  # load data
  files_load <- c(files_healthy, files_CN, files_CBF)
  files_load
  
  d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)
  
  # sample IDs, group IDs, patient IDs
  sample_IDs <- gsub("(_[0-9]+pc$)|(_0\\.[0-9]+pc$)", "", 
                     gsub("^AML_sim_", "", 
                          gsub("_randomseed[0-9]+\\.fcs$", "", basename(files_load))))
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
  
  cols_clustering <- cols_lineage
  
  
  # --------------
  # transform data
  # --------------
  
  # 'asinh' transform with 'cofactor' = 5 (see Bendall et al. 2011, Supp. Fig. S2)
  
  # note: 'cydar' vignette uses biexponential transform instead (using 'estimateLogicle'
  # function from 'flowCore' package)
  
  cofactor <- 5
  
  d_input <- lapply(d_input, function(d) {
    e <- exprs(d)
    e[, cols_markers] <- asinh(e[, cols_markers] / cofactor)
    flowFrame(e)
  })
  
  
  
  
  ################
  # cydar pipeline
  ################
  
  # following steps in Bioconductor vignette
  
  
  # -----------------------------
  # Pre-processing of intensities
  # -----------------------------
  
  runtime_pre <- system.time({
    
    # Subset markers and convert input data to required format ('ncdfFlowSet' object)
    
    d_input_cydar <- lapply(d_input, function(d) {
      e <- exprs(d)
      e <- e[, cols_clustering]
      flowFrame(e)
    })
    
    names(d_input_cydar) <- sample_IDs
    
    d_input_cydar <- flowSet(d_input_cydar)
    d_input_cydar <- ncdfFlowSet(d_input_cydar)
    
    
    # Pooling cells together
    # note: skip this since then 'prepareCellData' in next section doesn't work
    
    # Estimating transformation parameters
    # note: not required here, since we already applied 'asinh' transform above
    
    # Gating out uninteresting cells
    # note: not required for this data set
    
    # Normalizing intensities across batches
    # note: not required (see vignette)
    
  })
  
  
  # --------------------------------
  # Counting cells into hyperspheres
  # --------------------------------
  
  runtime_count <- system.time({
    
    set.seed(1234)
    
    # prepare 'CyData' object
    cd <- prepareCellData(d_input_cydar)
    
    # assign cells to hyperspheres
    cd <- countCells(cd)
    
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
  # note: test separately for each condition: CN vs. healthy, CBF vs. healthy
  
  
  runtime_test <- system.time({
    
    # create object
    y <- DGEList(assay(cd), lib.size = cd$totals)
    
    # filtering
    keep <- aveLogCPM(y) >= aveLogCPM(5, mean(cd$totals))
    cd <- cd[keep, ]
    y <- y[keep, ]
    
    # design matrix
    # note: including 'patient_IDs' for paired design
    design <- model.matrix(~ 0 + group_IDs + patient_IDs)
    
    # estimate dispersions and fit models
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust = TRUE)
    
  })
  
  
  # set up contrasts and test separately for each condition
  
  out_cydar_supp_benchmark_random_seed2[[th]] <- runtime_cydar_supp_benchmark_random_seed2[[th]] <- vector("list", length(cond_names))
  names(out_cydar_supp_benchmark_random_seed2[[th]]) <- names(runtime_cydar_supp_benchmark_random_seed2[[th]]) <- cond_names
  
  for (j in 1:length(cond_names)) {
    
    runtime_j <- system.time({
      
      # set up contrast
      contr_string <- paste0("group_IDs", cond_names[j], " - group_IDshealthy")
      contrast <- makeContrasts(contr_string, levels = design)
      
      # differential testing
      res_cydar <- glmQLFTest(fit, contrast = contrast)
      
      # raw p-values
      p_vals <- res_cydar$table$PValue
      
      # controlling the spatial false discovery rate (FDR)
      q_vals <- spatialFDR(intensities(cd), res_cydar$table$PValue)
      
    })
    
    # significant hyperspheres
    is.sig <- q_vals <= 0.9
    print(summary(is.sig))
    
    # # plots
    # sig.coords <- intensities(cd)[is.sig, ]
    # sig.res <- res_cydar$table[is.sig, ]
    # coords <- prcomp(sig.coords)
    # plotCellLogFC(coords$x[, 1], coords$x[, 2], sig.res$logFC)
    # 
    # par(mfrow = c(4, 4), mar = c(2.1, 1.1, 3.1, 1.1))
    # limits <- intensityRanges(cd, p = 0.05)
    # all.markers <- rownames(markerData(cd))
    # for (i in order(all.markers)) { 
    #   plotCellIntensity(coords$x[, 1], coords$x[, 2], sig.coords[, i], 
    #                     irange=limits[, i], main=all.markers[i])
    # }
    # par(mfrow = c(1, 1))
    
    # runtime
    runtime_total <- runtime_pre[["elapsed"]] + runtime_count[["elapsed"]] + runtime_test[["elapsed"]] + runtime_j[["elapsed"]]
    print(runtime_total)
    
    runtime_cydar_supp_benchmark_random_seed2[[th]][[j]] <- runtime_total
    
    
    
    
    ##############################
    # Return results at cell level
    ##############################
    
    # Note: cydar evaluates q-values at the hypersphere level. Since hyperspheres overlap,
    # the q-values are not unique for each cell. To evaluate performance at the cell
    # level, we assign a unique q-value to each cell, by selecting the smallest q-value
    # for any hypersphere containing that cell.
    
    
    # number of cells per sample (including spike-in cells)
    n_cells <- sapply(d_input, nrow)
    
    # spike-in status for each cell
    is_spikein <- unlist(sapply(d_input, function(d) exprs(d)[, "spikein"]))
    stopifnot(length(is_spikein) == sum(n_cells))
    
    # select samples for this condition
    ix_keep_cnd <- group_IDs == cond_names[j]
    
    
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
    stopifnot(length(p_vals_rep) == length(cells_rep), length(q_vals_rep) == length(cells_rep))
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
    stopifnot(length(p_vals_all[which_cnd]) == length(is_spikein_cnd))
    
    res_p_vals <- p_vals_all[which_cnd]
    res_q_vals <- q_vals_all[which_cnd]
    
    # replace any NAs to ensure same set of cells is returned for all methods
    res_p_vals[is.na(res_p_vals)] <- 1
    res_q_vals[is.na(res_q_vals)] <- 1
    
    # return values for this condition only
    res <- data.frame(p_vals = res_p_vals, 
                      q_vals = res_q_vals, 
                      spikein = is_spikein_cnd)
    
    # store results
    out_cydar_supp_benchmark_random_seed2[[th]][[j]] <- res
    
  }
}




#####################
# Save output objects
#####################

save(out_cydar_supp_benchmark_random_seed2, runtime_cydar_supp_benchmark_random_seed2, 
     file = file.path(DIR_RDATA, "outputs_AML_sim_cydar_supp_benchmark_random_seed2.RData"))




#####################
# Session information
#####################

sink(file.path(DIR_SESSION_INFO, "session_info_AML_sim_cydar_supp_benchmark_random_seed2.txt"))
sessionInfo()
sink()



