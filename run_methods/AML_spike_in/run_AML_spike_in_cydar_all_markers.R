##########################################################################################
# Script to run methods
# 
# - method: cydar-all-markers
# - data set: AML-spike-in
# 
# Lukas Weber, July 2017
##########################################################################################


library(flowCore)
library(ncdfFlow)
library(cydar)
library(edgeR)




################################
# Loop to run for each threshold
################################

# spike-in thresholds (must match filenames)
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")

# condition names
cond_names <- c("CN", "CBF")

# lists to store objects
n_cells_thresholds <- is_spikein_thresholds <- 
  out_cydar_data <- out_cydar_tests <- out_cydar_pvals <- out_cydar_qvals <- 
  vector("list", length(thresholds))

names(n_cells_thresholds) <- names(is_spikein_thresholds) <- 
  names(out_cydar_data) <- names(out_cydar_tests) <- names(out_cydar_pvals) <- names(out_cydar_qvals) <- 
  thresholds




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
  
  
  # -----------------------------------------
  # store number of cells and spike-in status
  # -----------------------------------------
  
  # number of cells
  n_cells_thresholds[[th]] <- sapply(d_input, nrow)
  
  # spike-in status for each cell
  is_spikein_thresholds[[th]] <- vector("list", length(sample_IDs))
  names(is_spikein_thresholds[[th]]) <- sample_IDs
  
  for (i in 1:length(sample_IDs)) {
    exprs_i <- exprs(d_input[[i]])
    is_spikein_thresholds[[th]][[i]] <- exprs_i[, "spikein"]
  }
  
  
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
  
  
  
  
  ##################
  # 'cydar' pipeline
  ##################
  
  # following steps in Bioconductor vignette
  
  
  # -----------------------------
  # Pre-processing of intensities
  # -----------------------------
  
  # Subset markers and convert input data to required format ('ncdfFlowSet' object)
  
  d_input <- lapply(d_input, function(d) {
    e <- exprs(d)
    e <- e[, cols_to_use]
    flowFrame(e)
  })
  
  names(d_input) <- sample_IDs
  
  d_input <- flowSet(d_input)
  d_input <- ncdfFlowSet(d_input)
  
  
  # Pooling cells together
  # note: skip this since then 'prepareCellData' in next section doesn't work
  
  # Estimating transformation parameters
  # note: not required here, since we already applied 'asinh' transform above
  
  # Gating out uninteresting cells
  # note: not required for this data set
  
  # Normalizing intensities across batches
  # note: not required (see vignette)
  
  
  # --------------------------------
  # Counting cells into hyperspheres
  # --------------------------------
  
  # prepare 'CyData' object
  cd <- prepareCellData(d_input)
  
  # assign cells to hyperspheres
  cd <- countCells(cd)
  
  # check hypersphere counts
  head(assay(cd))
  dim(assay(cd))
  
  # check hypersphere positions (median intensities)
  head(intensities(cd))
  
  
  # --------------------------------------------------------------------------
  # Test for differentially abundant (DA) hyperspheres and control spatial FDR
  # --------------------------------------------------------------------------
  
  # using 'edgeR' for testing
  # note: test separately for each condition: CN vs. healthy, CBF vs. healthy
  
  
  # create object
  y <- DGEList(assay(cd), lib.size = cd$totals)
  
  # filtering
  keep <- aveLogCPM(y) >= aveLogCPM(5, mean(cd$totals))
  cd <- cd[keep, ]
  y <- y[keep, ]
  
  # design matrix
  # note: including 'block_IDs' for paired design
  design <- model.matrix(~ 0 + group_IDs + block_IDs)
  
  # estimate dispersions and fit models
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  
  
  # set up contrasts and test separately for each condition
  out_cydar_data[[th]] <- out_cydar_tests[[th]] <- out_cydar_pvals[[th]] <- out_cydar_qvals[[th]] <- 
    vector("list", length(cond_names))
  names(out_cydar_data[[th]]) <- names(out_cydar_tests[[th]]) <- names(out_cydar_pvals[[th]]) <- names(out_cydar_qvals[[th]]) <- 
    cond_names
  
  
  for (j in 1:length(cond_names)) {
    
    # set up contrast
    contr_str <- paste0("group_IDs", cond_names[j], " - group_IDshealthy")
    contrast <- makeContrasts(contr_str, levels = design)
    
    runtime_cydar <- system.time({
      # differential testing
      res <- glmQLFTest(fit, contrast = contrast)
      
      # raw p-values
      pvals <- res$table$PValue
      
      # controlling the spatial false discovery rate (FDR)
      qvals <- spatialFDR(intensities(cd), res$table$PValue)
    })
    
    print(runtime_cydar)
    
    # significant hyperspheres
    is.sig <- qvals <= 0.9
    print(summary(is.sig))
    
    # # plots
    # sig.coords <- intensities(cd)[is.sig, ]
    # sig.res <- res$table[is.sig, ]
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
    
    # store results
    out_cydar_data[[th]][[j]] <- cd
    out_cydar_tests[[th]][[j]] <- res
    out_cydar_pvals[[th]][[j]] <- pvals
    out_cydar_qvals[[th]][[j]] <- qvals
  }
}




#####################
# Save output objects
#####################

save.image("../../../RData/AML_spike_in/outputs_AML_spike_in_cydar_all_markers.RData")




#####################
# Session information
#####################

sink("../../../session_info/AML_spike_in/session_info_AML_spike_in_cydar_all_markers.txt")
sessionInfo()
sink()



