##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: biological coefficient of variation (BCV) plots
# 
# - BCV plots for main benchmark data
# 
# Lukas Weber, April 2018
##########################################################################################


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)
library(edgeR)


DIR_BENCHMARK <- "../../../../../benchmark_data/AML_sim/data/main"


# load saved results
DIR_RDATA <- "../../../../RData/AML_sim/main"

load(file.path(DIR_RDATA, "outputs_AML_sim_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA, "out_clusters_AML_sim_diffcyt_DA_edgeR_main.RData"))
load(file.path(DIR_RDATA, "out_objects_AML_sim_diffcyt_DA_edgeR_main.RData"))


# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/data_BCV"




################
# Generate plots
################

# loop over thresholds (th)

# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc")




for (th in 1:length(thresholds)) {
  
  # -------------
  # load raw data
  # -------------
  
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
  
  
  
  # ----------------------------------------------
  # identify true spike-in status for each cluster
  # ----------------------------------------------
  
  # get objects
  d_se <- out_objects_diffcyt_DA_edgeR_main[[th]]$d_se
  d_counts <- out_objects_diffcyt_DA_edgeR_main[[th]]$d_counts
  # using results from one condition (CN) only, since clustering is same for both conditions
  res_clusters <- out_clusters_diffcyt_DA_edgeR_main[[th]][[1]]
  
  
  # spike-in status for each cell
  is_spikein <- unlist(sapply(d_input, function(d) exprs(d)[, "spikein"]))
  stopifnot(length(is_spikein) == sum(sapply(d_input, nrow)))
  
  
  # match cells to clusters
  ix_match <- match(rowData(d_se)$cluster_id, res_clusters$cluster_id)
  
  
  # calculate number of true spike-in cells per cluster
  stopifnot(length(is_spikein) == length(ix_match))
  
  n_spikein <- sapply(split(is_spikein, ix_match), sum)
  
  
  # total number of cells per cluster
  n_cells <- rowData(d_counts)$n_cells
  
  
  # clusters containing true spike-in cells
  stopifnot(length(n_spikein) == nlevels(rowData(d_se)$cluster_id), 
            length(n_cells) == nlevels(rowData(d_se)$cluster_id))
  
  prop_spikein <- n_spikein / n_cells
  spikein <- prop_spikein > 0.1
  
  
  
  # -------------
  # design matrix
  # -------------
  
  # set up design matrix
  # note: include fixed effects for 'patient_id'
  design <- createDesignMatrix(experiment_info, cols_design = 1:2)
  design
  
  
  
  # ---------------------
  # calculate dispersions
  # ---------------------
  
  y <- DGEList(assay(d_counts))
  
  # use GLM functions due to complex design
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  
  
  
  # --------
  # BCV plot
  # --------
  
  # with points to show clusters containing true spike-in cells
  
  fn <- file.path(DIR_PLOTS, paste0("AML_sim_BCV_", thresholds[th], ".pdf"))
  
  
  pdf(fn, width = 6, height = 6)
  
  plotBCV(y)
  points(y$AveLogCPM[spikein], sqrt(y$tagwise.dispersion)[spikein], pch = 20, col = 5)
  legend("top", ">10% spike-in", pch = 20, col = 5)
  
  dev.off()
  
}




###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



