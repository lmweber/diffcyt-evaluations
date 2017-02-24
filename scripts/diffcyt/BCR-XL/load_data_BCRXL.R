####################################
# Script to load data set 'BCR-XL' #
####################################


# first script in pipeline


library(diffcyt)
library(flowCore)
library(magrittr)


# filenames
files <- list.files("../../../../benchmark_data/Citrus_paper_data/experiment_15713_files", full.names = TRUE)

# groups: BCR-XL vs. reference
files_BCRXL <- files[grep("patient[1-8]_BCR-XL\\.fcs$", files)]
files_ref <- files[grep("patient[1-8]_Reference\\.fcs$", files)]

# load data
files_load <- c(files_BCRXL, files_ref)
data_BCRXL <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)

# sample IDs
basename(files_load) %>% 
  gsub("^PBMC8_30min_|\\.fcs$", "", .) %>% 
  gsub("BCR-XL", "BCRXL", .) %>% 
  gsub("Reference", "ref", .) -> 
  sample_IDs

sample_IDs

names(data_BCRXL) <- sample_IDs

# convert to flowSet (FlowSOM requires input data as flowSet)
data_BCRXL <- as(data_BCRXL, "flowSet")

# check 
checkInputData(data_BCRXL)


