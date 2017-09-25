##########################################################################################
# Generate plots
# 
# - data set: AML-sim
# - plot type: ROC curves
# - method: all methods
# 
# - supplementary results: randomized benchmark data sets
# 
# Lukas Weber, September 2017
##########################################################################################



# path to save plots
DIR_PLOTS <- "../../../../plots/AML_sim/supp_benchmark_random"






###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_PLOTS, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



