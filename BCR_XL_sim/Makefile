######################################################
# Makefile to run all analyses for BCR-XL-sim data set
######################################################

# ----------------
# Define variables
# ----------------

DIR_BENCHMARK = ../../../benchmark_data/BCR_XL_sim



# ------------
# Default rule
# ------------

#.PHONY: all
#all:



# --------------------------
# Prepare benchmark data set
# --------------------------

# 'nosim': copy .fcs files without modification
$(wildcard $(DIR_BENCHMARK)/data/*/*.fcs): 1_prepare_data/prepare_data_BCR_XL_sim.R
	cd 1_prepare_data && Rscript prepare_data_BCR_XL_sim.R

# reproduce cell population labels from Nowicka et al. (2017)
$(DIR_BENCHMARK)/population_IDs/BCR_XL_population_IDs.csv: 1_prepare_data/cell_population_labels_BCR_XL.R
	cd 1_prepare_data && Rscript $(notdir $<)
	cd 1_prepare_data && rm *.fcs *.xlsx *.pdf *.zip && rm consensus_plots/*.png && rmdir consensus_plots

