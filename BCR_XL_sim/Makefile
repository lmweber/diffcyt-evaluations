######################################################
# Makefile to run all analyses for BCR-XL-sim data set
######################################################


# ------------
# Default rule
# ------------

#.PHONY: all
#all:



# --------------
# Generate plots
# --------------



# -----------
# Run methods
# -----------



# --------------------------
# Prepare benchmark data set
# --------------------------

DIR_BENCHMARK = ../../../benchmark_data/BCR_XL_sim
DIR_PREPARE_DATA = 1_prepare_data

FILES_BENCHMARK = $(wildcard $(DIR_BENCHMARK)/data/*/*.fcs)
SCRIPT_BENCHMARK = $(DIR_PREPARE_DATA)/prepare_data_BCR_XL_sim.R

CSV_LABELS = $(DIR_BENCHMARK)/population_IDs/BCR_XL_population_IDs.csv
SCRIPT_LABELS = $(DIR_PREPARE_DATA)/cell_population_labels_BCR_XL.R


# 'nosim': copy .fcs files without modification
$(FILES_BENCHMARK): $(SCRIPT_BENCHMARK) $(CSV_LABELS)
	cd $(DIR_PREPARE_DATA) && Rscript $(notdir $<)

# reproduce cell population labels from Nowicka et al. (2017)
$(CSV_LABELS): $(SCRIPT_LABELS)
	cd $(DIR_PREPARE_DATA) && Rscript $(notdir $<)
	cd $(DIR_PREPARE_DATA) && rm *.fcs *.xlsx *.pdf *.zip && rm consensus_plots/*.png && rmdir consensus_plots


