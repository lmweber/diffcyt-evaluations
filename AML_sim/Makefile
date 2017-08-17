###################################################
# Makefile to run all analyses for AML-sim data set
###################################################


DIR_BENCHMARK = ../../../benchmark_data/AML_sim
DIR_PREPARE_DATA = 1_prepare_data
DIR_RUN_METHODS = 2_run_methods
DIR_GENERATE_PLOTS = 3_generate_plots
DIR_RDATA = ../../RData/AML_sim

TIMESTAMP_BENCHMARK = $(DIR_BENCHMARK)/data/timestamp.txt

FILES_RDATA_MAIN = $(DIR_RDATA)/main/*.RData



# ------------
# Default rule
# ------------

.PHONY: all
all: $(FILES_RDATA_MAIN)



# --------------
# Generate plots
# --------------



# -----------
# Run methods
# -----------

# note: can't use as default rule if rule contains '%'

# main results
$(DIR_RDATA)/main/outputs_%.RData: $(DIR_RUN_METHODS)/main/run_%.R $(TIMESTAMP_BENCHMARK)
	cd $(DIR_RUN_METHODS)/main && Rscript $(notdir $<)


# remove auto-generated files for Citrus
.PHONY: clean_Citrus
clean_Citrus:
	find ../../Citrus_files -type f -delete

# remove auto-generated files for CellCnn
.PHONY: clean_CellCnn
clean_CellCnn:
	find ../../CellCnn_files -type f -delete



# --------------------------
# Prepare benchmark data set
# --------------------------

# generate simulated data files (including 'nosim')
$(TIMESTAMP_BENCHMARK): $(DIR_PREPARE_DATA)/prepare_data_AML_sim.R
	cd $(DIR_PREPARE_DATA) && Rscript $(notdir $<)


