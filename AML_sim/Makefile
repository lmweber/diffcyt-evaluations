###################################################
# Makefile to run all analyses for AML-sim data set
###################################################


DIR_BENCHMARK = ../../../benchmark_data/AML_sim
DIR_PREPARE_DATA = 1_prepare_data
DIR_RUN_METHODS = 2_run_methods
DIR_GENERATE_PLOTS = 3_generate_plots
DIR_RDATA = ../../RData/AML_sim
DIR_PLOTS = ../../plots/AML_sim

TIMESTAMP_BENCHMARK = $(DIR_BENCHMARK)/data/timestamp.txt

TIMESTAMP_PLOTS = $(DIR_PLOTS)/timestamp.txt

SCRIPTS_PLOTS = $(shell find $(DIR_GENERATE_PLOTS) -type f -name '*.R')

FILES_RDATA = $(shell find $(DIR_RDATA) -type f -name '*.RData')



# ------------
# Default rule
# ------------

.PHONY: all
all: $(TIMESTAMP_PLOTS)



# --------------
# Generate plots
# --------------

# note: using automatic variables for dependencies, but typing in command for each script 
# manually (otherwise rules get too complex)
$(TIMESTAMP_PLOTS): $(SCRIPTS_PLOTS) $(FILES_RDATA)
	cd $(DIR_GENERATE_PLOTS)/main && Rscript plots_AML_sim_all_methods_main_ROC_and_TPR_FDR.R
	cd $(DIR_GENERATE_PLOTS)/supp_resolution && Rscript plots_AML_sim_all_methods_supp_resolution_ROC_and_TPR_FDR.R



# -----------
# Run methods
# -----------

# note: can't use as default rule if rule contains '%'

# main results
$(DIR_RDATA)/main/outputs_%.RData: $(DIR_RUN_METHODS)/main/run_%.R $(TIMESTAMP_BENCHMARK)
	cd $(DIR_RUN_METHODS)/main && Rscript $(notdir $<)

# supplementary: clustering resolution
$(DIR_RDATA)/supp_resolution/outputs_%.RData: $(DIR_RUN_METHODS)/supp_resolution/run_%.R $(TIMESTAMP_BENCHMARK)
	cd $(DIR_RUN_METHODS)/supp_resolution && Rscript $(notdir $<)



# --------------------------
# Prepare benchmark data set
# --------------------------

# generate simulated data files (including 'nosim')
$(TIMESTAMP_BENCHMARK): $(DIR_PREPARE_DATA)/prepare_data_AML_sim.R
	cd $(DIR_PREPARE_DATA) && Rscript $(notdir $<)



# -------------
# Miscellaneous
# -------------

# remove auto-generated files for Citrus
.PHONY: clean_Citrus
clean_Citrus:
	find ../../Citrus_files -type f -delete

# remove auto-generated files for CellCnn
.PHONY: clean_CellCnn
clean_CellCnn:
	find ../../CellCnn_files -type f -delete


