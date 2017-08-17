###################################################
# Makefile to run all analyses for AML-sim data set
###################################################


DIR_BENCHMARK      = ../../../benchmark_data/AML_sim
DIR_PREPARE_DATA   = 1_prepare_data
DIR_RUN_METHODS    = 2_run_methods
DIR_GENERATE_PLOTS = 3_generate_plots
DIR_RDATA          = ../../RData/AML_sim
DIR_PLOTS          = ../../plots/AML_sim

TIMESTAMP_BENCHMARK             = $(DIR_BENCHMARK)/data/timestamp.txt
TIMESTAMP_PLOTS_MAIN            = $(DIR_PLOTS)/main/timestamp.txt
TIMESTAMP_PLOTS_SUPP_RESOLUTION = $(DIR_PLOTS)/supp_resolution/timestamp.txt

SCRIPTS_PLOTS_MAIN            = $(wildcard $(DIR_GENERATE_PLOTS)/main/*.R)
SCRIPTS_PLOTS_SUPP_RESOLUTION = $(wildcard $(DIR_GENERATE_PLOTS)/supp_resolution/*.R)

FILES_RDATA_MAIN            = $(wildcard $(DIR_RDATA)/main/*.R)
FILES_RDATA_SUPP_RESOLUTION = $(wildcard $(DIR_RDATA)/supp_resolution/*.R)



# ------------
# Default rule
# ------------

.PHONY: all
all: timestamps_plots

.PHONY: timestamps_plots
timestamps_plots: $(TIMESTAMP_PLOTS_MAIN) $(TIMESTAMP_PLOTS_SUPP_RESOLUTION)



# --------------
# Generate plots
# --------------

# note: separate rule for each script (otherwise dependencies too complex)

$(TIMESTAMP_PLOTS_MAIN): $(SCRIPTS_PLOTS_MAIN) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main && Rscript $(notdir $<)

$(TIMESTAMP_PLOTS_SUPP_RESOLUTION): $(SCRIPTS_PLOTS_SUPP_RESOLUTION) $(FILES_RDATA_SUPP_RESOLUTION)
	cd $(DIR_GENERATE_PLOTS)/supp_resolution && Rscript $(notdir $<)



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


