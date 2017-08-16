######################################################
# Makefile to run all analyses for BCR-XL-sim data set
######################################################

# variables

DIR_BENCHMARK = ../../../benchmark_data/BCR_XL_sim

DIR_PREPARE_DATA = 1_prepare_data
DIR_RUN_METHODS = 2_run_methods
DIR_GENERATE_PLOTS = 3_generate_plots

DIR_RDATA = ../../RData/BCR_XL_sim
DIR_PLOTS = ../../plots/BCR_XL_sim


FILES_BENCHMARK = $(wildcard $(DIR_BENCHMARK)/data/*/*.fcs)
SCRIPT_BENCHMARK = $(DIR_PREPARE_DATA)/prepare_data_BCR_XL_sim.R

FILES_LABELS = $(wildcard $(DIR_BENCHMARK)/population_IDs/*.csv)
SCRIPT_LABELS = $(DIR_PREPARE_DATA)/cell_population_labels_BCR_XL.R

FILES_RDATA = $(DIR_RDATA)/*.RData

FILES_PLOTS = $(shell find $(DIR_PLOTS) -type f -name '*.pdf')
SCRIPTS_PLOTS = $(wildcard $(DIR_GENERATE_PLOTS)/*.R)



# ------------
# Default rule
# ------------

.PHONY: all
all: plot_files

.PHONY: plot_files
plot_files: $(FILES_PLOTS)



# --------------
# Generate plots
# --------------

# use 'secondary' variable so command is only run once
.SECONDARY: plot_files.sec
$(FILES_PLOTS): plot_files.sec

# note: include command for each script manually (otherwise dependencies get too complex)
plot_files.sec: $(SCRIPTS_PLOTS) $(FILES_RDATA)
	cd $(DIR_GENERATE_PLOTS) && Rscript plots_BCR_XL_sim_all_methods_ROC_and_TPR_FDR.R



# -----------
# Run methods
# -----------

$(DIR_RDATA)/outputs_%.RData: $(DIR_RUN_METHODS)/run_%.R $(FILES_BENCHMARK)
	cd $(DIR_RUN_METHODS) && Rscript $(notdir $<)



# --------------------------
# Prepare benchmark data set
# --------------------------

# generate simulated data files (including 'nosim')
$(FILES_BENCHMARK): $(SCRIPT_BENCHMARK) $(FILES_LABELS)
	cd $(DIR_PREPARE_DATA) && Rscript $(notdir $<)

# reproduce cell population labels from Nowicka et al. (2017); remove temporary files
$(FILES_LABELS): $(SCRIPT_LABELS)
	cd $(DIR_PREPARE_DATA) && Rscript $(notdir $<)
	cd $(DIR_PREPARE_DATA) && rm *.fcs *.xlsx *.pdf *.zip && rm consensus_plots/*.png && rmdir consensus_plots


