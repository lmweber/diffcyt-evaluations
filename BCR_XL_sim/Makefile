######################################################
# Makefile to run all analyses for BCR-XL-sim data set
######################################################


DIR_BENCHMARK = ../../../benchmark_data/BCR_XL_sim
DIR_PREPARE_DATA = 1_prepare_data
DIR_RUN_METHODS = 2_run_methods
DIR_GENERATE_PLOTS = 3_generate_plots
DIR_RDATA = ../../RData/BCR_XL_sim
DIR_PLOTS = ../../plots/BCR_XL_sim

TIMESTAMP_LABELS = $(DIR_BENCHMARK)/population_IDs/timestamp.txt
TIMESTAMP_BENCHMARK = $(DIR_BENCHMARK)/data/timestamp.txt
TIMESTAMP_PLOTS = $(DIR_PLOTS)/timestamp.txt

SCRIPTS_PLOTS = $(wildcard $(DIR_GENERATE_PLOTS)/*.R)

FILES_RDATA = $(DIR_RDATA)/*.RData



# ------------
# Default rule
# ------------

#.PHONY: all
#all: $(TIMESTAMP_PLOTS)



# --------------
# Generate plots
# --------------

# note: include command for each script manually (otherwise dependencies get too complex)
#$(TIMESTAMP_PLOTS): $(SCRIPTS_PLOTS) $(FILES_RDATA)
#	cd $(DIR_GENERATE_PLOTS) && Rscript plots_BCR_XL_sim_all_methods_ROC_and_TPR_FDR.R
#	cd $(DIR_GENERATE_PLOTS) && rm Rplots.pdf



# -----------
# Run methods
# -----------

# note: can't use as default rule if rule contains '%'
#$(DIR_RDATA)/outputs_%.RData: $(DIR_RUN_METHODS)/run_%.R $(TIMESTAMP_BENCHMARK)
#	cd $(DIR_RUN_METHODS) && Rscript $(notdir $<)



# --------------------------
# Prepare benchmark data set
# --------------------------

# generate simulated data files (including 'nosim')
#$(TIMESTAMP_BENCHMARK): $(DIR_PREPARE_DATA)/prepare_data_BCR_XL_sim.R $(TIMESTAMP_LABELS)
#	cd $(DIR_PREPARE_DATA) && Rscript $(notdir $<)

# reproduce cell population labels from Nowicka et al. (2017); remove temporary files
#$(TIMESTAMP_LABELS): $(DIR_PREPARE_DATA)/cell_population_labels_BCR_XL.R
#	cd $(DIR_PREPARE_DATA) && Rscript $(notdir $<)
#	cd $(DIR_PREPARE_DATA) && rm *.fcs *.xlsx *.pdf *.zip && rm consensus_plots/*.png && rmdir consensus_plots


