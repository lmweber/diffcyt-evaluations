######################################################
# Makefile to run all analyses for BCR-XL-sim data set
######################################################


DIR_BENCHMARK      = ../../../benchmark_data/BCR_XL_sim/data

DIR_PREPARE_DATA   = 1_prepare_data
DIR_RUN_METHODS    = 2_run_methods
DIR_GENERATE_PLOTS = 3_generate_plots

DIR_RDATA          = ../../RData/BCR_XL_sim
DIR_PLOTS          = ../../plots/BCR_XL_sim


TIMESTAMP_BENCHMARK_MAIN            = $(DIR_BENCHMARK)/main/timestamp.txt
TIMESTAMP_BENCHMARK_NULL_SIMULATION = $(DIR_BENCHMARK)/null_simulation/timestamp.txt

TIMESTAMP_PLOTS_MAIN_PERFORMANCE = $(DIR_PLOTS)/main_performance/timestamp.txt
TIMESTAMP_PLOTS_MAIN_TSNE        = $(DIR_PLOTS)/main_tSNE/timestamp.txt
TIMESTAMP_PLOTS_MAIN_HEATMAPS    = $(DIR_PLOTS)/main_heatmaps/timestamp.txt
TIMESTAMP_PLOTS_NULL_SIMULATION  = $(DIR_PLOTS)/null_simulation/timestamp.txt

SCRIPTS_PLOTS_MAIN_PERFORMANCE = $(wildcard ($(DIR_GENERATE_PLOTS)/main_performance/*.R))
SCRIPTS_PLOTS_MAIN_TSNE        = $(wildcard ($(DIR_GENERATE_PLOTS)/main_tSNE/*.R))
SCRIPTS_PLOTS_MAIN_HEATMAPS    = $(wildcard ($(DIR_GENERATE_PLOTS)/main_heatmaps/*.R))
SCRIPTS_PLOTS_NULL_SIMULATION  = $(wildcard ($(DIR_GENERATE_PLOTS)/null_simulation/*.R))

FILES_RDATA_MAIN            = $(wildcard $(DIR_RDATA)/main/*.RData)
FILES_RDATA_NULL_SIMULATION = $(wildcard $(DIR_RDATA)/null_simulation/*.RData)




# ------------
# Default rule
# ------------

.PHONY: all
all: timestamps_plots

.PHONY: timestamps_plots
timestamps_plots: $(TIMESTAMP_PLOTS_MAIN_PERFORMANCE) \
                  $(TIMESTAMP_PLOTS_MAIN_TSNE) \
                  $(TIMESTAMP_PLOTS_MAIN_HEATMAPS) \
                  $(TIMESTAMP_PLOTS_NULL_SIMULATION)




# --------------
# Generate plots
# --------------

# note: separate rule for each set of results (otherwise dependencies get too complex)

# main results: performance metrics
$(TIMESTAMP_PLOTS_MAIN_PERFORMANCE): $(SCRIPTS_PLOTS_MAIN_PERFORMANCE) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_performance && bash run_all.sh

# main results: tSNE plots
$(TIMESTAMP_PLOTS_MAIN_TSNE): $(SCRIPTS_PLOTS_MAIN_TSNE) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_tSNE && bash run_all.sh

# main results: heatmaps
$(TIMESTAMP_PLOTS_MAIN_HEATMAPS): $(SCRIPTS_PLOTS_MAIN_HEATMAPS) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_heatmaps && bash run_all.sh

# null simulation
$(TIMESTAMP_PLOTS_NULL_SIMULATION): $(SCRIPTS_PLOTS_NULL_SIMULATION) $(FILES_RDATA_NULL_SIMULATION)
	cd $(DIR_GENERATE_PLOTS)/null_simulation && bash run_all.sh




# -----------
# Run methods
# -----------

# note: can't use as default rule if target contains '%'

# main results: performance metrics
$(DIR_RDATA)/main/outputs_%.RData: $(DIR_RUN_METHODS)/main/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/main && Rscript $(notdir $<)

# null simulation
$(DIR_RDATA)/null_simulation/outputs_%.RData: $(DIR_RUN_METHODS)/null_simulation/run_%.R $(TIMESTAMP_BENCHMARK_NULL_SIMULATION)
	cd $(DIR_RUN_METHODS)/null_simulation && Rscript $(notdir $<)




# ---------------------------
# Prepare benchmark data sets
# ---------------------------

# generate benchmark data files: main
$(TIMESTAMP_BENCHMARK_MAIN): $(DIR_PREPARE_DATA)/main/prepare_data_BCR_XL_sim_main.R
	cd $(DIR_PREPARE_DATA)/main && Rscript $(notdir $<)

# generate benchmark data files: main
$(TIMESTAMP_BENCHMARK_NULL_SIMULATION): $(DIR_PREPARE_DATA)/null_simulation/prepare_data_BCR_XL_sim_null.R
	cd $(DIR_PREPARE_DATA)/null_simulation && Rscript $(notdir $<)



