#################################
# AML-SPIKE-IN DATA: MAIN RESULTS
#################################

# variables
DIR_BENCHMARK_DATA = prepare_data
DIR_RUN_SCRIPTS = run_methods/AML_spike_in/main
DIR_RDATA = ../RData/AML_spike_in/main
DIR_PLOT_SCRIPTS = generate_plots/AML_spike_in/main
DIR_PLOTS = ../plots/AML_spike_in/all_methods/main

SCRIPT_BENCHMARK_DATA = $(DIR_BENCHMARK_DATA)/prepare_data_AML_spike_in.R
SCRIPTS_PLOT = $(wildcard $(DIR_PLOT_SCRIPTS)/*.R)

FILES_BENCHMARK_DATA = $(wildcard ../../benchmark_data/AML_spike_in/data/*/*.fcs)
FILES_RDATA = $(wildcard $(DIR_RDATA)/*.RData)
FILES_PLOTS = $(shell find $(DIR_PLOTS) -type f -name '*.pdf')


# default rule
.PHONY: all
all: plots


# generate plots
.PHONY: plots
.SECONDARY: plots.sec

plots: $(FILES_PLOTS)
$(FILES_PLOTS): plots.sec
plots.sec: $(SCRIPTS_PLOT) $(FILES_RDATA)
	( cd $(DIR_PLOT_SCRIPTS) && Rscript $(notdir $<) )


# run methods
.PHONY: RData_files

RData_files: $(FILES_RDATA)
$(DIR_RDATA)/outputs_%.RData: $(DIR_RUN_SCRIPTS)/run_%.R $(SCRIPT_BENCHMARK_DATA) $(FILES_BENCHMARK_DATA)
	( cd $(DIR_RUN_SCRIPTS) && Rscript run_$*.R )


# generate benchmark data
.PHONY: benchmark_data
.SECONDARY: benchmark_data.sec

benchmark_data: $(FILES_BENCHMARK_DATA)
$(FILES_BENCHMARK_DATA): benchmark_data.sec
benchmark_data.sec: $(SCRIPT_BENCHMARK_DATA)
	( cd $(DIR_BENCHMARK_DATA) && Rscript $(notdir $<) )



###############
# MISCELLANEOUS
###############

# remove auto-generated files for Citrus
.PHONY: clean_Citrus
clean_Citrus:
	find ../Citrus_files -type f -delete

# remove auto-generated files for CellCnn
.PHONY: clean_CellCnn
clean_CellCnn:
	find ../CellCnn_files -type f -delete


# show variables
.PHONY: variables
variables:
	@echo SCRIPTS_PLOT: $(SCRIPTS_PLOT)


