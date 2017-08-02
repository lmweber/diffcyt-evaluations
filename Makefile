# define variables

# benchmark data
BENCHMARK_DATA_FILES_AML = $(wildcard ../../benchmark_data/AML_spike_in/data/*/*.fcs)

# scripts
DIR_RUN_SCRIPTS_AML_MAIN = run_methods/AML_spike_in/main
DIR_PLOT_SCRIPTS_AML_MAIN = generate_plots/AML_spike_in
# RData files
DIR_RDATA_AML_MAIN = ../RData/AML_spike_in/main
RDATA_FILES_AML_MAIN = $(wildcard $(DIR_RDATA_AML_MAIN)/*.RData)
# plots
DIR_PLOTS_AML_MAIN = ../plots/AML_spike_in/all_methods/main
PLOTS_AML_MAIN = $(shell find $(DIR_PLOTS_AML_MAIN) -type f -name '*.pdf')


# all
.PHONY : all
all : $(RDATA_FILES_AML_MAIN)


# generate plots: ROC and TPR-FDR curves
#.PHONY : plots_ROC_TPRFDR
#plots_ROC_TPRFDR : $(PLOTS_AML_MAIN)

## not working yet
#$(PLOTS_AML_MAIN) : $(DIR_PLOT_SCRIPTS_AML_MAIN)/%.R $(RDATA_FILES_AML_MAIN)
#	( cd $(DIR_PLOT_SCRIPTS_AML_MAIN) && Rscript $*.R )


# run methods to generate RData output files
# note: run R scripts in directory where they are saved (using 'cd' and parentheses for sub-shell)
.PHONY : RData_files
RData_files : $(RDATA_FILES_AML_MAIN)

$(DIR_RDATA_AML_MAIN)/outputs_%.RData : $(DIR_RUN_SCRIPTS_AML_MAIN)/run_%.R $(BENCHMARK_DATA_FILES_AML)
	( cd $(DIR_RUN_SCRIPTS_AML_MAIN) && Rscript run_$*.R )


# generate benchmark data set
.PHONY : benchmark_data_AML
benchmark_data_AML : prepare_data/prepare_data_AML_spike_in.R
	( cd prepare_data && Rscript prepare_data_AML_spike_in.R )

$(BENCHMARK_DATA_FILES_AML) : prepare_data/prepare_data_AML_spike_in.R


# remove auto-generated files for Citrus
.PHONY : clean_Citrus
clean_Citrus :
	find ../Citrus_files -type f -delete

# remove auto-generated files for CellCnn
.PHONY : clean_CellCnn
clean_CellCnn :
	find ../CellCnn_files -type f -delete


# show variables
.PHONY : variables
variables :
	@echo DIR_PLOT_SCRIPTS_AML_MAIN: $(DIR_PLOT_SCRIPTS_AML_MAIN)

