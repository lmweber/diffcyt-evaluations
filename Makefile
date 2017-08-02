# define variables

# benchmark data
BENCHMARK_DATA_AML = $(wildcard ../../benchmark_data/AML_spike_in/data/*/*.fcs)

# scripts
DIR_RUN_SCRIPTS_AML_MAIN = run_methods/AML_spike_in/main

# RData files
DIR_RDATA_AML_MAIN = ../RData/AML_spike_in/main
RDATA_FILES_AML_MAIN = $(wildcard $(DIR_RDATA_AML_MAIN)/*.RData)


#FILES_SCRIPTS_PLOTS_AML = $(wildcard $(DIR_SCRIPTS_PLOTS_AML)/*.R)
#FILES_SCRIPTS_PLOTS_AML_ALL_METHODS = $(wildcard $(DIR_SCRIPTS_PLOTS_AML)/plots_AML_spike_in_all_methods*.R)

#FILES_RDATA_AML = $(wildcard $(DIR_RDATA_AML)/*.RData)
#FILES_RDATA_AML = $(patsubst $(DIR_SCRIPTS_RUN_AML)/run_%.R, $(DIR_RDATA_AML)/outputs_%.RData, $(FILES_SCRIPTS_RUN_AML))

#FILES_PLOTS_AML_ALL_METHODS = $(wildcard $(DIR_PLOTS_AML)/all_methods/*/*/*/*.pdf)


# all
.PHONY : all
all : $(RDATA_FILES_AML_MAIN)


# generate plots for all methods combined
#.PHONY : plots_all
#plots_all : $(FILES_PLOTS_AML_ALL_METHODS)

#$(FILES_PLOTS_AML_ALL_METHODS) : $(wildcard $(DIR_SCRIPTS_PLOTS_AML)/plots_AML_spike_in_all_methods_%.R) $(FILES_RDATA_AML)
#	( cd $(DIR_SCRIPTS_PLOTS_AML) && Rscript plots_AML_spike_in_all_methods_%.R )


# run methods to generate RData output files
# note: run R scripts in directory where they are saved (using 'cd' and parentheses for sub-shell)
.PHONY : RData_files
RData_files : $(RDATA_FILES_AML_MAIN)

$(DIR_RDATA_AML_MAIN)/outputs_%.RData : $(DIR_RUN_SCRIPTS_AML_MAIN)/run_%.R $(BENCHMARK_DATA_AML)
	( cd $(DIR_RUN_SCRIPTS_AML_MAIN) && Rscript run_$*.R )


# show variables
.PHONY : variables
variables :
	@echo RDATA_FILES_AML_MAIN: $(RDATA_FILES_AML_MAIN)
