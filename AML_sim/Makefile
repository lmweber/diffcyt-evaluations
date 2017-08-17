###################################################
# Makefile to run all analyses for AML-sim data set
###################################################


DIR_BENCHMARK = ../../../benchmark_data/AML_sim
DIR_PREPARE_DATA = 1_prepare_data

TIMESTAMP_BENCHMARK = $(DIR_BENCHMARK)/data/timestamp.txt



# ------------
# Default rule
# ------------

#.PHONY: all
#all: $(TIMESTAMP_PLOTS)



# --------------
# Generate plots
# --------------



# -----------
# Run methods
# -----------



# --------------------------
# Prepare benchmark data set
# --------------------------

# generate simulated data files (including 'nosim')
$(TIMESTAMP_BENCHMARK): $(DIR_PREPARE_DATA)/prepare_data_AML_sim.R
	cd $(DIR_PREPARE_DATA) && Rscript $(notdir $<)


