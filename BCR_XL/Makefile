##################################################
# Makefile to run all analyses for BCR-XL data set
##################################################


DIR_RUN_METHODS    = run_methods
DIR_GENERATE_PLOTS = generate_plots

DIR_RDATA          = ../../RData/BCR_XL
DIR_PLOTS          = ../../plots/BCR_XL


TIMESTAMP_PLOTS_MAIN_PHENOTYPES = $(DIR_PLOTS)/main_phenotypes/timestamp.txt
TIMESTAMP_PLOTS_MAIN_PROPORTIONS = $(DIR_PLOTS)/main_proportions/timestamp.txt

SCRIPTS_PLOTS_MAIN_PHENOTYPES = $(wildcard $(DIR_GENERATE_PLOTS)/main_phenotypes/*.R)
SCRIPTS_PLOTS_MAIN_PROPORTIONS = $(wildcard $(DIR_GENERATE_PLOTS)/main_proportions/*.R)

FILES_RDATA_MAIN = $(wildcard $(DIR_RDATA)/main/*.RData)




# ------------
# Default rule
# ------------

.PHONY: all
all: timestamps_plots

.PHONY: timestamps_plots
timestamps_plots: $(TIMESTAMP_PLOTS_MAIN_PHENOTYPES) \
                  $(TIMESTAMP_PLOTS_MAIN_PROPORTIONS)



# --------------
# Generate plots
# --------------

# main: phenotypes
$(TIMESTAMP_PLOTS_MAIN_PHENOTYPES): $(SCRIPTS_PLOTS_MAIN_PHENOTYPES) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_phenotypes && bash run_all.sh

# main: proportions
$(TIMESTAMP_PLOTS_MAIN_PROPORTIONS): $(SCRIPTS_PLOTS_MAIN_PROPORTIONS) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_proportions && bash run_all.sh



# -----------
# Run methods
# -----------

# note: can't use as default rule if target contains '%'

# main
$(DIR_RDATA)/main/outputs_%.RData: $(DIR_RUN_METHODS)/main/run_%.R
	cd $(DIR_RUN_METHODS)/main && Rscript $(notdir $<)



