#####################################################
# Makefile to run all analyses for Anti-PD-1 data set
#####################################################


DIR_RUN_METHODS    = run_methods
DIR_GENERATE_PLOTS = generate_plots

DIR_RDATA          = ../../RData/Anti_PD_1
DIR_PLOTS          = ../../plots/Anti_PD_1


TIMESTAMP_PLOTS_MAIN          = $(DIR_PLOTS)/main/timestamp.txt
TIMESTAMP_PLOTS_SUPPLEMENTARY = $(DIR_PLOTS)/supplementary/timestamp.txt

SCRIPTS_PLOTS_MAIN          = $(wildcard $(DIR_GENERATE_PLOTS)/main/*.R)
SCRIPTS_PLOTS_SUPPLEMENTARY = $(wildcard $(DIR_GENERATE_PLOTS)/supplementary/*.R)

FILES_RDATA_MAIN          = $(wildcard $(DIR_RDATA)/main/*.RData)
FILES_RDATA_SUPPLEMENTARY = $(wildcard $(DIR_RDATA)/supplementary/*.RData)




# ------------
# Default rule
# ------------

.PHONY: all
all: timestamps_plots

.PHONY: timestamps_plots
timestamps_plots: $(TIMESTAMP_PLOTS_MAIN)
                  #$(TIMESTAMP_PLOTS_SUPPLEMENTARY)



# --------------
# Generate plots
# --------------

# main
$(TIMESTAMP_PLOTS_MAIN): $(SCRIPTS_PLOTS_MAIN) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main && bash run_all.sh

# supplementary
$(TIMESTAMP_PLOTS_SUPPLEMENTARY): $(SCRIPTS_PLOTS_SUPPLEMENTARY) $(FILES_RDATA_SUPPLEMENTARY)
	cd $(DIR_GENERATE_PLOTS)/supplementary && bash run_all.sh




# -----------
# Run methods
# -----------

# note: can't use as default rule if target contains '%'

# main
$(DIR_RDATA)/main/outputs_%.RData: $(DIR_RUN_METHODS)/main/run_%.R
	cd $(DIR_RUN_METHODS)/main && Rscript $(notdir $<)

# supplementary
$(DIR_RDATA)/supplementary/outputs_%.RData: $(DIR_RUN_METHODS)/supplementary/run_%.R
	cd $(DIR_RUN_METHODS)/supplementary && Rscript $(notdir $<)



