######################################################
# Makefile to run all analyses for BCR-XL-sim data set
######################################################

# define variables
FILES_BENCHMARK_DATA = $(wildcard ../../../benchmark_data/BCR_XL_sim/data/*/*.fcs)


# default rule
#.PHONY: all
#all:


# reproduce cell population labels from Nowicka et al. (2017)
FILES_BENCHMARK_DATA: 1_prepare_data/cell_population_labels_BCR_XL.R
	cd 1_prepare_data && Rscript $(notdir $<)
	cd 1_prepare_data && rm *.fcs *.xlsx *.pdf *.zip && rm consensus_plots/*.png && rmdir consensus_plots


