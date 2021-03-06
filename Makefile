#############################################################
# Top-level Makefile to run other Makefiles in subdirectories
#############################################################


all:
	$(MAKE) -C AML_sim
	$(MAKE) -C BCR_XL_sim
	$(MAKE) -C Anti_PD_1
	$(MAKE) -C BCR_XL



# -------------
# Miscellaneous
# -------------

# remove auto-generated files for Citrus
.PHONY: clean_Citrus
clean_Citrus:
	find ../Citrus_files -type f -delete

# remove auto-generated files for CellCnn
.PHONY: clean_CellCnn
clean_CellCnn:
	find ../CellCnn_files -type f -delete

# remove auto-generated files for cydar
.PHONY: clean_cydar
clean_cydar:
	find ../cydar_files -type f -delete


# remove auto-generated 'Rplots.pdf' files
# (from plotting scripts and voom diagnostic plots)
.PHONY: clean_Rplots
clean_Rplots:
	find AML_sim/2_run_methods -name "Rplots.pdf" -delete
	find AML_sim/3_generate_plots -name "Rplots.pdf" -delete
	find BCR_XL_sim/2_run_methods -name "Rplots.pdf" -delete
	find BCR_XL_sim/3_generate_plots -name "Rplots.pdf" -delete
	find Anti_PD_1/run_methods -name "Rplots.pdf" -delete
	find Anti_PD_1/generate_plots -name "Rplots.pdf" -delete
	find BCR_XL/run_methods -name "Rplots.pdf" -delete
	find BCR_XL/generate_plots -name "Rplots.pdf" -delete


# Save lists of benchmark data files and paths
.PHONY: benchmark_data_paths
benchmark_data_paths:
	cd benchmark_data_paths && Rscript benchmark_data_paths.R



