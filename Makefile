#############################################################
# Top-level Makefile to run other Makefiles in subdirectories
#############################################################


all:
	$(MAKE) -C AML_sim



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

# remove auto-generated 'Rplots.pdf' files
# (from heatmaps plotting script, and voom diagnostic plots from diffcyt-DA-limma)
.PHONY: clean_Rplots
clean_Rplots:
	find AML_sim/2_run_methods -name "Rplots.pdf" -delete
	find AML_sim/3_generate_plots -name "Rplots.pdf" -delete
	find BCR_XL_sim/2_run_methods -name "Rplots.pdf" -delete
	find BCR_XL_sim/3_generate_plots -name "Rplots.pdf" -delete



