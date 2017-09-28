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

# remove auto-generated diagnostic plots (voom) for diffcyt-DA-limma
.PHONY: clean_voom
clean_voom:
	find AML_sim/2_run_methods -name "Rplots.pdf" -delete


