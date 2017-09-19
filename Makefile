###########################################################
# Top-level Makefile to run all Makefiles in subdirectories
###########################################################


all:
#	$(MAKE) -C BCR_XL_sim
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


