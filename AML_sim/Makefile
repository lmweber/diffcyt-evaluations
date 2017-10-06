###################################################
# Makefile to run all analyses for AML-sim data set
###################################################


DIR_BENCHMARK      = ../../../benchmark_data/AML_sim/data

DIR_PREPARE_DATA   = 1_prepare_data
DIR_RUN_METHODS    = 2_run_methods
DIR_GENERATE_PLOTS = 3_generate_plots

DIR_RDATA          = ../../RData/AML_sim
DIR_PLOTS          = ../../plots/AML_sim


TIMESTAMP_BENCHMARK_MAIN          = $(DIR_BENCHMARK)/main/timestamp.txt
TIMESTAMP_BENCHMARK_RANDOM        = $(DIR_BENCHMARK)/random_seeds/timestamp.txt
TIMESTAMP_BENCHMARK_LESS_DISTINCT = $(DIR_BENCHMARK)/less_distinct/timestamp.txt

TIMESTAMP_PLOTS_MAIN_ROC_TPRFDR            = $(DIR_PLOTS)/main_ROC_TPRFDR/timestamp.txt
TIMESTAMP_PLOTS_MAIN_RUNTIME               = $(DIR_PLOTS)/main_runtime/timestamp.txt
TIMESTAMP_PLOTS_MAIN_CLUSTERING_QUALITY    = $(DIR_PLOTS)/main_clustering_quality/timestamp.txt
TIMESTAMP_PLOTS_MAIN_MST                   = $(DIR_PLOTS)/main_MST/timestamp.txt
TIMESTAMP_PLOTS_MAIN_TSNE                  = $(DIR_PLOTS)/main_tSNE/timestamp.txt
TIMESTAMP_PLOTS_MAIN_HEATMAPS              = $(DIR_PLOTS)/main_heatmaps/timestamp.txt
TIMESTAMP_PLOTS_SUPP_CLUSTERING_RESOLUTION = $(DIR_PLOTS)/supp_clustering_resolution/timestamp.txt
TIMESTAMP_PLOTS_SUPP_FIXED_EFFECTS         = $(DIR_PLOTS)/supp_fixed_effects/timestamp.txt
TIMESTAMP_PLOTS_SUPP_ALL_MARKERS           = $(DIR_PLOTS)/supp_all_markers/timestamp.txt
# note no separate plots for 'supp_lineage_markers' (included in 'supp_all_markers')
TIMESTAMP_PLOTS_SUPP_SEPARATE_CLUSTERING   = $(DIR_PLOTS)/supp_separate_clustering/timestamp.txt
TIMESTAMP_PLOTS_SUPP_RANDOM_STARTS         = $(DIR_PLOTS)/supp_random_starts/timestamp.txt
TIMESTAMP_PLOTS_SUPP_BENCHMARK_RANDOM      = $(DIR_PLOTS)/supp_benchmark_random/timestamp.txt

SCRIPTS_PLOTS_MAIN_ROC_TPRFDR            = $(wildcard $(DIR_GENERATE_PLOTS)/main_ROC_TPRFDR/*.R)
SCRIPTS_PLOTS_MAIN_RUNTIME               = $(wildcard $(DIR_GENERATE_PLOTS)/main_runtime/*.R)
SCRIPTS_PLOTS_MAIN_CLUSTERING_QUALITY    = $(wildcard $(DIR_GENERATE_PLOTS)/main_clustering_quality/*.R)
SCRIPTS_PLOTS_MAIN_MST                   = $(wildcard $(DIR_GENERATE_PLOTS)/main_MST/*.R)
SCRIPTS_PLOTS_MAIN_TSNE                  = $(wildcard $(DIR_GENERATE_PLOTS)/main_tSNE/*.R)
SCRIPTS_PLOTS_MAIN_HEATMAPS              = $(wildcard $(DIR_GENERATE_PLOTS)/main_heatmaps/*.R)
SCRIPTS_PLOTS_SUPP_CLUSTERING_RESOLUTION = $(wildcard $(DIR_GENERATE_PLOTS)/supp_clustering_resolution/*.R)
SCRIPTS_PLOTS_SUPP_FIXED_EFFECTS         = $(wildcard $(DIR_GENERATE_PLOTS)/supp_fixed_effects/*.R)
SCRIPTS_PLOTS_SUPP_ALL_MARKERS           = $(wildcard $(DIR_GENERATE_PLOTS)/supp_all_markers/*.R)
SCRIPTS_PLOTS_SUPP_LINEAGE_MARKERS       = $(wildcard $(DIR_GENERATE_PLOTS)/supp_lineage_markers/*.R)
SCRIPTS_PLOTS_SUPP_SEPARATE_CLUSTERING   = $(wildcard $(DIR_GENERATE_PLOTS)/supp_separate_clustering/*.R)
SCRIPTS_PLOTS_SUPP_RANDOM_STARTS         = $(wildcard $(DIR_GENERATE_PLOTS)/supp_random_starts/*.R)
SCRIPTS_PLOTS_SUPP_BENCHMARK_RANDOM      = $(wildcard $(DIR_GENERATE_PLOTS)/supp_benchmark_random/*.R)

FILES_RDATA_MAIN                       = $(wildcard $(DIR_RDATA)/main/*.RData)
FILES_RDATA_SUPP_CLUSTERING_RESOLUTION = $(wildcard $(DIR_RDATA)/supp_clustering_resolution/*.RData)
FILES_RDATA_SUPP_FIXED_EFFECTS         = $(wildcard $(DIR_RDATA)/supp_fixed_effects/*.RData)
FILES_RDATA_SUPP_ALL_MARKERS           = $(wildcard $(DIR_RDATA)/supp_all_markers/*.RData)
FILES_RDATA_SUPP_LINEAGE_MARKERS       = $(wildcard $(DIR_RDATA)/supp_lineage_markers/*.RData)
FILES_RDATA_SUPP_SEPARATE_CLUSTERING   = $(wildcard $(DIR_RDATA)/supp_separate_clustering/*.RData)
FILES_RDATA_SUPP_RANDOM_STARTS         = $(wildcard $(DIR_RDATA)/supp_random_starts/*.RData)
# note subdirectories (e.g. 'supp_benchmark_random/seed1')
FILES_RDATA_SUPP_BENCHMARK_RANDOM      = $(wildcard $(DIR_RDATA)/supp_benchmark_random/*/*.RData)




# ------------
# Default rule
# ------------

.PHONY: all
all: timestamps_plots

.PHONY: timestamps_plots
timestamps_plots: $(TIMESTAMP_PLOTS_MAIN_ROC_TPRFDR) \
                  $(TIMESTAMP_PLOTS_MAIN_RUNTIME) \
                  $(TIMESTAMP_PLOTS_MAIN_CLUSTERING_QUALITY) \
                  $(TIMESTAMP_PLOTS_MAIN_MST) \
                  $(TIMESTAMP_PLOTS_MAIN_TSNE) \
                  $(TIMESTAMP_PLOTS_MAIN_HEATMAPS) \
                  $(TIMESTAMP_PLOTS_SUPP_CLUSTERING_RESOLUTION) \
                  $(TIMESTAMP_PLOTS_SUPP_FIXED_EFFECTS) \
                  $(TIMESTAMP_PLOTS_SUPP_ALL_MARKERS) \
                  $(TIMESTAMP_PLOTS_SUPP_SEPARATE_CLUSTERING) \
                  $(TIMESTAMP_PLOTS_SUPP_RANDOM_STARTS) \
                  $(TIMESTAMP_PLOTS_SUPP_BENCHMARK_RANDOM)




# --------------
# Generate plots
# --------------

# note: separate rule for each set of results (otherwise dependencies get too complex)

# main results: ROC and TPR-FDR curves
$(TIMESTAMP_PLOTS_MAIN_ROC_TPRFDR): $(SCRIPTS_PLOTS_MAIN_ROC_TPRFDR) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_ROC_TPRFDR && bash run_all.sh

# main results: runtime
$(TIMESTAMP_PLOTS_MAIN_RUNTIME): $(SCRIPTS_PLOTS_MAIN_RUNTIME) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_runtime && bash run_all.sh

# main results: clustering quality
$(TIMESTAMP_PLOTS_MAIN_CLUSTERING_QUALITY): $(SCRIPTS_PLOTS_MAIN_CLUSTERING_QUALITY) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_clustering_quality && bash run_all.sh

# main results: MST
$(TIMESTAMP_PLOTS_MAIN_MST): $(SCRIPTS_PLOTS_MAIN_MST) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_MST && bash run_all.sh

# main results: tSNE
$(TIMESTAMP_PLOTS_MAIN_TSNE): $(SCRIPTS_PLOTS_MAIN_TSNE) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_tSNE && bash run_all.sh

# main results: heatmaps
$(TIMESTAMP_PLOTS_MAIN_HEATMAPS): $(SCRIPTS_PLOTS_MAIN_HEATMAPS) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_heatmaps && bash run_all.sh

# supplementary: clustering resolution
$(TIMESTAMP_PLOTS_SUPP_CLUSTERING_RESOLUTION): $(SCRIPTS_PLOTS_SUPP_CLUSTERING_RESOLUTION) $(FILES_RDATA_SUPP_CLUSTERING_RESOLUTION)
	cd $(DIR_GENERATE_PLOTS)/supp_clustering_resolution && bash run_all.sh

# supplementary: diffcyt-DA-limma fixed effects
$(TIMESTAMP_PLOTS_SUPP_FIXED_EFFECTS): $(SCRIPTS_PLOTS_SUPP_FIXED_EFFECTS) $(FILES_RDATA_SUPP_FIXED_EFFECTS)
	cd $(DIR_GENERATE_PLOTS)/supp_fixed_effects && bash run_all.sh

# supplementary: all markers
# (note: also depends on 'lineage markers' since CellCnn uses 'all markers' for main results)
$(TIMESTAMP_PLOTS_SUPP_ALL_MARKERS): $(SCRIPTS_PLOTS_SUPP_ALL_MARKERS) $(FILES_RDATA_SUPP_ALL_MARKERS) $(SCRIPTS_PLOTS_SUPP_LINEAGE_MARKERS) $(FILES_RDATA_SUPP_LINEAGE_MARKERS)
	cd $(DIR_GENERATE_PLOTS)/supp_all_markers && bash run_all.sh

# supplementary: separate_clustering
$(TIMESTAMP_PLOTS_SUPP_SEPARATE_CLUSTERING): $(SCRIPTS_PLOTS_SUPP_SEPARATE_CLUSTERING) $(FILES_RDATA_SUPP_SEPARATE_CLUSTERING)
	cd $(DIR_GENERATE_PLOTS)/supp_separate_clustering && bash run_all.sh

# supplementary: random starts for clustering
$(TIMESTAMP_PLOTS_SUPP_RANDOM_STARTS): $(SCRIPTS_PLOTS_SUPP_RANDOM_STARTS) $(FILES_RDATA_SUPP_RANDOM_STARTS)
	cd $(DIR_GENERATE_PLOTS)/supp_random_starts && bash run_all.sh

# supplementary: randomized benchmark data sets
$(TIMESTAMP_PLOTS_SUPP_BENCHMARK_RANDOM): $(SCRIPTS_PLOTS_SUPP_BENCHMARK_RANDOM) $(FILES_RDATA_SUPP_BENCHMARK_RANDOM)
	cd $(DIR_GENERATE_PLOTS)/supp_benchmark_random && bash run_all.sh




# -----------
# Run methods
# -----------

# note: can't use as default rule if target contains '%'

# main results
$(DIR_RDATA)/main/outputs_%.RData: $(DIR_RUN_METHODS)/main/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/main && Rscript $(notdir $<)

# supplementary: clustering resolution
$(DIR_RDATA)/supp_clustering_resolution/outputs_%.RData: $(DIR_RUN_METHODS)/supp_clustering_resolution/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_clustering_resolution && Rscript $(notdir $<)

# supplementary: diffcyt-DA-limma fixed effects
$(DIR_RDATA)/supp_fixed_effects/outputs_%.RData: $(DIR_RUN_METHODS)/supp_fixed_effects/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_fixed_effects && Rscript $(notdir $<)

# supplementary: all markers
$(DIR_RDATA)/supp_all_markers/outputs_%.RData: $(DIR_RUN_METHODS)/supp_all_markers/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_all_markers && Rscript $(notdir $<)

# supplementary: lineage markers
$(DIR_RDATA)/supp_lineage_markers/outputs_%.RData: $(DIR_RUN_METHODS)/supp_lineage_markers/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_lineage_markers && Rscript $(notdir $<)

# supplementary: separate clustering
$(DIR_RDATA)/supp_separate_clustering/outputs_%.RData: $(DIR_RUN_METHODS)/supp_separate_clustering/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_separate_clustering && Rscript $(notdir $<)

# supplementary: random starts for clustering
$(DIR_RDATA)/supp_random_starts/outputs_%.RData: $(DIR_RUN_METHODS)/supp_random_starts/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_random_starts && Rscript $(notdir $<)

# supplementary: randomized benchmark data sets
$(DIR_RDATA)/supp_benchmark_random/seed1/outputs_%.RData: $(DIR_RUN_METHODS)/supp_benchmark_random/seed1/run_%.R $(TIMESTAMP_BENCHMARK_RANDOM)
	cd $(DIR_RUN_METHODS)/supp_benchmark_random/seed1 && Rscript $(notdir $<)

$(DIR_RDATA)/supp_benchmark_random/seed2/outputs_%.RData: $(DIR_RUN_METHODS)/supp_benchmark_random/seed2/run_%.R $(TIMESTAMP_BENCHMARK_RANDOM)
	cd $(DIR_RUN_METHODS)/supp_benchmark_random/seed2 && Rscript $(notdir $<)

$(DIR_RDATA)/supp_benchmark_random/seed3/outputs_%.RData: $(DIR_RUN_METHODS)/supp_benchmark_random/seed3/run_%.R $(TIMESTAMP_BENCHMARK_RANDOM)
	cd $(DIR_RUN_METHODS)/supp_benchmark_random/seed3 && Rscript $(notdir $<)

# supplementary: less distinct populations
$(DIR_RDATA)/supp_less_distinct/50pc/outputs_%.RData: $(DIR_RUN_METHODS)/supp_less_distinct/50pc/run_%.R $(TIMESTAMP_BENCHMARK_LESS_DISTINCT)
	cd $(DIR_RUN_METHODS)/supp_less_distinct/50pc && Rscript $(notdir $<)

$(DIR_RDATA)/supp_less_distinct/75pc/outputs_%.RData: $(DIR_RUN_METHODS)/supp_less_distinct/75pc/run_%.R $(TIMESTAMP_BENCHMARK_LESS_DISTINCT)
	cd $(DIR_RUN_METHODS)/supp_less_distinct/75pc && Rscript $(notdir $<)




# ---------------------------
# Prepare benchmark data sets
# ---------------------------

# generate benchmark data files: main
$(TIMESTAMP_BENCHMARK_MAIN): $(DIR_PREPARE_DATA)/main/prepare_data_AML_sim_main.R
	cd $(DIR_PREPARE_DATA)/main && Rscript $(notdir $<)

# generate benchmark data files: random seeds
$(TIMESTAMP_BENCHMARK_RANDOM): $(DIR_PREPARE_DATA)/random_seeds/prepare_data_AML_sim_random_seeds.R
	cd $(DIR_PREPARE_DATA)/random_seeds && Rscript $(notdir $<)

# generate benchmark data files: less distinct populations
$(TIMESTAMP_BENCHMARK_LESS_DISTINCT): $(DIR_PREPARE_DATA)/less_distinct/prepare_data_AML_sim_less_distinct.R
	cd $(DIR_PREPARE_DATA)/less distinct && Rscript $(notdir $<)



