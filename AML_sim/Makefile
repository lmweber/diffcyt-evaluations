###################################################
# Makefile to run all analyses for AML-sim data set
###################################################


DIR_BENCHMARK      = ../../../benchmark_data/AML_sim/data

DIR_PREPARE_DATA   = 1_prepare_data
DIR_RUN_METHODS    = 2_run_methods
DIR_GENERATE_PLOTS = 3_generate_plots

DIR_RDATA          = ../../RData/AML_sim
DIR_PLOTS          = ../../plots/AML_sim


TIMESTAMP_BENCHMARK_MAIN             = $(DIR_BENCHMARK)/main/timestamp.txt
TIMESTAMP_BENCHMARK_RANDOM_SEEDS     = $(DIR_BENCHMARK)/random_seeds/timestamp.txt
TIMESTAMP_BENCHMARK_LESS_DISTINCT    = $(DIR_BENCHMARK)/less_distinct/timestamp.txt
TIMESTAMP_BENCHMARK_NULL_SIMULATIONS = $(DIR_BENCHMARK)/null_simulations/timestamp.txt

TIMESTAMP_PLOTS_DATA_DISTRIBUTIONS           = $(DIR_PLOTS)/data_distributions/timestamp.txt
TIMESTAMP_PLOTS_MAIN_PERFORMANCE             = $(DIR_PLOTS)/main_performance/timestamp.txt
TIMESTAMP_PLOTS_NULL_SIMULATIONS             = $(DIR_PLOTS)/null_simulations/timestamp.txt
TIMESTAMP_PLOTS_MAIN_TSNE                    = $(DIR_PLOTS)/main_tSNE/timestamp.txt
TIMESTAMP_PLOTS_MAIN_HEATMAPS                = $(DIR_PLOTS)/main_heatmaps/timestamp.txt
TIMESTAMP_PLOTS_MAIN_CLUSTERING_PERFORMANCE  = $(DIR_PLOTS)/main_clustering_performance/timestamp.txt
TIMESTAMP_PLOTS_SUPP_RANDOM_EFFECTS_LIMMA    = $(DIR_PLOTS)/supp_random_effects_limma/timestamp.txt
TIMESTAMP_PLOTS_SUPP_METACLUSTERING          = $(DIR_PLOTS)/supp_metaclustering/timestamp.txt
TIMESTAMP_PLOTS_SUPP_SEPARATE_CLUSTERING     = $(DIR_PLOTS)/supp_separate_clustering/timestamp.txt
TIMESTAMP_PLOTS_SUPP_RANDOM_SEEDS_CLUSTERING = $(DIR_PLOTS)/supp_random_seeds_clustering/timestamp.txt
TIMESTAMP_PLOTS_SUPP_RANDOM_SEEDS_DATA       = $(DIR_PLOTS)/supp_random_seeds_data/timestamp.txt
TIMESTAMP_PLOTS_SUPP_LESS_DISTINCT           = $(DIR_PLOTS)/supp_less_distinct/timestamp.txt
TIMESTAMP_PLOTS_SUPP_SAMPLE_SIZES            = $(DIR_PLOTS)/supp_sample_sizes/timestamp.txt
TIMESTAMP_PLOTS_SUPP_CLUSTERING_RESOLUTION   = $(DIR_PLOTS)/supp_clustering_resolution/timestamp.txt
TIMESTAMP_PLOTS_COMPARISONS_CITRUS           = $(DIR_PLOTS)/comparisons_Citrus/timestamp.txt
TIMESTAMP_PLOTS_COMPARISONS_CELLCNN          = $(DIR_PLOTS)/comparisons_CellCnn/timestamp.txt
TIMESTAMP_PLOTS_COMPARISONS_CYDAR            = $(DIR_PLOTS)/comparisons_cydar/timestamp.txt
TIMESTAMP_PLOTS_COMPARISONS_RUNTIMES         = $(DIR_PLOTS)/comparisons_runtimes/timestamp.txt

SCRIPTS_PLOTS_DATA_DISTRIBUTIONS           = $(wildcard $(DIR_GENERATE_PLOTS)/data_distributions/*.R)
SCRIPTS_PLOTS_MAIN_PERFORMANCE             = $(wildcard $(DIR_GENERATE_PLOTS)/main_performance/*.R)
SCRIPTS_PLOTS_NULL_SIMULATIONS             = $(wildcard $(DIR_GENERATE_PLOTS)/null_simulations/*.R)
SCRIPTS_PLOTS_MAIN_TSNE                    = $(wildcard $(DIR_GENERATE_PLOTS)/main_tSNE/*.R)
SCRIPTS_PLOTS_MAIN_HEATMAPS                = $(wildcard $(DIR_GENERATE_PLOTS)/main_heatmaps/*.R)
SCRIPTS_PLOTS_MAIN_CLUSTERING_PERFORMANCE  = $(wildcard $(DIR_GENERATE_PLOTS)/main_clustering_performance/*.R)
SCRIPTS_PLOTS_SUPP_RANDOM_EFFECTS_LIMMA    = $(wildcard $(DIR_GENERATE_PLOTS)/supp_random_effects_limma/*.R)
SCRIPTS_PLOTS_SUPP_METACLUSTERING          = $(wildcard $(DIR_GENERATE_PLOTS)/supp_metaclustering/*.R)
SCRIPTS_PLOTS_SUPP_SEPARATE_CLUSTERING     = $(wildcard $(DIR_GENERATE_PLOTS)/supp_separate_clustering/*.R)
SCRIPTS_PLOTS_SUPP_RANDOM_SEEDS_CLUSTERING = $(wildcard $(DIR_GENERATE_PLOTS)/supp_random_seeds_clustering/*.R)
SCRIPTS_PLOTS_SUPP_RANDOM_SEEDS_DATA       = $(wildcard $(DIR_GENERATE_PLOTS)/supp_random_seeds_data/*.R)
SCRIPTS_PLOTS_SUPP_LESS_DISTINCT           = $(wildcard $(DIR_GENERATE_PLOTS)/supp_less_distinct/*.R)
SCRIPTS_PLOTS_SUPP_SAMPLE_SIZES            = $(wildcard $(DIR_GENERATE_PLOTS)/supp_sample_sizes/*.R)
SCRIPTS_PLOTS_SUPP_CLUSTERING_RESOLUTION   = $(wildcard $(DIR_GENERATE_PLOTS)/supp_clustering_resolution/*.R)
SCRIPTS_PLOTS_COMPARISONS_CITRUS           = $(wildcard $(DIR_GENERATE_PLOTS)/comparisons_Citrus/*.R)
SCRIPTS_PLOTS_COMPARISONS_CELLCNN          = $(wildcard $(DIR_GENERATE_PLOTS)/comparisons_CellCnn/*.R)
SCRIPTS_PLOTS_COMPARISONS_CYDAR            = $(wildcard $(DIR_GENERATE_PLOTS)/comparisons_cydar/*.R)
SCRIPTS_PLOTS_COMPARISONS_RUNTIMES         = $(wildcard $(DIR_GENERATE_PLOTS)/comparisons_runtimes/*.R)

FILES_RDATA_MAIN                         = $(wildcard $(DIR_RDATA)/main/*.RData)
FILES_RDATA_NULL_SIMULATIONS             = $(wildcard $(DIR_RDATA)/null_simulations/*.RData)
FILES_RDATA_SUPP_RANDOM_EFFECTS_LIMMA    = $(wildcard $(DIR_RDATA)/supp_random_effects_limma/*.RData)
FILES_RDATA_SUPP_METACLUSTERING          = $(wildcard $(DIR_RDATA)/supp_metaclustering/*.RData)
FILES_RDATA_SUPP_SEPARATE_CLUSTERING     = $(wildcard $(DIR_RDATA)/supp_separate_clustering/*.RData)
FILES_RDATA_SUPP_RANDOM_SEEDS_CLUSTERING = $(wildcard $(DIR_RDATA)/supp_random_seeds_clustering/*.RData)
FILES_RDATA_SUPP_RANDOM_SEEDS_DATA       = $(wildcard $(DIR_RDATA)/supp_random_seeds_data/*.RData)
FILES_RDATA_SUPP_LESS_DISTINCT           = $(wildcard $(DIR_RDATA)/supp_less_distinct/*.RData)
FILES_RDATA_SUPP_SAMPLE_SIZES            = $(wildcard $(DIR_RDATA)/supp_sample_sizes/*.RData)
FILES_RDATA_SUPP_CLUSTERING_RESOLUTION   = $(wildcard $(DIR_RDATA)/supp_clustering_resolution/*.RData)
FILES_RDATA_COMPARISONS_CITRUS           = $(wildcard $(DIR_RDATA)/comparisons_Citrus/*.RData)
FILES_RDATA_COMPARISONS_CELLCNN          = $(wildcard $(DIR_RDATA)/comparisons_CellCnn/*.RData)
FILES_RDATA_COMPARISONS_CYDAR            = $(wildcard $(DIR_RDATA)/comparisons_cydar/*.RData)




# ------------
# Default rule
# ------------

.PHONY: all
all: timestamps_plots

.PHONY: timestamps_plots
timestamps_plots: $(TIMESTAMP_PLOTS_DATA_DISTRIBUTIONS) \
                  $(TIMESTAMP_PLOTS_MAIN_PERFORMANCE) \
                  $(TIMESTAMP_PLOTS_NULL_SIMULATIONS) \
                  $(TIMESTAMP_PLOTS_MAIN_TSNE) \
                  $(TIMESTAMP_PLOTS_MAIN_HEATMAPS) \
                  $(TIMESTAMP_PLOTS_MAIN_CLUSTERING_PERFORMANCE) \
                  $(TIMESTAMP_PLOTS_SUPP_RANDOM_EFFECTS_LIMMA) \
                  $(TIMESTAMP_PLOTS_SUPP_METACLUSTERING) \
                  $(TIMESTAMP_PLOTS_SUPP_SEPARATE_CLUSTERING) \
                  $(TIMESTAMP_PLOTS_SUPP_RANDOM_SEEDS_CLUSTERING) \
                  $(TIMESTAMP_PLOTS_SUPP_RANDOM_SEEDS_DATA) \
                  $(TIMESTAMP_PLOTS_SUPP_LESS_DISTINCT) \
                  $(TIMESTAMP_PLOTS_SUPP_SAMPLE_SIZES) \
                  $(TIMESTAMP_PLOTS_SUPP_CLUSTERING_RESOLUTION) \
                  $(TIMESTAMP_PLOTS_COMPARISONS_CITRUS) \
                  $(TIMESTAMP_PLOTS_COMPARISONS_CELLCNN) \
                  $(TIMESTAMP_PLOTS_COMPARISONS_CYDAR) \
                  $(TIMESTAMP_PLOTS_COMPARISONS_RUNTIMES)



# --------------
# Generate plots
# --------------

# note: separate rule for each set of results (otherwise dependencies get too complex)

# data distributions
$(TIMESTAMP_PLOTS_DATA_DISTRIBUTIONS): $(SCRIPTS_PLOTS_DATA_DISTRIBUTIONS) $(TIMESTAMP_BENCHMARK_MAIN) $(TIMESTAMP_BENCHMARK_LESS_DISTINCT)
	cd $(DIR_GENERATE_PLOTS)/data_distributions && bash run_all.sh

# main results: performance metrics
$(TIMESTAMP_PLOTS_MAIN_PERFORMANCE): $(SCRIPTS_PLOTS_MAIN_PERFORMANCE) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_performance && bash run_all.sh

# null simulations
$(TIMESTAMP_PLOTS_NULL_SIMULATIONS): $(SCRIPTS_PLOTS_NULL_SIMULATIONS) $(FILES_RDATA_NULL_SIMULATIONS)
	cd $(DIR_GENERATE_PLOTS)/null_simulations && bash run_all.sh

# main results: t-SNE plots
$(TIMESTAMP_PLOTS_MAIN_TSNE): $(SCRIPTS_PLOTS_MAIN_TSNE) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_tSNE && bash run_all.sh

# main results: heatmaps
$(TIMESTAMP_PLOTS_MAIN_HEATMAPS): $(SCRIPTS_PLOTS_MAIN_HEATMAPS) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_heatmaps && bash run_all.sh

# main results: clustering performance
$(TIMESTAMP_PLOTS_MAIN_CLUSTERING_PERFORMANCE): $(SCRIPTS_PLOTS_MAIN_CLUSTERING_PERFORMANCE) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_clustering_performance && bash run_all.sh


# supplementary: random effects for diffcyt-DA-limma
$(TIMESTAMP_PLOTS_SUPP_RANDOM_EFFECTS_LIMMA): $(SCRIPTS_PLOTS_SUPP_RANDOM_EFFECTS_LIMMA) $(FILES_RDATA_SUPP_RANDOM_EFFECTS_LIMMA)
	cd $(DIR_GENERATE_PLOTS)/supp_random_effects_limma && bash run_all.sh

# supplementary: meta-clustering
$(TIMESTAMP_PLOTS_SUPP_METACLUSTERING): $(SCRIPTS_PLOTS_SUPP_METACLUSTERING) $(FILES_RDATA_SUPP_METACLUSTERING)
	cd $(DIR_GENERATE_PLOTS)/supp_metaclustering && bash run_all.sh

# supplementary: separate clustering
$(TIMESTAMP_PLOTS_SUPP_SEPARATE_CLUSTERING): $(SCRIPTS_PLOTS_SUPP_SEPARATE_CLUSTERING) $(FILES_RDATA_SUPP_SEPARATE_CLUSTERING)
	cd $(DIR_GENERATE_PLOTS)/supp_separate_clustering && bash run_all.sh

# supplementary: varying random seeds for clustering
$(TIMESTAMP_PLOTS_SUPP_RANDOM_SEEDS_CLUSTERING): $(SCRIPTS_PLOTS_SUPP_RANDOM_SEEDS_CLUSTERING) $(FILES_RDATA_SUPP_RANDOM_SEEDS_CLUSTERING) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/supp_random_seeds_clustering && bash run_all.sh

# supplementary: varying random seeds for data generation
$(TIMESTAMP_PLOTS_SUPP_RANDOM_SEEDS_DATA): $(SCRIPTS_PLOTS_SUPP_RANDOM_SEEDS_DATA) $(FILES_RDATA_SUPP_RANDOM_SEEDS_DATA) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/supp_random_seeds_data && bash run_all.sh

# supplementary: 'less distinct' benchmark data
$(TIMESTAMP_PLOTS_SUPP_LESS_DISTINCT): $(SCRIPTS_PLOTS_SUPP_LESS_DISTINCT) $(FILES_RDATA_SUPP_LESS_DISTINCT) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/supp_less_distinct && bash run_all.sh

# supplementary: sample sizes
$(TIMESTAMP_PLOTS_SUPP_SAMPLE_SIZES): $(SCRIPTS_PLOTS_SUPP_SAMPLE_SIZES) $(FILES_RDATA_SUPP_SAMPLE_SIZES)
	cd $(DIR_GENERATE_PLOTS)/supp_sample_sizes && bash run_all.sh

# supplementary: clustering resolution
$(TIMESTAMP_PLOTS_SUPP_CLUSTERING_RESOLUTION): $(SCRIPTS_PLOTS_SUPP_CLUSTERING_RESOLUTION) $(FILES_RDATA_SUPP_CLUSTERING_RESOLUTION)
	cd $(DIR_GENERATE_PLOTS)/supp_clustering_resolution && bash run_all.sh


# comparisons: Citrus
$(TIMESTAMP_PLOTS_COMPARISONS_CITRUS): $(SCRIPTS_PLOTS_COMPARISONS_CITRUS) $(FILES_RDATA_COMPARISONS_CITRUS) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/comparisons_Citrus && bash run_all.sh

# comparisons: CellCnn
$(TIMESTAMP_PLOTS_COMPARISONS_CELLCNN): $(SCRIPTS_PLOTS_COMPARISONS_CELLCNN) $(FILES_RDATA_COMPARISONS_CELLCNN) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/comparisons_CellCnn && bash run_all.sh

# comparisons: cydar
$(TIMESTAMP_PLOTS_COMPARISONS_CYDAR): $(SCRIPTS_PLOTS_COMPARISONS_CYDAR) $(FILES_RDATA_COMPARISONS_CYDAR) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/comparisons_cydar && bash run_all.sh


# comparisons: runtimes
$(TIMESTAMP_PLOTS_COMPARISONS_RUNTIMES): $(SCRIPTS_PLOTS_COMPARISONS_RUNTIMES) $(FILES_RDATA_MAIN) $(FILES_RDATA_COMPARISONS_CITRUS) $(FILES_RDATA_COMPARISONS_CELLCNN) $(FILES_RDATA_COMPARISONS_CYDAR)
	cd $(DIR_GENERATE_PLOTS)/comparisons_runtimes && bash run_all.sh




# -----------
# Run methods
# -----------

# note: can't use as default rule if target contains '%'

# main results: performance metrics
$(DIR_RDATA)/main/outputs_%.RData: $(DIR_RUN_METHODS)/main/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/main && Rscript $(notdir $<)

# null simulations
$(DIR_RDATA)/null_simulations/outputs_%.RData: $(DIR_RUN_METHODS)/null_simulations/run_%.R $(TIMESTAMP_BENCHMARK_NULL_SIMULATIONS)
	cd $(DIR_RUN_METHODS)/null_simulations && Rscript $(notdir $<)


# supplementary: random effects for diffcyt-DS-limma
$(DIR_RDATA)/supp_random_effects_limma/outputs_%.RData: $(DIR_RUN_METHODS)/supp_random_effects_limma/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_random_effects_limma && Rscript $(notdir $<)

# supplementary: meta-clustering
$(DIR_RDATA)/supp_metaclustering/outputs_%.RData: $(DIR_RUN_METHODS)/supp_metaclustering/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_metaclustering && Rscript $(notdir $<)

# supplementary: separate clustering
$(DIR_RDATA)/supp_separate_clustering/outputs_%.RData: $(DIR_RUN_METHODS)/supp_separate_clustering/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_separate_clustering && Rscript $(notdir $<)

# supplementary: varying random seeds for clustering
$(DIR_RDATA)/supp_random_seeds_clustering/outputs_%.RData: $(DIR_RUN_METHODS)/supp_random_seeds_clustering/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_random_seeds_clustering && Rscript $(notdir $<)

# supplementary: varying random seeds for data generation
$(DIR_RDATA)/supp_random_seeds_data/outputs_%.RData: $(DIR_RUN_METHODS)/supp_random_seeds_data/run_%.R $(TIMESTAMP_BENCHMARK_RANDOM_SEEDS)
	cd $(DIR_RUN_METHODS)/supp_random_seeds_data && Rscript $(notdir $<)

# supplementary: 'less distinct' benchmark data
$(DIR_RDATA)/supp_less_distinct/outputs_%.RData: $(DIR_RUN_METHODS)/supp_less_distinct/run_%.R $(TIMESTAMP_BENCHMARK_LESS_DISTINCT)
	cd $(DIR_RUN_METHODS)/supp_less_distinct && Rscript $(notdir $<)

# supplementary: sample sizes
$(DIR_RDATA)/supp_sample_sizes/outputs_%.RData: $(DIR_RUN_METHODS)/supp_sample_sizes/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_sample_sizes && Rscript $(notdir $<)

# supplementary: clustering resolution
$(DIR_RDATA)/supp_clustering_resolution/outputs_%.RData: $(DIR_RUN_METHODS)/supp_clustering_resolution/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_clustering_resolution && Rscript $(notdir $<)


# comparisons: Citrus
$(DIR_RDATA)/comparisons_Citrus/outputs_%.RData: $(DIR_RUN_METHODS)/comparisons_Citrus/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/comparisons_Citrus && Rscript $(notdir $<)

# comparisons: CellCnn
$(DIR_RDATA)/comparisons_CellCnn/outputs_%.RData: $(DIR_RUN_METHODS)/comparisons_CellCnn/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/comparisons_CellCnn && Rscript $(notdir $<)

# comparisons: cydar
$(DIR_RDATA)/comparisons_cydar/outputs_%.RData: $(DIR_RUN_METHODS)/comparisons_cydar/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/comparisons_cydar && Rscript $(notdir $<)




# ---------------------------
# Prepare benchmark data sets
# ---------------------------

# generate benchmark data files: main
$(TIMESTAMP_BENCHMARK_MAIN): $(DIR_PREPARE_DATA)/main/prepare_data_AML_sim_main.R
	cd $(DIR_PREPARE_DATA)/main && Rscript $(notdir $<)

# generate benchmark data files: random seeds
$(TIMESTAMP_BENCHMARK_RANDOM_SEEDS): $(DIR_PREPARE_DATA)/random_seeds/prepare_data_AML_sim_random_seeds.R
	cd $(DIR_PREPARE_DATA)/random_seeds && Rscript $(notdir $<)

# generate benchmark data files: 'less distinct' data sets
$(TIMESTAMP_BENCHMARK_LESS_DISTINCT): $(DIR_PREPARE_DATA)/less_distinct/prepare_data_AML_sim_less_distinct.R
	cd $(DIR_PREPARE_DATA)/less_distinct && Rscript $(notdir $<)

# generate benchmark data files: null simulations
$(TIMESTAMP_BENCHMARK_NULL_SIMULATIONS): $(DIR_PREPARE_DATA)/null_simulations/prepare_data_AML_sim_null.R
	cd $(DIR_PREPARE_DATA)/null_simulations && Rscript $(notdir $<)



