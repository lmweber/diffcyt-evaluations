######################################################
# Makefile to run all analyses for BCR-XL-sim data set
######################################################


DIR_BENCHMARK      = ../../../benchmark_data/BCR_XL_sim/data

DIR_PREPARE_DATA   = 1_prepare_data
DIR_RUN_METHODS    = 2_run_methods
DIR_GENERATE_PLOTS = 3_generate_plots

DIR_RDATA          = ../../RData/BCR_XL_sim
DIR_PLOTS          = ../../plots/BCR_XL_sim


TIMESTAMP_BENCHMARK_MAIN            = $(DIR_BENCHMARK)/main/timestamp.txt
TIMESTAMP_BENCHMARK_NULL_SIMULATION = $(DIR_BENCHMARK)/null_simulation/timestamp.txt
TIMESTAMP_BENCHMARK_RANDOM_SEEDS    = $(DIR_BENCHMARK)/random_seeds/timestamp.txt
TIMESTAMP_BENCHMARK_LESS_DISTINCT   = $(DIR_BENCHMARK)/less_distinct/timestamp.txt

TIMESTAMP_PLOTS_MAIN_PERFORMANCE             = $(DIR_PLOTS)/main_performance/timestamp.txt
TIMESTAMP_PLOTS_MAIN_RUNTIME                 = $(DIR_PLOTS)/main_runtime/timestamp.txt
TIMESTAMP_PLOTS_MAIN_TSNE                    = $(DIR_PLOTS)/main_tSNE/timestamp.txt
TIMESTAMP_PLOTS_MAIN_ONE_SENSE               = $(DIR_PLOTS)/main_One_SENSE/timestamp.txt
TIMESTAMP_PLOTS_MAIN_HEATMAPS                = $(DIR_PLOTS)/main_heatmaps/timestamp.txt
TIMESTAMP_PLOTS_NULL_SIMULATION              = $(DIR_PLOTS)/null_simulation/timestamp.txt
TIMESTAMP_PLOTS_SUPP_RANDOM_EFFECTS          = $(DIR_PLOTS)/supp_random_effects/timestamp.txt
TIMESTAMP_PLOTS_SUPP_METACLUSTERING          = $(DIR_PLOTS)/supp_metaclustering/timestamp.txt
TIMESTAMP_PLOTS_SUPP_RANDOM_SEEDS_CLUSTERING = $(DIR_PLOTS)/supp_random_seeds_clustering/timestamp.txt
TIMESTAMP_PLOTS_SUPP_RANDOM_SEEDS_DATA       = $(DIR_PLOTS)/supp_random_seeds_data/timestamp.txt
TIMESTAMP_PLOTS_SUPP_CLUSTERING_RESOLUTION   = $(DIR_PLOTS)/supp_clustering_resolution/timestamp.txt
TIMESTAMP_PLOTS_SUPP_SAMPLE_SIZES            = $(DIR_PLOTS)/supp_sample_sizes/timestamp.txt
TIMESTAMP_PLOTS_SUPP_LESS_DISTINCT           = $(DIR_PLOTS)/supp_less_distinct/timestamp.txt
TIMESTAMP_PLOTS_COMPARISONS_CYDAR            = $(DIR_PLOTS)/comparisons_cydar/timestamp.txt

SCRIPTS_PLOTS_MAIN_PERFORMANCE             = $(wildcard $(DIR_GENERATE_PLOTS)/main_performance/*.R)
SCRIPTS_PLOTS_MAIN_RUNTIME                 = $(wildcard $(DIR_GENERATE_PLOTS)/main_runtime/*.R)
SCRIPTS_PLOTS_MAIN_TSNE                    = $(wildcard $(DIR_GENERATE_PLOTS)/main_tSNE/*.R)
SCRIPTS_PLOTS_MAIN_ONE_SENSE               = $(wildcard $(DIR_GENERATE_PLOTS)/main_One_SENSE/*.R)
SCRIPTS_PLOTS_MAIN_HEATMAPS                = $(wildcard $(DIR_GENERATE_PLOTS)/main_heatmaps/*.R)
SCRIPTS_PLOTS_NULL_SIMULATION              = $(wildcard $(DIR_GENERATE_PLOTS)/null_simulation/*.R)
SCRIPTS_PLOTS_SUPP_RANDOM_EFFECTS          = $(wildcard $(DIR_GENERATE_PLOTS)/supp_random_effects/*.R)
SCRIPTS_PLOTS_SUPP_METACLUSTERING          = $(wildcard $(DIR_GENERATE_PLOTS)/supp_metaclustering/*.R)
SCRIPTS_PLOTS_SUPP_RANDOM_SEEDS_CLUSTERING = $(wildcard $(DIR_GENERATE_PLOTS)/supp_random_seeds_clustering/*.R)
SCRIPTS_PLOTS_SUPP_RANDOM_SEEDS_DATA       = $(wildcard $(DIR_GENERATE_PLOTS)/supp_random_seeds_data/*.R)
SCRIPTS_PLOTS_SUPP_CLUSTERING_RESOLUTION   = $(wildcard $(DIR_GENERATE_PLOTS)/supp_clustering_resolution/*.R)
SCRIPTS_PLOTS_SUPP_SAMPLE_SIZES            = $(wildcard $(DIR_GENERATE_PLOTS)/supp_sample_sizes/*.R)
SCRIPTS_PLOTS_SUPP_LESS_DISTINCT           = $(wildcard $(DIR_GENERATE_PLOTS)/supp_less_distinct/*.R)
SCRIPTS_PLOTS_COMPARISONS_CYDAR            = $(wildcard $(DIR_GENERATE_PLOTS)/comparisons_cydar/*.R)

FILES_RDATA_MAIN                         = $(wildcard $(DIR_RDATA)/main/*.RData)
FILES_RDATA_NULL_SIMULATION              = $(wildcard $(DIR_RDATA)/null_simulation/*.RData)
FILES_RDATA_SUPP_RANDOM_EFFECTS          = $(wildcard $(DIR_RDATA)/supp_random_effects/*.RData)
FILES_RDATA_SUPP_METACLUSTERING          = $(wildcard $(DIR_RDATA)/supp_metaclustering/*.RData)
FILES_RDATA_SUPP_RANDOM_SEEDS_CLUSTERING = $(wildcard $(DIR_RDATA)/supp_random_seeds_clustering/*.RData)
FILES_RDATA_SUPP_RANDOM_SEEDS_DATA       = $(wildcard $(DIR_RDATA)/supp_random_seeds_data/*.RData)
FILES_RDATA_SUPP_CLUSTERING_RESOLUTION   = $(wildcard $(DIR_RDATA)/supp_clustering_resolution/*.RData)
FILES_RDATA_SUPP_SAMPLE_SIZES            = $(wildcard $(DIR_RDATA)/supp_sample_sizes/*.RData)
FILES_RDATA_SUPP_LESS_DISTINCT           = $(wildcard $(DIR_RDATA)/supp_less_distinct/*.RData)
FILES_RDATA_COMPARISONS_CYDAR       = $(wildcard $(DIR_RDATA)/comparisons_cydar/*.RData)




# ------------
# Default rule
# ------------

.PHONY: all
all: timestamps_plots

.PHONY: timestamps_plots
timestamps_plots: $(TIMESTAMP_PLOTS_MAIN_PERFORMANCE) \
                  $(TIMESTAMP_PLOTS_MAIN_RUNTIME) \
                  $(TIMESTAMP_PLOTS_MAIN_TSNE) \
                  $(TIMESTAMP_PLOTS_MAIN_ONE_SENSE) \
                  $(TIMESTAMP_PLOTS_MAIN_HEATMAPS) \
                  $(TIMESTAMP_PLOTS_NULL_SIMULATION) \
                  $(TIMESTAMP_PLOTS_SUPP_RANDOM_EFFECTS) \
                  $(TIMESTAMP_PLOTS_SUPP_METACLUSTERING) \
                  $(TIMESTAMP_PLOTS_SUPP_RANDOM_SEEDS_CLUSTERING) \
                  $(TIMESTAMP_PLOTS_SUPP_RANDOM_SEEDS_DATA) \
                  $(TIMESTAMP_PLOTS_SUPP_CLUSTERING_RESOLUTION) \
                  $(TIMESTAMP_PLOTS_SUPP_SAMPLE_SIZES) \
                  $(TIMESTAMP_PLOTS_SUPP_LESS_DISTINCT) \
                  $(TIMESTAMP_PLOTS_COMPARISONS_CYDAR)




# --------------
# Generate plots
# --------------

# note: separate rule for each set of results (otherwise dependencies get too complex)

# main results: performance metrics
$(TIMESTAMP_PLOTS_MAIN_PERFORMANCE): $(SCRIPTS_PLOTS_MAIN_PERFORMANCE) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_performance && bash run_all.sh

# main results: runtime
$(TIMESTAMP_PLOTS_MAIN_RUNTIME): $(SCRIPTS_PLOTS_MAIN_RUNTIME) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_runtime && bash run_all.sh

# main results: tSNE plots
$(TIMESTAMP_PLOTS_MAIN_TSNE): $(SCRIPTS_PLOTS_MAIN_TSNE) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_tSNE && bash run_all.sh

# main results: One-SENSE plots
$(TIMESTAMP_PLOTS_MAIN_ONE_SENSE): $(SCRIPTS_PLOTS_MAIN_ONE_SENSE) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_One_SENSE && bash run_all.sh

# main results: heatmaps
$(TIMESTAMP_PLOTS_MAIN_HEATMAPS): $(SCRIPTS_PLOTS_MAIN_HEATMAPS) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/main_heatmaps && bash run_all.sh

# null simulation
$(TIMESTAMP_PLOTS_NULL_SIMULATION): $(SCRIPTS_PLOTS_NULL_SIMULATION) $(FILES_RDATA_NULL_SIMULATION)
	cd $(DIR_GENERATE_PLOTS)/null_simulation && bash run_all.sh

# supplementary: random effects
$(TIMESTAMP_PLOTS_SUPP_RANDOM_EFFECTS): $(SCRIPTS_PLOTS_SUPP_RANDOM_EFFECTS) $(FILES_RDATA_SUPP_RANDOM_EFFECTS)
	cd $(DIR_GENERATE_PLOTS)/supp_random_effects && bash run_all.sh

# supplementary: meta-clustering
$(TIMESTAMP_PLOTS_SUPP_METACLUSTERING): $(SCRIPTS_PLOTS_SUPP_METACLUSTERING) $(FILES_RDATA_SUPP_METACLUSTERING)
	cd $(DIR_GENERATE_PLOTS)/supp_metaclustering && bash run_all.sh

# supplementary: varying random seeds for clustering
$(TIMESTAMP_PLOTS_SUPP_RANDOM_SEEDS_CLUSTERING): $(SCRIPTS_PLOTS_SUPP_RANDOM_SEEDS_CLUSTERING) $(FILES_RDATA_SUPP_RANDOM_SEEDS_CLUSTERING)
	cd $(DIR_GENERATE_PLOTS)/supp_random_seeds_clustering && bash run_all.sh

# supplementary: varying random seeds for data generation
$(TIMESTAMP_PLOTS_SUPP_RANDOM_SEEDS_DATA): $(SCRIPTS_PLOTS_SUPP_RANDOM_SEEDS_DATA) $(FILES_RDATA_SUPP_RANDOM_SEEDS_DATA)
	cd $(DIR_GENERATE_PLOTS)/supp_random_seeds_data && bash run_all.sh

# supplementary: clustering resolution
$(TIMESTAMP_PLOTS_SUPP_CLUSTERING_RESOLUTION): $(SCRIPTS_PLOTS_SUPP_CLUSTERING_RESOLUTION) $(FILES_RDATA_SUPP_CLUSTERING_RESOLUTION)
	cd $(DIR_GENERATE_PLOTS)/supp_clustering_resolution && bash run_all.sh

# supplementary: sample sizes
$(TIMESTAMP_PLOTS_SUPP_SAMPLE_SIZES): $(SCRIPTS_PLOTS_SUPP_SAMPLE_SIZES) $(FILES_RDATA_SUPP_SAMPLE_SIZES)
	cd $(DIR_GENERATE_PLOTS)/supp_sample_sizes && bash run_all.sh

# supplementary: 'less distinct' benchmark data
$(TIMESTAMP_PLOTS_SUPP_LESS_DISTINCT): $(SCRIPTS_PLOTS_SUPP_LESS_DISTINCT) $(FILES_RDATA_SUPP_LESS_DISTINCT)
	cd $(DIR_GENERATE_PLOTS)/supp_less_distinct && bash run_all.sh

# comparisons: cydar
$(TIMESTAMP_PLOTS_COMPARISONS_CYDAR): $(SCRIPTS_PLOTS_COMPARISONS_CYDAR) $(FILES_RDATA_COMPARISONS_CYDAR) $(FILES_RDATA_MAIN)
	cd $(DIR_GENERATE_PLOTS)/comparisons_cydar && bash run_all.sh




# -----------
# Run methods
# -----------

# note: can't use as default rule if target contains '%'

# main results: performance metrics
$(DIR_RDATA)/main/outputs_%.RData: $(DIR_RUN_METHODS)/main/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/main && Rscript $(notdir $<)

# null simulation
$(DIR_RDATA)/null_simulation/outputs_%.RData: $(DIR_RUN_METHODS)/null_simulation/run_%.R $(TIMESTAMP_BENCHMARK_NULL_SIMULATION)
	cd $(DIR_RUN_METHODS)/null_simulation && Rscript $(notdir $<)

# supplementary: random effects
$(DIR_RDATA)/supp_random_effects/outputs_%.RData: $(DIR_RUN_METHODS)/supp_random_effects/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_random_effects && Rscript $(notdir $<)

# supplementary: meta-clustering
$(DIR_RDATA)/supp_metaclustering/outputs_%.RData: $(DIR_RUN_METHODS)/supp_metaclustering/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_metaclustering && Rscript $(notdir $<)

# supplementary: varying random seeds for clustering
$(DIR_RDATA)/supp_random_seeds_clustering/outputs_%.RData: $(DIR_RUN_METHODS)/supp_random_seeds_clustering/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_random_seeds_clustering && Rscript $(notdir $<)

# supplementary: varying random seeds for data generation
$(DIR_RDATA)/supp_random_seeds_data/outputs_%.RData: $(DIR_RUN_METHODS)/supp_random_seeds_data/run_%.R $(TIMESTAMP_BENCHMARK_RANDOM_SEEDS)
	cd $(DIR_RUN_METHODS)/supp_random_seeds_data && Rscript $(notdir $<)

# supplementary: clustering resolution
$(DIR_RDATA)/supp_clustering_resolution/outputs_%.RData: $(DIR_RUN_METHODS)/supp_clustering_resolution/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_clustering_resolution && Rscript $(notdir $<)

# supplementary: sample sizes
$(DIR_RDATA)/supp_sample_sizes/outputs_%.RData: $(DIR_RUN_METHODS)/supp_sample_sizes/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/supp_sample_sizes && Rscript $(notdir $<)

# supplementary: 'less distinct' benchmark data
$(DIR_RDATA)/supp_less_distinct/outputs_%.RData: $(DIR_RUN_METHODS)/supp_less_distinct/run_%.R $(TIMESTAMP_BENCHMARK_LESS_DISTINCT)
	cd $(DIR_RUN_METHODS)/supp_less_distinct && Rscript $(notdir $<)

# comparisons: cydar
$(DIR_RDATA)/comparisons_cydar/outputs_%.RData: $(DIR_RUN_METHODS)/comparisons_cydar/run_%.R $(TIMESTAMP_BENCHMARK_MAIN)
	cd $(DIR_RUN_METHODS)/comparisons_cydar && Rscript $(notdir $<)




# ---------------------------
# Prepare benchmark data sets
# ---------------------------

# generate benchmark data files: main
$(TIMESTAMP_BENCHMARK_MAIN): $(DIR_PREPARE_DATA)/main/prepare_data_BCR_XL_sim_main.R
	cd $(DIR_PREPARE_DATA)/main && Rscript $(notdir $<)

# generate benchmark data files: null simulation
$(TIMESTAMP_BENCHMARK_NULL_SIMULATION): $(DIR_PREPARE_DATA)/null_simulation/prepare_data_BCR_XL_sim_null.R
	cd $(DIR_PREPARE_DATA)/null_simulation && Rscript $(notdir $<)

# generate benchmark data files: random seeds
$(TIMESTAMP_BENCHMARK_RANDOM_SEEDS): $(DIR_PREPARE_DATA)/random_seeds/prepare_data_BCR_XL_sim_random_seeds.R
	cd $(DIR_PREPARE_DATA)/random_seeds && Rscript $(notdir $<)

# generate benchmark data files: 'less distinct' data sets
$(TIMESTAMP_BENCHMARK_LESS_DISTINCT): $(DIR_PREPARE_DATA)/less_distinct/prepare_data_BCR_XL_sim_less_distinct.R
	cd $(DIR_PREPARE_DATA)/less_distinct && Rscript $(notdir $<)



