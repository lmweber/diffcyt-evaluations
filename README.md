# diffcyt-evaluations

This repository contains scripts to reproduce all evaluations and figures in our paper introducing `diffcyt`, a new computational framework for differential discovery in high-dimensional flow and mass cytometry based on high-resolution clustering and empirical Bayes moderated tests.

The `diffcyt` R package is available from the [diffcyt](https://github.com/lmweber/diffcyt) repository, and will be submitted to [Bioconductor](http://bioconductor.org/).

The package includes documentation and examples of usage, including extended workflow examples available in the package vignette.


## Contents

Scripts are organized by data set, with the following sub-directories for each data set:

- 1_prepare_data: Scripts to prepare benchmark data sets
- 2_run_methods: Scripts to run each method, including performance comparisons against other methods
- 3_generate_plots: Scripts to generate figures shown in main text and supplementary material

