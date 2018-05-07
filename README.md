# diffcyt-evaluations

This repository contains scripts to reproduce all performance evaluations, comparisons, and figures in our paper introducing the `diffcyt` framework.

The `diffcyt` R package implements statistical methods for differential discovery analyses in high-dimensional cytometry data (including flow cytometry, mass cytometry or CyTOF, and oligonucleotide-tagged cytometry), based on (i) high-resolution clustering and (ii) moderated tests adapted from transcriptomics.


## Contents

Scripts are organized by data set (`AML-sim`, `BCR-XL-sim`, `Anti-PD-1`, and `BCR-XL`).

Within each data set, scripts are organized into sub-directories to:

- prepare data (`AML-sim` and `BCR-XL-sim` only)
- run methods
- generate plots

Code comments are included to explain the purpose of each script.


## Preprint

A preprint of the paper will be made available from bioRxiv.


## Data files

Data files will be made available from FlowRepository.


## Installation of `diffcyt` package

The `diffcyt` package is available from [Bioconductor](http://bioconductor.org/packages/diffcyt). It can be installed using the Bioconductor installer (`biocLite`):

```{r}
# Download the Bioconductor installer
source("https://bioconductor.org/biocLite.R")

# Install 'diffcyt' package
biocLite("diffcyt")
```

## Tutorial and examples

For a tutorial and examples of usage, see the Bioconductor package vignette (accessible under 'Documentation' on the [Bioconductor page](http://bioconductor.org/packages/diffcyt)).


