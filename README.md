# diffcyt-evaluations

This repository contains scripts to reproduce all performance evaluations, comparisons, and figures in our paper introducing the `diffcyt` framework.

The `diffcyt` R package implements statistical methods for differential discovery analyses in high-dimensional cytometry data (including flow cytometry, mass cytometry or CyTOF, and oligonucleotide-tagged cytometry), based on (i) high-resolution clustering and (ii) empirical Bayes moderated tests adapted from transcriptomics.


## Contents

Scripts are organized by dataset: 'AML-sim', 'BCR-XL-sim', 'Anti-PD-1', and 'BCR-XL'.

Within each dataset, scripts are organized into sub-directories to:

- prepare data (semi-simulated datasets 'AML-sim' and 'BCR-XL-sim' only)
- run methods
- generate plots

Code comments are included to explain the purpose of each script.


## Preprint

A preprint of our paper introducing the `diffcyt` framework is available from bioRxiv:

- Weber L. M. et al. (2018), *diffcyt: Differential discovery in high-dimensional cytometry via high-resolution clustering*, bioRxiv preprint. [Available here.](https://www.biorxiv.org/content/early/2018/06/18/349738)


## Data files

Data files for the benchmarking datasets are available from FlowRepository under accession number [FR-FCM-ZYL8](http://flowrepository.org/id/FR-FCM-ZYL8).


## Installation of `diffcyt` package

The `diffcyt` package is freely available from [Bioconductor](http://bioconductor.org/packages/diffcyt). The stable release version can be installed using the Bioconductor installer as follows. Note that installation requires R version 3.5.0 or later.

```{r}
# Install Bioconductor installer from CRAN
install.packages("BiocManager")

# Install 'diffcyt' package from Bioconductor
BiocManager::install("diffcyt")
```

To run the examples in the package vignette, the `HDCytoData` and `CATALYST` packages from Bioconductor are also required.

```{r}
BiocManager::install("HDCytoData")
BiocManager::install("CATALYST")
```

For details on the development version of the `diffcyt` package, see the [GitHub page](https://github.com/lmweber/diffcyt).


## Tutorial and examples

For a tutorial and examples of usage, see the Bioconductor [package vignette](http://bioconductor.org/packages/release/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html) (link also available via the main Bioconductor page for the [diffcyt package](http://bioconductor.org/packages/diffcyt)).

<p> <img src="diffcyt.png" width="130"/> </p>


