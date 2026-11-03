---
title: Setup
---

Please follow the steps below and install the required software **before** the scheduled workshop.

## RStudio Setup

We use RStudio for coding in R.

[Click here and follow the instructions](https://posit.co/download/rstudio-desktop/) to install RStudio Desktop in your system.

::::::::::::::::: discussion

### R packages

Most workshops using R will require the installation of specific packages. Make sure to check in advance with the workshop organisers what packages need to be installed.

You can install packages from CRAN using:

```r
install.packages("package_name")
```

If your package is in a different R repository, such as Bioconductor or GitHub, you may need the [BiocManager](https://www.bioconductor.org/install/) or [devtools](https://devtools.r-lib.org/) packages to install them. To install BiocManager:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

For devtools, you can simply do:

```r
install.packages("devtools")
```

You can then install packages directly from GitHub with:

```r
devtools::install_github("username/reponame")
```

::::::::::::::::::::::::::::

### Install required R packages

Run the following commands in the **Console** pane of RStudio to install the packages for this workshop. Required dependencies are installed automatically.

```r
install.packages(c("dplyr", "circlize"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "minfi", "limma", "GEOquery",
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
  "IlluminaHumanMethylationEPICv2manifest",
  "AnnotationHub", "ComplexHeatmap", "DMRcate", "missMethyl"
), update = FALSE, ask = FALSE)
```

For EPIC v2 arrays, the matching packages are [`IlluminaHumanMethylationEPICv2manifest`](https://bioconductor.org/packages/IlluminaHumanMethylationEPICv2manifest/) and [`IlluminaHumanMethylationEPICv2anno.20a1.hg38`](https://bioconductor.org/packages/IlluminaHumanMethylationEPICv2anno.20a1.hg38/). 
Manifest and annotation package required will change depending on array type used. 
