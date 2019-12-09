# `mirage`

An R package for MIxture model based Rare variant Analysis on GEnes. 

## Description

`mirage` is a new Bayesian statistical method for rare variant (RV) association testing that better accounts for heterogeneity of variant effects within a gene using 
external annotation information. It models variants in a gene as a mixture of risk and non-risk variants, with a prior probability of being a risk variant determined by functional annotations of the variant such as
conservation score and impact on protein structure. Since in general external annotations alone have limited accuracy in predicting functional effects, a simple filter based on such annotations (as commonly performed in 
many RV association analysis) may result in both false positive and negatives. Instead, by incorporating such information as prior and using a hierarchical model to pool information across genes, `mirage` is able to better
characterize the inclusion probability of rare variants for different functional categories, thus improving the power to detect an association.

## Quick Start

1. Follow the setup instructions below.

2. See the [Quick Start Example](https://xinhe-lab.github.io/mirage/articles/mwe.html) for a toy example to run `mirage`.

## Setup

To install the package,

```R
library("devtools")
install_github('xinhe-lab/mirage')
```

## Developer notes

To build documentation for the package,

```R
setwd("/path/to/package/root")
devtools::document()
```
Please **do not** manually edit any `.Rd` files under `man` folder!

To add tutorials, you write them as `.Rmd` files and put them under
`vignettes` folder. Then edit [this list](https://github.com/xinhe-lab/mirage/blob/fc6e9f664740996def10e6b35f2b754d91e4c329/_pkgdown.yml#L21)
simply adding to it the names of your `.Rmd` file (without `.Rmd` extension).

To build the documentation website (make sure you are
connected to Internet while running these commands):

```R
setwd("/path/to/package/root")
devtools::document()
pkgdown::build_site()
```

To install locally from source code via `devtools`, 

```R
setwd("/path/to/package/root")
devtools::document()
devtools::install()
```

To auto-format the code,
```R
setwd("/path/to/package/root")
formatR::tidy_dir('./R', wrap = 120)
```
