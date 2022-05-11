# metapredict

[![check-deploy](https://github.com/hugheylab/metapredict/workflows/check-deploy/badge.svg)](https://github.com/hugheylab/metapredict/actions)
[![codecov](https://codecov.io/gh/hugheylab/metapredict/branch/master/graph/badge.svg?token=nRgjANwZ6s)](https://codecov.io/gh/hugheylab/metapredict)
[![Netlify Status](https://api.netlify.com/api/v1/badges/6ff1bdb7-ea07-4c5f-9a52-f27f9eae8554/deploy-status)](https://app.netlify.com/sites/infallible-spence-3f03b0/deploys)
[![drat version](https://raw.githubusercontent.com/hugheylab/drat/gh-pages/badges/metapredict_drat_badge.svg)](https://github.com/hugheylab/drat/tree/gh-pages/src/contrib)

`metapredict` enables meta-analysis of gene expression using the elastic net (i.e., `glmnet`). For technical details, please check out the [paper](https://doi.org/10.1093/nar/gkv229).

**If you just want to process microarray data and don't need to merge the expression data to perform a meta-analysis, then check out [seeker](https://seeker.hugheylab.org), which makes fetching, processing, and normalizing transcriptome data even simpler.**

## Installation

1. Install [`BiocManager`](https://cran.r-project.org/package=BiocManager).

    ```r
    if (!requireNamespace('BiocManager', quietly = TRUE))
      install.packages('BiocManager')
    ```

1. If you use RStudio, go to Tools → Global Options... → Packages → Add... (under Secondary repositories), then enter:

    - Name: hugheylab
    - Url: https://hugheylab.github.io/drat/

    You only have to do this once. Then you can install or update the package by entering:

    ```r
    BiocManager::install('metapredict')
    ```

    Alternatively, you can install or update the package by entering:

    ```r
    BiocManager::install('metapredict', site_repository = 'https://hugheylab.github.io/drat/')
    ```

## Usage

Go through the vignettes to

- [prepare the data for the example meta-analysis](https://metapredict.hugheylab.org/articles/prepare_example.html),
- [run the example meta-analysis](https://metapredict.hugheylab.org/articles/run_example.html), and
- [see how to run your own meta-analysis](https://metapredict.hugheylab.org/articles/guide.html).

For help on using specific functions, consult the [reference documentation](https://metapredict.hugheylab.org/reference/index.html).
