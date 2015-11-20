# metapredict

`metapredict` enables meta-analysis of gene expression using the elastic net (i.e., `glmnet`). For technical details, please check out the [paper](http://nar.oxfordjournals.org/content/43/12/e79.full).

## Installation

```R
source('http://bioconductor.org/biocLite.R')
biocLite(c('affy', 'annotate', 'GEOquery', 'org.Hs.eg.db', 'org.Mm.eg.db', 'impute', 'sva'))

install.packages(c('cba', 'devtools', 'dplyr', 'foreach', 'glmnet', 'gridExtra', 'pheatmap', 'RColorBrewer', 'reshape2', 'matrixStats'))

devtools::install_github('hadley/ggplot2')

devtools::install_github('jakejh/metapredict')
```

## Getting started

First download and prepare the data for the example meta-analysis: `vignette('prepare_example', package='metapredict')`. Then run the example meta-analysis: `vignette('run_example', package='metapredict')`. For detailed instructions about running your own meta-analysis, please read `vignette('guidelines', package='metapredict')`.
