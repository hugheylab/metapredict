# metapredict

`metapredict` enables meta-analysis of gene expression using the elastic net (i.e., `glmnet`). For technical details, please check out the [paper](http://nar.oxfordjournals.org/content/43/12/e79.full).

## Installation

```R
source('https://bioconductor.org/biocLite.R')
install.packages('devtools')
devtools::install_github('hadley/ggplot2')
devtools::install_github('jakejh/metapredict', repos=BiocInstaller::biocinstallRepos())
```

## Getting started

First follow the instructions to download and prepare the data for the example meta-analysis.
```R
file.show(system.file('extdata', 'prepare_example.html', package='metapredict'))
```

Then work through the example meta-analysis.
```R
file.show(system.file('extdata', 'run_example.html', package='metapredict'))
```

Finally, consult the detailed guide to running your own meta-analysis.
```R
file.show(system.file('extdata', 'guide.html', package='metapredict'))
```
