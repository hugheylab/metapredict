# metapredict

`metapredict` enables meta-analysis of gene expression using the elastic net (i.e., `glmnet`). Even if you plan to run your own analysis and not the elastic net, `metapredict` makes it straightforward to normalize all your microarray data and map probes from various platforms to Entrez Gene IDs. For technical details, please check out the [paper](http://nar.oxfordjournals.org/content/43/12/e79.full).

## Installation

Install the package in your local version of R:
```R
source('https://bioconductor.org/biocLite.R')
biocLite(c('org.Hs.eg.db', 'org.Mm.eg.db', 'org.Dr.eg.db'))
install.packages('devtools')
devtools::install_github('jakejh/metapredict', repos=BiocInstaller::biocinstallRepos())
```

Or use a pre-built [docker image](https://hub.docker.com/r/jakejh/hugheyverse), which has all dependencies already installed:
```
docker pull jakejh/hugheyverse
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
