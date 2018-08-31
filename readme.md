# metapredict
`metapredict` enables meta-analysis of gene expression using the elastic net (i.e., `glmnet`). Even if you don't want or need to run the elastic net, `metapredict` makes it straightforward to normalize all your microarray data and to map probes from various platforms to Entrez Gene IDs. For technical details, please check out the [paper](https://doi.org/10.1093/nar/gkv229).

## Installation
First install drat.
```R
install.packages('drat')
```

Then add the following line to your `.Rprofile` file (located at "~/.Rprofile"), which gets run every time R starts. See [here](https://csgillespie.github.io/efficientR/3-3-r-startup.html#r-startup) for details.
```R
drat::addRepo('hugheylab')
```

Now you can install the package.
```R
setRepositories(ind = c(1:5, 9))
install.packages('metapredict', type = 'source')
```
You can update the package using `update.packages()`.

## Docker
You can also use a pre-built [docker image](https://hub.docker.com/r/hugheylab/hugheyverse), which has all dependencies installed.
```bash
docker pull hugheylab/hugheyverse
```

## Getting started
First follow the instructions to download and prepare the data for the example meta-analysis.
```R
vignette('prepare_example', package = 'metapredict'))
```

Then work through the example meta-analysis.
```R
vignette('run_example', package = 'metapredict'))
```

Finally, consult the detailed guide to running your own meta-analysis.
```R
vignette('guide', package = 'metapredict'))
```
