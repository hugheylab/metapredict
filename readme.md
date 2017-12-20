# metapredict
`metapredict` enables meta-analysis of gene expression using the elastic net (i.e., `glmnet`). Even if you don't want or need to run the elastic net, `metapredict` makes it straightforward to normalize all your microarray data and map probes from various platforms to Entrez Gene IDs. For technical details, please check out the [paper](http://nar.oxfordjournals.org/content/43/12/e79.full).

## Install using drat
```R
install.packages('drat')
setRepositories(ind=1:5)
drat::addRepo('hugheylab')
install.packages('metapredict', type='source')
```
You can then update the package using `update.packages`.

## Install using devtools
```R
install.packages('devtools')
setRepositories(ind=1:5)
devtools::install_github('hugheylab/metapredict')
```
You can update the package using these same three lines.

## Install using docker
You can use a pre-built [docker image](https://hub.docker.com/r/hugheylab/hugheyverse), which has all dependencies already installed:
```
docker pull hugheylab/hugheyverse
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
