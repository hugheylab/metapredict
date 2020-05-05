# metapredict

`metapredict` enables meta-analysis of gene expression using the elastic net (i.e., `glmnet`). Even if you don't want or need to run the elastic net, `metapredict` makes it straightforward to normalize all your microarray data and to map probes from various platforms to Entrez Gene IDs. For technical details, please check out the [paper](https://doi.org/10.1093/nar/gkv229).

## Installation

First add the hugheylab repository to your repos. There are multiple ways to do this.

If you use RStudio, go to Tools -> Global Options... -> Packages -> Add... (under Secondary repositories), then enter the following values.

- Name: hugheylab
- Url: https://hugheylab.github.io/drat/

You only have to do this once.

Alternatively, enter the following command each time you want to install or update the package.

```R
options(repos = c(getOption('repos'), 'https://hugheylab.github.io/drat/'))
```

Now you can install the package.

```R
setRepositories(ind = c(1:5, 9))
install.packages('metapredict', type = 'source')
```

You can update the package using `update.packages()`.

There's also a pre-built [docker image](https://hub.docker.com/r/hugheylab/hugheyverse), which has all dependencies installed.

```bash
docker pull hugheylab/hugheyverse
```

## Usage

Go through the vignettes to

- [prepare the data for the example meta-analysis](articles/prepare_example.html),
- [run the example meta-analysis](articles/run_example.html), and
- [see how to run your own meta-analysis](articles/guide.html).

For help on using specific functions, consult the [reference documentation](reference/index.html).
