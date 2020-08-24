# artemis: an R package for eDNA analysis  

![artemis logo](man/figures/logo.png)


#### 

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/fishsciences/artemis.svg?branch=main)](https://travis-ci.org/fishsciences/artemis)
[![Codecov test coverage](https://codecov.io/gh/fishsciences/artemis/branch/main/graph/badge.svg)](https://codecov.io/gh/fishsciences/artemis?branch=main)

<!-- badges: end -->

The `artemis` package was created to aid in the design and analysis of
eDNA survey studies by offering a custom suite of models for eDNA
sampling and qPCR data. It implements a set of Bayesian
latent-variable, truncated data models which are fit using
[Stan](https://mc-stan.org/). 

## Installation

`artemis` can be installed on most platforms using one of the pre-built
binaries available from the [releases
page](https://github.com/fishsciences/artemis/releases), with more detailed installation instructions found [here](https://fishsciences.github.io/artemis/articles/artemis-installation-guide.html).

To install from github using devtools:

```
devtools::install_github("fishsciences/artemis", ref = "main")
```

## Basic use

Please refer to the [Getting Started with the `artemis` package](https://fishsciences.github.io/artemis/articles/artemis-overview.html) vignette, which covers most of the functionality of artemis.

Additional vignettes are forthcoming!


## Reporting bugs

Please report all bugs via an issue at the package
[repo](https://github.com/fishsciences/artemis/issues).

