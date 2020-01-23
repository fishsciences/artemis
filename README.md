# artemis: an R package for eDNA analysis  

![artemis logo](man/figures/logo.png)

#### 

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/fishsciences/artemis.svg?branch=master)](https://travis-ci.org/fishsciences/artemis)
[![Codecov test coverage](https://codecov.io/gh/fishsciences/artemis/branch/master/graph/badge.svg)](https://codecov.io/gh/fishsciences/artemis?branch=master)

<!-- badges: end -->


The `artemis` package was created to aid in the design and analysis of
eDNA survey studies by offering a custom suite of models for eDNA
sampling and qPCR data. It implements a set of Bayesian
latent-variable, truncated data models which are fit using
[Stan](https://mc-stan.org/). 

## Installation

Currently, artemis can be installed using one of the pre-built
binaries available from the [releases
page](https://github.com/fishsciences/artemis/releases). Alternatively,
you can install the package from github using the `devtools` package,

```
devtools::install_github("fishsciences/artemis", build_vignettes = TRUE)
```

Please note, installing the package from github requires a C++
compiler (e.g. gcc, clang).

## Basic use

Please refer to the `artemis-overview` vignette, which covers most of the functionality of artemis,

```
vignette("artemis-overview", package = "artemis")
```

Additional vignettes are forthcoming!


## Reporting bugs

Please report all bugs via an issue at the package
[repo](https://github.com/fishsciences/artemis/issues).

