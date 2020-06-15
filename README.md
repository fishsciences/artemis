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

## Experimental Branch!

This is an experimental branch which uses cmdstanr as the backend for
fitting the models instead of rstan. Because of this, there are
several pieces behind the scenes which are different from the master
branch. 

## Installation

`artemis` can be installed on most platforms using one of the pre-built
binaries available from the [releases
page](https://github.com/fishsciences/artemis/releases), with more detailed installation instructions found [here](https://fishsciences.github.io/artemis/articles/artemis-installation-guide.html).

## Basic use

Please refer to the [Getting Started with the `artemis` package](https://fishsciences.github.io/artemis/articles/artemis-overview.html) vignette, which covers most of the functionality of artemis.

Additional vignettes are forthcoming!


## Reporting bugs

Please report all bugs via an issue at the package
[repo](https://github.com/fishsciences/artemis/issues).

## cmdstanr backend: experimental

If you would prefer to use cmdstan as the backend for fitting the
models, there is an experimental branch, cmdstan. Using cmdstan as the
backend reduces the overhead in R associated with fitting the model,
and allows for more robust parallelization of chains.

However, this branch is experimental, so use at your own risk.

