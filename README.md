# artemis: an R package for eDNA analysis

The `artemis` package was created to aid in the design and analysis of
eDNA survey studies by offering a custom suite of models for eDNA
sampling and qPCR data. It implements a set of Bayesian
latent-variable, truncated data models which are fit using
[Stan](mc-stan.org). 

## Installation

Currently, artemis can be installed using one of the pre-built
binaries available from the [release
page](https://github.com/fishsciences/artemis/releases). Alternatively,
you can install the package from github using the `devtools` package,

```
devtools::install_github("fishsciences/artemis")
```

Please note, installing the package from github requires a C++
compiler (e.g. gcc).

## Basic use

Please refer to the Getting Started vignette,

```
vignette("Getting_started", package = "artemis")
```

## Reporting bugs

Please report all bugs via an issue at the package [repo](https://github.com/fishsciences/artemis).
