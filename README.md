# artemis: an R package for eDNA analysis  

![artemis logo](man/figures/logo.png)

#### 

The `artemis` package was created to aid in the design and analysis of
eDNA survey studies by offering a custom suite of models for eDNA
sampling and qPCR data. It implements a set of Bayesian
latent-variable, truncated data models which are fit using
[Stan](mc-stan.org). 

## Installation

Currently, artemis can be installed using one of the pre-built
binaries available from the [releases
page](https://github.com/fishsciences/artemis/releases). Alternatively,
you can install the package from github using the `devtools` package,

```
devtools::install_github("fishsciences/artemis")
```

Please note, installing the package from github requires a C++
compiler (e.g. gcc, clang).

## Basic use

Please refer to the `artemis-overview` vignette, which covers most of the
major functionality of artemis,

```
vignette("artemis-overview", package = "artemis")
```

## Reporting bugs

Please report all bugs via an issue at the package
[repo](https://github.com/fishsciences/artemis/issues).

