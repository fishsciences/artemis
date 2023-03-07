# artemis: an R package for eDNA analysis  

![artemis logo](man/figures/logo.png)


#### 

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/fishsciences/artemis/branch/main/graph/badge.svg)](https://app.codecov.io/gh/fishsciences/artemis?branch=main)

[![R-CMD-check](https://github.com/fishsciences/artemis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fishsciences/artemis/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `artemis` package was created to aid in the design and analysis of
eDNA survey studies by offering a custom suite of models for eDNA
sampling and qPCR data. It implements a set of Bayesian
latent-variable, truncated data models which are fit using
[Stan](https://mc-stan.org/). 

## Installation

<br>
The artemis package requires the `cmdstanr` package (as of 2022-03-24). This can be installed via

```
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```

Once that is installed, the easiest way to install `artemis` is with the `devtools` or `remotes` package:

```
devtools::install_github("fishsciences/artemis")

```
<br>

### Testing your installation

Before you can run the models in `artemis`, the backend `cmdstan` must be installed. This can be done with,

```
cmdstanr::install_cmdstan()
```

If your installation of `artemis` and its dependencies was successful, the following code should run without error (although you may see warning messages from `Stan` about Bulk/Tail Effective Samples Sizes being too low). 

<!-- If the first or second model returns an error that seems to have something to do with your `c++` compiler, you may need to [follow instructions to edit your `Makevars` or `Makevars.win` file](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started). -->

```
library(artemis)

model_fit = eDNA_lm(Cq ~ scale(Distance_m) + scale(Volume_mL), 
                    data = eDNA_data,
                    std_curve_alpha = 21.2, std_curve_beta = -1.5)

model_fit2 = eDNA_lmer(Cq ~ scale(Distance_m) + scale(Volume_mL) + (1|FilterID),
                       eDNA_data,
                       std_curve_alpha = 21.2, std_curve_beta = -1.5)

```

The first time these models are run the model code will need to be compiled. Thereafter, they should run without needing to be re-compiled.


## Basic use

Please refer to the [Getting Started with the `artemis` package](https://fishsciences.github.io/artemis/articles/artemis-overview.html) vignette, which covers most of the functionality of artemis.

Additional vignettes are forthcoming!


## Reporting bugs

Please report all bugs via an issue at the package
[repo](https://github.com/fishsciences/artemis/issues).

