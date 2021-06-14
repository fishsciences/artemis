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

<br>
The easiest way to install `artemis` is from [CRAN](https://cran.r-project.org/package=artemis) with the `install.packages` function:


```{r, eval=FALSE, include=TRUE}

install.packages("artemis")

```
<br>

### Testing your installation

If your installation of `artemis` and its dependencies was successful, the following code should run without error (although you may see warning messages from `rstan` about Bulk/Tail Effective Samples Sizes being too low). If the first or second model returns an error that seems to have something to do with your `c++` compiler, you may need to [follow instructions to edit your `Makevars` or `Makevars.win` file](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

```{r, eval=FALSE}
library(artemis)

model_fit = eDNA_lm(Cq ~ scale(Distance_m) + scale(Volume_mL), 
                    data = eDNA_data,
                    std_curve_alpha = 21.2, std_curve_beta = -1.5)

model_fit2 = eDNA_lmer(Cq ~ scale(Distance_m) + scale(Volume_mL) + (1|FilterID),
                       eDNA_data,
                       std_curve_alpha = 21.2, std_curve_beta = -1.5)

```


### Installing `artemis` from source

Installing `artemis` from source on Windows is not currently well-supported; we recommend installing from the pre-compiled binary if you're on Windows.  If you're on MacOS or Linux and you prefer to install from source, then go ahead and do that with your function/utility of choice (`devtools::install_github()`, `utils::install.packages(type = "source")`, `R CMD INSTALL`, etc.).  
<br>
If you have sub-architecture you're really in to customizing, the source code is [here](https://github.com/fishsciences/artemis), go nuts.

<br>

## Basic use

Please refer to the [Getting Started with the `artemis` package](https://fishsciences.github.io/artemis/articles/artemis-overview.html) vignette, which covers most of the functionality of artemis.

Additional vignettes are forthcoming!


## Reporting bugs

Please report all bugs via an issue at the package
[repo](https://github.com/fishsciences/artemis/issues).

