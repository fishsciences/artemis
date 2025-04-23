# artemis: ChangeLog

# v3.0

## 3.0.0

Updates the models in `artemis` to include a Poisson distributed GLM
(count model) and zero-inflated models. Added to address zero-inflated
data.

# v2.0

## 2.0.1

Fixed bug when data has no censored values. Corrected the behavior of `compile_models()`.

## 2.0.0

Complete overhaul to use the R package `cmdstanr` as the backend for fitting models, rather than `rstan`. 

# v1.0

## 1.0.6

Fixes missing import for rstan::plot. Moves/removes vignette material
to make the vignettes more streamlined.

## 1.0.5

Addresses review comments from CRAN.

## 1.0.4

Fixed issue with predict methods when random effects are present in
the formula.

## 1.0.0-3

Various fixes to reduce NOTES for CRAN submission.

# Beta

WARNING: Any versions here are still under active development, and
will change in ways that will likely break your existing code. Until
formal release v1.0.0, do not expect stability!

## v0.18.1

Adds check to prevent running model functions with NA values in the
response. Normally, these values are dropped silently by `model.frame`
or `lFormula`, but for our purposes this can create issues due to
mismatches between the standard curve parameters and the model matrix.
We decide to force conscious removal by the user. 

## v0.18.0

Changes the eDNA_data and data documentation to match what is in Espe
et al. 2021 (in press). Data has completely changed; column names,
where applicable, are the same. Year column is no longer present.

## v0.17.0

Moves the main development branch back to using the `rstan` backend
for running the models. This was done because `cmdstanr` is not yet on
CRAN, and we cannot have a strong dependency on a github package if
`artemis` is to be on CRAN.

## v0.16.0

Changed model and functions to apply the measurement error on ln(eDNA)
rather than on Cq. This is because when there are multiple standard
curves in use, effectively there will also be multiple sigma_Cqs. This
is atypical for a regression model, where the residual error is
assumed to be I.I.D. Changing the error to apply to the ln(eDNA)
avoids this, as the observations are on the same scale.

## v0.15.0

Refactored the model code to make increase clarity and
readability. This might reduce model performance slightly, but should
make the code easier to maintain.

Separated the zero-inflated models into their own Stan models and
associated R functions. This should will force users to explicitely
specify that they want the zero-inflated model. When this model does
not fit well, it should be clearer.

## v0.14.0

## v0.13.0

Changed zero-inflated probability to a user-supplied value rather than
an estimate in the model. This is because a large amount of data is
required to reliably estimate the probability of zeros. In situations
with lower amounts of data, the model was failing to fit properly.

## v0.12.0

Added a zero-inflated probability to the model after observing that
even with fairly consistent concentrations of eDNA, filters sometimes
experience a zero value. 

## v0.11.0

Added the ability to provide multiple standard curves, which is useful
when analyzing data from multiple experiments in a single model.

## v0.10.0

Added the ability to allow measurement error to vary across different
Cq values. This was driven by the observation that the measurement
error tended to increase as Cq increased (i.e., fewer particles of
eDNA to detect). This is enabled in the modeling functions by setting
`Cq_error_type = varying`. By default, the error is assumed to be
stationary, i.e. fixed for all values of Cq.

WARNING: In general, more data is required to reliably estimate a
varying measurement error. Enabling this with small datasets is likely
to result in unreliable estimates. 

## v0.9.5

Refactoring of how the package handled intercepts. Previously,
intercepts were just part of the model matrix. While conceptually
simple and easy to implement, this arrangement created issues when
trying to add other features, primarily priors over the betas such as
the Regularized Hierarchical Shrinkage prior.

- Intercepts have priors applied separately from other predictors.
- Each column of the model matrix is centered internally for more
  efficient computations
- By default, weakly informative priors are put on all of the
  "betas". We followed the general recommendations of the developers
  of
  [Stan/rstanarm](https://cran.r-project.org/web/packages/rstanarm/vignettes/priors.html).
  There is a `normal()` function which facilitates specifying these priors.
- Priors are no longer specified by providing the `mu` and `sd`
  values. Rather, similar to rstanarm, the priors are now specified
  using a named list. This will help us implement more prior
  distributions in the future.
- summary and predict methods automatically detect the intercept and handle it if present.

## v0.9.0

First "public" beta.
