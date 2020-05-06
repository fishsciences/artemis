# artemis: ChangeLog

# v 0.9.5

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

# v 0.9.0

First "public" beta.
