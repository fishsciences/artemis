## From https://stackoverflow.com/questions/72559243/is-there-a-way-method-to-assign-a-r6-object-to-an-s4-object-slot
## Required to put an R6 object inside an S4 object
setClass("CmdStanMCMC_CSV")


################################################################################
## simulations

##' eDNA simulation results
##'
##' An S4 object holding the results of \code{sim_eDNA_lm*}.
##' 
##' @slot ln_conc_hat matrix, simulated log(eDNA concentration)
##'     including measurement error. One row per observation, one
##'     column per simulation.
##' @slot ln_conc_star matrix, expected log(eDNA concentration) from
##'     the predictors and random effects
##' @slot Cq_star matrix, simulated Cq values after censoring at
##'     \code{upper_Cq} and zero inflation
##' @slot formula the formula used
##' @slot variable_levels the variable levels used
##' @slot betas the effect levels used
##' @slot x the model matrix used 
##' @slot std_curve_alpha the alpha for the std. curve conversion formula used
##' @slot std_curve_beta the beta for the std. curve conversion formula used
##' @slot upper_Cq the upper limit for CQ
##' 
##' @export
setClass("eDNA_simulation", 
         slots = c(ln_conc_hat = "matrix",ln_conc_star = "matrix", Cq_star = "matrix",
                   formula = "formula", variable_levels = "list",
                   betas = "numeric", x = "data.frame",
                   std_curve_alpha = "numeric", std_curve_beta = "numeric",
                   upper_Cq = "numeric"))

##' @rdname eDNA_simulation-class
##' @export
setClass("eDNA_simulation_lmer", contains = "eDNA_simulation")

##' @rdname eDNA_simulation-class
##' @export
setClass("eDNA_simulation_lm", contains = "eDNA_simulation")
setAs("CmdStanMCMC_CSV", "eDNA_simulation_lmer", function(from) callNextMethod())
setAs("CmdStanMCMC_CSV", "eDNA_simulation_lm", function(from) callNextMethod())

setAs("CmdStanMCMC_CSV", "eDNA_simulation", function(from){
  tmp = unclass(from$draws(format = "draws_matrix"))

  new("eDNA_simulation",
      ln_conc_hat = sim_draws(tmp, "ln_conc_hat"),
      ln_conc_star = sim_draws(tmp, "ln_conc_star"),
      Cq_star = sim_draws(tmp, "Cq_star"))

})

sim_draws = function(draws, var)
  # one row per observation, one column per simulation
{
  i = startsWith(colnames(draws), paste0(var, "["))
  t(draws[, i, drop = FALSE])
}

setAs("eDNA_simulation", "data.frame",
      function(from){
          if(dim(from@ln_conc_hat)[2] != 1)
              stop("Only single simulations can be converted at this time")
          
          to = cbind(data.frame(ln_conc_hat = as.vector(from@ln_conc_hat),
                                Cq_star = as.vector(from@Cq_star)),
                     from@x)
          to
      })

##' Methods for eDNA simulations
##'
##' as.data.frame methods for eDNA simulations. This allows the
##' conversion of the simulations to a form suitable for additional
##' operations, e.g. plotting.
##' 
##' @title  Methods for eDNA simulations
##' @param x object of class eDNA_simulation
##' @param row.names ignored
##' @param optional ignored
##' @param ... ignored
##' @return data.frame 
##' @author Matt Espe
##' @method as.data.frame eDNA_simulation
##' @export 
as.data.frame.eDNA_simulation = function(x, row.names, optional, ...) {
    as(x, "data.frame")
}

get_marginals = function(y, X, fun, ...)
{
    ans = lapply(X, function(x) tapply(y, x, fun, ...))
    if(all(sapply(ans, is.list)))
        ans = lapply(ans, function(x) do.call(rbind, x))

    ans       
}

p_detect = function(y, thresh)
{
    sum(y < thresh) / length(y)
}

################################################################################
## Model results

##' eDNA model fit results
##'
##' An S4 object holding the results of a model fit by \code{eDNA_lm*}
##' 
##' @slot intercept matrix, the posterior draws of the intercept, if
##'     the model has one
##' @slot betas matrix, the posterior estimate for each beta in the model
##' @slot sigma_ln_eDNA the estimated measurement error on ln_eDNA
##' @slot formula the formula used to fit the model
##' @slot x the model.matrix used to fit the model
##' @slot xz the model.matrix for the zero-inflation component, for
##'     zero-inflated models
##' @slot std_curve_alpha the alpha for the std. curve conversion formula used
##' @slot std_curve_beta the beta for the std. curve conversion formula used
##' @slot upper_Cq the upper limit for CQ
##' @slot random_names names of the random effects, one per level
##' @slot random_groups names of the random effect groups, one per
##'     estimated standard deviation
##' @slot fit the result from cmdstanr::sample
##'
##' @export
setClass("eDNA_model",
         slots = c(intercept = "array",
                   betas = "array",
                   sigma_ln_eDNA = "array",
                   formula = "formula",
                   x = "data.frame",
                   xz = "data.frame", # for zero-infl
                   random_names = "character", # for mixed models
                   random_groups = "character",
                   std_curve_alpha = "numeric", std_curve_beta = "numeric",
                   upper_Cq = "numeric",
                   fit = "ANY"))

##' @rdname eDNA_model-class
##' @export
setClass("eDNA_model_lm", contains = "eDNA_model")

##' @rdname eDNA_model-class
##' @export
setClass("eDNA_model_zip", contains = "eDNA_model_lm")

##' @rdname eDNA_model-class
##' @export
setClass("eDNA_model_zipr", contains = "eDNA_model_zip")

##' @rdname eDNA_model-class
##' @export
setClass("eDNA_model_count", contains = "eDNA_model")

##' @rdname eDNA_model-class
##' @export
setClass("eDNA_model_countr", contains = "eDNA_model_count")


setAs("eDNA_model_zip", "eDNA_model", function(from){
  from
})

setAs("eDNA_model_zipr", "eDNA_model", function(from){
  from
})

setAs("eDNA_model_count", "eDNA_model", function(from){
  from
})

setAs("eDNA_model_countr", "eDNA_model", function(from){
  from
})

##' @rdname eDNA_model-class
##' @export
setClass("eDNA_model_lmer", contains = "eDNA_model_lm")

## models with random effects, for ranef() dispatch
setClassUnion("eDNA_model_mixed",
              c("eDNA_model_lmer", "eDNA_model_zipr", "eDNA_model_countr"))

setAs("eDNA_model_lmer", "eDNA_model", function(from) {
  
  ## fix later
  from
})

setAs("CmdStanMCMC_CSV", "eDNA_model",
      function(from){
        tmp = as.data.frame(from$draws(format = "draws_df"))
        nms = names(tmp)
        i = grepl("^betas\\[", nms)
        betas = if(any(i)) as.matrix(tmp[i]) else array(numeric())

        j = grepl("intercept", nms)
        intercept = if(any(j)) as.matrix(tmp[j]) else array(numeric())
        
        k = grepl("sigma_ln_eDNA", nms)
        sigma = if(any(k)) as.matrix(tmp[k]) else array(numeric())
        
        new("eDNA_model",
            intercept = intercept,
            betas = betas,
            sigma_ln_eDNA = sigma,
            fit = from)
        
      })

