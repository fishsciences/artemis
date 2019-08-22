################################################################################
## simulations

##' eDNA simulation results
##'
##' Placeholder
##' 
##' @slot ln_conc matrix, simulated log(eDNA concentration)
##' @slot Cq_star matrx, simulated CQ star
##' @slot formula the formula used
##' @slot variable_levels the variable levels used
##' @slot betas the effect levels used
##' @slot x the model matrix used 
##' @slot std_curve_alpha the alpha for the std. curve conversion formual used
##' @slot std_curve_beta the beta for the std. curve conversion formula used
##' @slot upper_Cq the upper limit for CQ
##' 
##' @export
setClass("eDNA_simulation", 
         slots = c(ln_conc = "matrix", Cq_star = "matrix",
                   formula = "formula", variable_levels = "list",
                   betas = "numeric", x = "data.frame",
                   std_curve_alpha = "numeric", std_curve_beta = "numeric",
                   upper_Cq = "numeric"))

##' eDNA simulation results
##'
##' Placeholder
##' 
##' @slot groups the grouping variables used
##' @slot random_sd the stdev of the random effects
##' 
##' @rdname eDNA_simulation
##' @export
setClass("eDNA_simulation_lmer", contains = "eDNA_simulation",
         slots = c(groups = "data.frame", random_sd = "numeric"))

##' eDNA simulation results
##'
##' Placeholder
##' 
##' @rdname eDNA_simulation
##' @export
setClass("eDNA_simulation_lm", contains = "eDNA_simulation")

setAs("stanfit", "eDNA_simulation_lmer", function(from) callNextMethod())
setAs("stanfit", "eDNA_simulation_lm", function(from) callNextMethod())

setAs("stanfit", "eDNA_simulation",
      function(from){
          tmp = extract(from)
          new("eDNA_simulation", 
              ln_conc = tmp$ln_conc,
              Cq_star = tmp$Cq_star)
       
      })



setAs("eDNA_simulation", "data.frame",
      function(from){
          if(dim(from@ln_conc)[1] != 1)
              stop("Only single simulations can be converted at this time")
          
          to = cbind(data.frame(ln_conc = as.vector(from@ln_conc),
                                Cq_star = as.vector(from@Cq_star)),
                     from@x)
          to
      })

as.data.frame.eDNA_simulation = function(x) {
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
##' Placeholder.
##' 
##' @slot ln_conc matrix, the estimated latent variable, log(eDNA concentration)
##' @slot Cq_star matrix, the predicted CQ value for each obs
##' @slot betas matrix, the posterior estimate for each beta in the model
##' @slot sigma_Cq the estimated measurement error on CQ
##' @slot formula the formula used to fit the model
##' @slot x the model.matrix used to fit the model
##' @slot std_curve_alpha the alpha for the std. curve conversion formual used
##' @slot std_curve_beta the beta for the std. curve conversion formula used
##' @slot upper_Cq the upper limit for CQ
##' @slot stanfit the result from rstan::sampling
##'
##' @export
setClass("eDNA_model",
         slots = c(ln_conc = "matrix", Cq_star = "matrix",
                   betas = "array", sigma_Cq = "array",
                   formula = "formula", x = "data.frame",
                   std_curve_alpha = "numeric", std_curve_beta = "numeric",
                   upper_Cq = "numeric",
                   stanfit = "stanfit"))
##' eDNA simulation results
##'
##' Placeholder
##' 
##' @rdname eDNA_model
##' @export
setClass("eDNA_model_lm", contains = "eDNA_model")

##' eDNA simulation results
##'
##' Placeholder
##' 
##' @slot random_x data.frame of the grouping variables used
##' @slot random_sd the estimated stdev. of each of the random effects
##' @rdname eDNA_model
##' @export
setClass("eDNA_model_lmer", contains = "eDNA_model_lm",
         slots = c(random_x = "data.frame", random_sd = "array"))


setAs("stanfit", "eDNA_model_lmer", function(from) callNextMethod())
setAs("stanfit", "eDNA_model_lm", function(from) callNextMethod())

setAs("stanfit", "eDNA_model",
      function(from){
          tmp = extract(from)
          new("eDNA_model", 
              betas = tmp$betas,
              sigma_Cq = tmp$sigma_Cq,
              stanfit = from)
       
      })

