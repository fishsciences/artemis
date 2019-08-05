################################################################################
## simulations

##' @export
setClass("eDNA_simulation",
         slots = c(ln_conc = "matrix", Cq_star = "matrix",
                   formula = "formula", variable_levels = "list",
                   betas = "numeric", x = "data.frame",
                   std_curve_alpha = "numeric", std_curve_beta = "numeric",
                   upper_Cq = "numeric"))

##' @export
setClass("eDNA_simulation_lmer", contains = "eDNA_simulation",
         slots = c(groups = "data.frame", random_sd = "numeric"))

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

##' @export
setClass("eDNA_model",
         slots = c(ln_conc = "matrix", Cq_star = "matrix",
                   betas = "array", sigma_Cq = "array",
                   formula = "formula", x = "data.frame",
                   std_curve_alpha = "numeric", std_curve_beta = "numeric",
                   upper_Cq = "numeric",
                   stanfit = "stanfit"))

##' @export
setClass("eDNA_model_lmer", contains = "eDNA_model",
         slots = c(groups = "data.frame", random_sd = "numeric"))

##' @export
setClass("eDNA_model_lm", contains = "eDNA_model")

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
