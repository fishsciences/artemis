setClass("eDNA_simulation",
         slots = c(ln_conc = "matrix", Cq_star = "matrix",
                   formula = "formula", variable_levels = "list",
                   betas = "numeric", x = "data.frame",
                   std_curve_alpha = "numeric", std_curve_beta = "numeric",
                   upper_Cq = "numeric"))

setClass("eDNA_simulation_lmer", contains = "eDNA_simulation",
         slots = c(groups = "data.frame", random_sd = "numeric"))

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


setMethod("print", "eDNA_simulation",
          function(x, FUN = summary, digits = 3) {
              cat("\nformula: "); print(x@formula)
              cat("\nvariable levels:\n")
              
              vars = sapply(seq_along(x@variable_levels), function(i)
                  paste0(names(x@variable_levels)[i],  " :",
                         paste(x@variable_levels[[i]], collapse = " ")))

              cat(vars, sep = "\n")
              cat("\n ln concentration: \n")
              print(get_marginals(x@ln_conc, x@x, FUN, digits))

              cat("\n simulated Cq: \n")
              print(get_marginals(x@Cq_star, x@x, FUN, digits))
              
              return(invisible(NULL))
          })

get_marginals = function(y, X, fun, n_digits = getOption("digits"), ...)
{
    ans = lapply(X, function(x) tapply(y, x, function(dd) round(fun(dd, ...), n_digits)))
    if(all(sapply(ans, is.list)))
        ans = lapply(ans, function(x) do.call(rbind, x))

    ans       
}

setMethod("summary", "eDNA_simulation",
          function(object, ...) {
              qtl = get_marginals(object@Cq_star, object@x, quantile)
              dt = get_marginals(object@Cq_star, object@x, p_detect,
                                 thresh = object@upper_Cq)
              mapply(function(a, b) cbind(a,p_detect = b), qtl, dt)
          })


p_detect = function(y, thresh)
{
    sum(y <= thresh) / length(y)
}

setAs("eDNA_simulation", "data.frame",
      function(from){
          if(dim(from@ln_conc)[1] != 1)
              stop("Only single simulations can be converted at this time")
          
          to = cbind(data.frame(ln_conc = as.vector(from@ln_conc),
                                Cq_star = as.vector(from@Cq_star)),
                     from@x)
          to
      })


