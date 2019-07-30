setClass("eDNA_predict_lm",
         slots = c(ln_conc = "matrix",
                   Cq_hat = "matrix"))

setClass("eDNA_predict_lmer", contains = "eDNA_predict_lm")

setMethod("predict", "eDNA_model",
          function(object, newdata = NULL, 
                   type = c(), ...)
          {
              if(!is.null(newdata)){
                  X = newdata
              } else {
                  X = object@x
              }
              
              ln_conc = apply(object@betas, 1, function(x) as.matrix(X) %*% x)
              Cq_hat = object@std_curve_alpha + object@std_curve_beta * ln_conc
              new("eDNA_predict_lm",
                  ln_conc = ln_conc,
                  Cq_hat = Cq_hat)
          })

setMethod("predict", "eDNA_model_lmer",
          function(object, newdata = NULL, 
                   type = c(), ...)
          {
              ans = callNextMethod()
              ## Do something with the random effects here
              
              return(ans)
          })
