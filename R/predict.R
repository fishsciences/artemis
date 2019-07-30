setMethod("predict", "eDNA_model",
          function(object, newdata = NULL, 
                   type = c(), ...)
          {
              
          })

setMethod("predict", "eDNA_model_lmer",
          function(object, newdata = NULL, 
                   type = c(), ...)
          {
              ans = callNextMethod()
              ## Do something with the random effects here

              return(ans)
          })
