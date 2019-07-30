
setMethod("summary", "eDNA_simulation",
          function(object, ...) {
              qtl = get_marginals(object@Cq_star, object@x, quantile)
              dt = get_marginals(object@Cq_star, object@x, p_detect,
                                 thresh = object@upper_Cq)
              mapply(function(a, b) cbind(a,p_detect = b), qtl, dt)
          })
