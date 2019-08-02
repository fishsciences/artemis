setMethod("plot", "eDNA_simulation",
          function(x, y, ...) {})

setMethod("plot", "eDNA_model",
          function(x, y, ...) {})

## Might do something slightly different/expanded with lmer sims and fits
setMethod("plot", "eDNA_simulation_lmer",
          function(x, y, ...) {})

setMethod("plot", "eDNA_model_lmer",
          function(x, y, ...) {})

plot.eDNA_p_detect = function(x, y, probs = c(0.025, 0.975),
                              ylim = c(0,1), ...)
{
    if(length(probs) != 2)
        stop("Please provide a lower and upper bound in probs")
    
    reps = attr(x, "reps")

    if(!is.null(dim(x))) {
        mn = apply(x, 2, mean)
    
        ci = apply(x, 2, quantile, probs)
        plot(reps, mn, ylab = "p(detect)", type = "b",
             ylim = ylim, xlim = range(reps), ...)

        segments(reps, y0 = ci[1,], y1 = ci[2,])
        lines(reps, ci[1,], lty = 2)
        lines(reps, ci[2,], lty = 2)
    } else {
        plot(reps, x, ylab = "p(detect)", ylim = ylim, ...)
    }
    

    invisible(x)
}
