##' Plot method for eDNA simulations
##'
##' Plot method for eDNA simulations
##' @title Plot method for eDNA simulations
##' @param x object of class eDNA_simulation
##' @param y ignored
##' @param response the response variable to plot
##' @param probs the probability for plotting CIs
##' @param ... ignored
##' @return a gtable object from marrangeGrob
##' @author Matt Espe
##' @method plot eDNA_simulation
##' @export
plot.eDNA_simulation = function(x, y,
                                response = "Cq_star",
                                probs = c(0.025, 0.975),
                                ...)
{
    y = slot(x, response)
    if(dim(y)[1] == 1){
        ymn = as.vector(y)
        yci = rep(0, length(ymn))
    } else {
        ymn = apply(y, 2, mean)
        yci = apply(y, 2, quantile, probs)
    }
    vars = colnames(x@x)
    p = lapply(vars, function(v){
        ggplot(, aes(y = ymn, x = x@x[[v]] )) +
            geom_point() +
            geom_smooth() +
            xlab(v) +
            ylab(response) +
            theme_bw()
    })
    marrangeGrob(p, ncol = 3, nrow = 1)
}


##' Plot eDNA model results
##'
##' Plot eDNA model results. Currently, this is just a wrapper for
##' \code{rstan::plot}, which produces a "catapillar" plot.
##' 
##' @title Plot eDNA model results
##' @param x object of class eDNA_model_*
##' @param y ignored
##' @param pars parameters to plot
##' @param ... additional args passed to plot.stanfit
##' @return ggplot of posterior estimates
##' @author Matt Espe
##' @method plot eDNA_model
##' @export
plot.eDNA_model = function(x, y, pars = "betas",  ...) {
    rstan::plot(x@stanfit, pars = pars, ...)
}



##' Plot the p(detect)
##'
##' Plot the p(detect)
##' @title Plot p(detect)
##' @param x an object of class "eDNA_p_detect", produced by \code{est_p_detect}
##' @param y ignored
##' @param probs numeric vector of length 2, the lower and upper probabilities for
##'     quantiles displayed
##' @param ylim the y limits on the plot
##' @param ... additional args passed to \code{plot}
##' @return NULL
##' 
##' @author Matt Espe
##' @method plot eDNA_p_detect
##' @export
plot.eDNA_p_detect = function(x, y, probs = c(0.025, 0.975),
                              ylim = c(0,1), ...)
{
    if(length(probs) != 2)
        stop("Please provide a lower and upper bound in probs")
    
    reps = attr(x, "reps")

    if(is.matrix(x)) {
        tmp = summary(x, probs)
        plot(tmp$n_rep, tmp$mean,
             ylab = "p(detect)", xlab = "N replicates",
             type = "b",
             ylim = ylim, xlim = range(reps),
             lty = 2, ...)

        segments(tmp$n_rep, y0 = tmp[, 3], y1 = tmp[, 4])
        lines(tmp$n_rep, tmp[, 3], lty = 3)
        lines(tmp$n_rep, tmp[, 4], lty = 3)
    } else {
        plot(reps, x, ylab = "p(detect)", ylim = ylim, ...)
    }
    
    invisible(x)
}
