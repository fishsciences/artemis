##' Plot method for eDNA simulations
##'
##' Plot method for eDNA simulations
##' @title Plot method for eDNA simulations
##' @param x object of class eDNA_simulation
##' @param y ignored
##' @param response the response variable to plot
##' @param probs the probability for plotting CIs
##' @param alpha the alpha value, i.e. transparancy, of the points
##' @param jitter_width the width of the jitter applied to the points
##' @param ... ignored
##' @return a list of ggplots
##' @author Matt Espe
##' @method plot eDNA_simulation
##' @export
plot.eDNA_simulation = function(x, y,
                                response = "Cq_star",
                                probs = c(0.025, 0.975),
                                alpha = 0.1, jitter_width = 0.35, 
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
            geom_jitter(alpha = alpha, width = jitter_width) +
            xlab(v) +
            ylab(response) +
            theme_bw()
    })
    p
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
    plot(x@stanfit, pars = pars, ...)
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
##' @param point_size point size for the estimate or the mean value
##' @param n_breaks passed to \code{pretty()}
##' @param ... additional args, currently ignored
##' @return a ggplot object
##' 
##' @author Matt Espe
##' @method plot eDNA_p_detect
##' @export
plot.eDNA_p_detect = function(x, y, probs = c(0.025, 0.975),
                              ylim = c(0,1),
                              point_size = 2.5,
                              n_breaks = 5,
                              ...)
{
    if(length(probs) != 2)
        stop("Please provide a lower and upper bound in probs")
    
    reps = attr(x, "reps")

    if(is.matrix(x)) {
        tmp = summary(x, probs)
        p = ggplot(, aes(x = as.integer(tmp$n_rep), y = tmp$mean)) + #avoids R CMD check NOTE
            geom_pointrange(aes(ymin = tmp[,3], # lower
                                ymax = tmp[,4]), # upper
                            size = point_size) +
            geom_line(lty = 3) 
    } else {
        p = ggplot( , aes(x = as.integer(reps), y = x)) +
            geom_point(size = point_size) 
    }
        p + ylab("p(detect)") + xlab("N replicates") +
            scale_x_continuous(breaks = break_fun(as.integer(reps) , n = n_breaks),
                expand = expand_scale(mult = 0.015)) +
            scale_y_continuous(limits = ylim, expand = expand_scale(add = 0.01)) +
            theme_bw()
}

break_fun = function(x, n) {
    pretty(x, n)[round(pretty(x, n),1) %% 1 == 0]
}
