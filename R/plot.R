##' Plot method for eDNA simulations
##'
##' Plot method for eDNA simulations
##' @title Plot method for eDNA simulations
##' @param x object of class eDNA_simulation
##' @param y ignored
##' @param response the response variable to plot
##' @param probs the probability for plotting CIs
##' @param ... ignored
##' @return a list of ggplots
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
            geom_jitter(alpha = 1/10, width = 0.35) +
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
##' @param ... additional args, currently ignored
##' @return a ggplot object
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
        p = ggplot(tmp, aes(x = as.integer(n_rep), y = tmp$mean)) + 
            geom_point(size = 2.5) +
            geom_line(lty = 3) + 
            geom_errorbar(aes(ymin = tmp[,3], # lower
                              ymax = tmp[,4]), # upper
                            width = rel(0.15)) 


     } else {
        p = ggplot( , aes(x = as.integer(reps), y = x)) +
            geom_point(size = 2.5) 
    }
        p + ylab("p(detect)") + xlab("N replicates") +
            scale_x_continuous(breaks = function(x, n = 5)  
                pretty(x, n)[round(pretty(x, n),1) %% 1 == 0],
                                   expand = expand_scale(mult = 0.015)) +
            scale_y_continuous(limits = ylim, expand = expand_scale(add = 0.01)) +

            coord_flip()  +
            theme_bw()
        
}
