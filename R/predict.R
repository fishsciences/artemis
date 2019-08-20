##' Predict values for eDNA model
##'
##' Predict methods.
##' 
##' @title Predict eDNA model
##' 
##' @param object an object of class eDNA_lm{er}
##' @param newdata optional, data.frame of new observations to predict
##'     values for
##' @param include_sigma logical, should the predictions include
##'     measurement error?
##' @param interval logical, should the raw predictions be returned
##'     (\code{interval = FALSE}) or should an interval be computed
##'     (\code{interval = TRUE})
##' @param interval_fun a function which computes an interval given a
##'     vector of posterior samples
##' @param ... additional arguments passed to the interval function
##' @return either a vector of predictions, or a matrix with the
##'     prediction plus interval, depending on the value of
##'     \code{interval}
##' @author Matt Espe
##' @method predict eDNA_model
##' @export
predict.eDNA_model = function(object, newdata = NULL, include_sigma = FALSE,
                              interval = FALSE,
                              interval_fun = posterior_interval,
                              ...)
{

    
    if(!is.null(newdata)){
        if(ncol(object@x) != ncol(object@betas))
            stop("Please provide the same number of predictors as the original data, including the intercept")
        X = newdata
    } else {
        X = object@x
    }
    
    ln_conc = apply(object@betas, 1, function(x) as.matrix(X) %*% x)
    Cq_hat = object@std_curve_alpha + object@std_curve_beta * ln_conc
    if(include_sigma) {
        Cq_star = sapply(seq(ncol(Cq_hat)), function(i)
            Cq_hat[,i] + rnorm(nrow(Cq_hat)) * object@sigma_Cq[i])
        Cq_star[Cq_star > object@upper_Cq] = object@upper_Cq
    } else {
        Cq_star = NULL
    }
    
    ans = list(ln_conc = ln_conc,
               Cq_hat = Cq_hat,
               Cq_star = Cq_star)

    if(interval){
        isNull = sapply(ans, is.null)
        ans[!isNull] = lapply(ans[!isNull], function(x) interval_fun(t(x), ...))
    }
    structure(ans,
              interval = interval,
              class = c("eDNA_predict_lm", class(ans)))
}

##' @param type the type of the prediction, MORE HERE
##' @rdname predict.eDNA_model
##' @export
predict.eDNA_model_lmer = function(object, newdata = NULL, 
                                   type = c(), ...)
{
    ans = NextMethod()
    ## Do something with the random effects here
    
    return(ans)
}
