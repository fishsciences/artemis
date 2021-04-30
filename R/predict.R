##' Predict values for eDNA model
##'
##' Predict methods for an eDNA model fit. Currently, these functions
##' are quite limited.
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
        
        ## if(ncol(newdata) != ncol(object@betas))
            ## stop("Please provide the same number of predictors as the original data, including the intercept")
        X = form_newdata(object, newdata)
    } else {
        X = object@x
    }
    inter = if(length(object@intercept)) as.vector(object@intercept) else 0
    
    ln_conc = apply(object@betas, 1, function(x) (as.matrix(X) %*% x)) + inter
    Cq_hat = object@std_curve_alpha + object@std_curve_beta * ln_conc
    if(include_sigma) {
        
        ln_conc_tmp = ln_conc + rnorm(nrow(ln_conc), 0, object@sigma_ln_eDNA)
        Cq_star = object@std_curve_alpha + object@std_curve_beta * ln_conc_tmp
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

form_newdata = function(object, newdata, fm = object@formula)
{
    fm = formula(delete.response(terms(fm)))
    fm = nobars(fm)
    model.frame(fm, newdata)
}

is_rand_eff = function(fm)
{
    sapply(fm, function(x) any(grep("\\|", as.character(x))))
}
