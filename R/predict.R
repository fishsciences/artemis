##' Predict values for eDNA model
##'
##' Predict methods for an eDNA model fit. Currently, these functions
##' are quite limited.
##' 
##' @title Predict eDNA model
##' 
##' @param object an object of class eDNA_lm\{er\}
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
##' @return a list with elements \code{ln_conc}, \code{Cq_hat}, and
##'     \code{Cq_star}. Each is a matrix with one row per observation.
##'     With \code{interval = FALSE}, each column is one posterior
##'     draw. With \code{interval = TRUE}, the columns are the output
##'     of \code{interval_fun}. \code{Cq_star} is \code{NULL} unless
##'     \code{include_sigma = TRUE}.
##' @author Matt Espe
##' @method predict eDNA_model
##' @export
predict.eDNA_model = function(object, newdata = NULL, include_sigma = FALSE,
                              interval = FALSE,
                              interval_fun = function(x) quantile(x, probs = c(0.05, 0.95)),
                              ...)
{
    X = if(is.null(newdata)) object@x else form_newdata(object, newdata)
    alpha = object@std_curve_alpha
    beta = object@std_curve_beta
    if(!is.null(newdata)){
        alpha = unique(alpha)
        beta = unique(beta)
        if(length(alpha) > 1 || length(beta) > 1)
            stop("Model was fit with multiple standard curves; predictions for newdata are not supported")
    }

    # one row per observation, one column per draw
    ln_conc = lin_pred(as.matrix(X), object@betas, object@intercept)
    Cq_hat = alpha + beta * ln_conc
    if(include_sigma) {
        sigma = rep(as.vector(object@sigma_ln_eDNA), each = nrow(ln_conc))
        ln_conc_tmp = ln_conc + rnorm(length(ln_conc), 0, sigma)
        Cq_star = alpha + beta * ln_conc_tmp
        Cq_star[Cq_star > object@upper_Cq] = object@upper_Cq
    } else {
        Cq_star = NULL
    }

    ans = list(ln_conc = ln_conc,
               Cq_hat = Cq_hat,
               Cq_star = Cq_star)

    if(interval){
        isNull = sapply(ans, is.null)
        ans[!isNull] = lapply(ans[!isNull], function(x) t(apply(x, 1, interval_fun, ...)))
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

lin_pred = function(X, betas, intercept)
{
    n_draw = max(NROW(betas), NROW(intercept))
    b = matrix(betas, n_draw, ncol(X))
    inter = if(length(intercept)) as.vector(intercept) else 0
    sweep(X %*% t(b), 2, inter, "+")
}

form_newdata = function(object, newdata, fm = object@formula)
{
    X = remove_intercept(model.matrix(delete.response(terms(nobars(fm))), newdata))
    if(!identical(colnames(X), colnames(object@x)))
        stop("Predictors in newdata do not match the model: expected ",
             paste(colnames(object@x), collapse = ", "))
    X
}

is_rand_eff = function(fm)
{
    sapply(fm, function(x) any(grep("\\|", as.character(x))))
}
