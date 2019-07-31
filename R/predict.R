
predict.eDNA_model = function(object, newdata = NULL, include_sigma = FALSE,
                              interval = FALSE, ...)
{

    
    if(!is.null(newdata)){
        if(ncol(x) != ncol(object@betas))
            stop("Please provide the same number of predictors as the original data, including the intercept")
        X = newdata
    } else {
        X = object@x
    }
    
    ln_conc = apply(object@betas, 1, function(x) as.matrix(X) %*% x)
    Cq_hat = object@std_curve_alpha + object@std_curve_beta * ln_conc
    Cq_star = if(include_sigma) {
                  sapply(seq(ncol(Cq_hat)), function(i)
                      Cq_hat[,i] + rnorm(nrow(Cq_hat)) * object@sigma_Cq[i])
    } else {
        NULL
    }
    ans = list(ln_conc = ln_conc,
               Cq_hat = Cq_hat,
               Cq_star = Cq_star)

    if(interval){
        isNull = sapply(ans, is.null)
        ans[!isNull] = lapply(ans[!isNull], posterior_interval, ...)
    }
    structure(ans,
              class = "eDNA_predict_lm")
}

predict.eDNA_model_lmer = function(object, newdata = NULL, 
                                   type = c(), ...)
{
    ans = NextMethod()
    ## Do something with the random effects here
    
    return(ans)
}
