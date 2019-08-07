
##' @rdname eDNA_lmer
##' @export
eDNA_lm = function(formula, data, 
                   std_curve_alpha, std_curve_beta,
                   upper_Cq = 40,
                   n_chain = 4L, iters = 1000L, verbose = FALSE, ...)
{
    
    ml = gen_model_list_lm(formula, data)
       
    # This works because Stan ignores extra input data
    md = prep_data(ml, std_curve_alpha, std_curve_beta,
                    Cq_upper = upper_Cq, type = "model")
    md$y = ml$y
    fit = run_model(data = md, n_chain = n_chain, iters = iters, verbose = verbose, ...)
    fit = load_slots_model(fit)
    class(fit) = "eDNA_model_lm"
    
    return(fit)
}

##' Fit eDNA model
##'
##' Fit eDNA model. 
##' @title Fit eDNA Model
##'
##' @aliases eDNA_lm
##' 
##' @param formula a formula, specifying the relationship between the
##'     predictors and the latent variable eDNA concentration. 
##' @param data data.frame, with the response and predictors
##' @param std_curve_alpha the alpha value for the formula for converting between
##'     log(eDNA concentration) and CQ value
##' @param std_curve_beta the beta value for the formula for converting between
##'     log(eDNA concentration) and CQ value
##' @param upper_Cq numeric, the upper limit on CQ detection. Any value of
##'     log(concentration) which would result in a value greater than this limit is
##'     instead recorded as the limit.
##' @param n_chain integer, the number of chains to use. Please note that less than
##'     two is not recommended.
##' @param iters integer, the number of iterations per chain
##' @param verbose logical, when TRUE output from \code{rstan::sampling} is written
##'     to the console. 
##' @param ... additional arguments passed to \code{\link[rstan]{sampling}}
##' @return S4 object, with the following slots:
##' \describe{
##'   \item{ln_conc}{matrix, the posterior samples for the latent variable,
##'       eDNA concentration}
##'   \item{Cq_star}{matrix, the posterior prediction for the observed response}
##'   \item{betas}{array, the posterior estimates for the betas for the linear model}
##'   \item{sigma_Cq}{array, the posterior estimates for the measurement error of CQ}
##'   \item{formula}{formula, the original formula used in the model}
##'   \item{x}{data.frame, the model matrix used in the model}
##'   \item{std_curve_alpha}{numeric, the std. curve value used}
##'   \item{std_curve_beta}{numeric, the std. curve value used}
##'   \item{upper_Cq}{numeric, the upper limit for observed CQ used}
##'   \item{stanfit}{stanfit, the original results from \code{rstan::sampling}}
##' }
##' @author Matt Espe
##'
##' @examples
##'
##' ans = eDNA_lm(Cq ~ Distance, eDNA_samples,
##'               std_curve_alpha = 21.2, std_curve_beta = -1.5)
##'
##' ans2 = eDNA_lmer(Cq ~ Distance + Volume + (1|TechnicalRep), eDNA_samples,
##'                  std_curve_alpha = 21.2, std_curve_beta = -1.5)
##'
##'
##' @export
eDNA_lmer = function(formula, data, 
                     std_curve_alpha, std_curve_beta,
                     upper_Cq = 40,
                     n_chain = 4L, iters = 1000L, verbose = FALSE, ...)
{
    
    ml = gen_model_list_lmer(formula, data)
    
    md = prep_data(ml, std_curve_alpha, std_curve_beta,
                   Cq_upper = upper_Cq, type = "model")
    md$y = ml$y

    fit = run_model(model = stanmodels$eDNA_lmer,
                    data = md, n_chain = n_chain, iters = iters, verbose = verbose, ...)
    fit = load_slots_model(fit)
    class(fit) = "eDNA_model_lmer"
    return(fit)
}

run_model = function(model = stanmodels$eDNA_lm, data, n_chain,
                  iters, verbose, ...)

{
    if(!verbose) sink(tempfile())
    
    fit = sampling(model, data, chains = n_chain,
                  iter = iters,
                  refresh = ifelse(verbose, 100, -1), show_messages = verbose, ...)

    if(!verbose) sink()
    
    fit = as(fit, "eDNA_model")
    return(fit)
}

load_slots_model = function(obj)
    {
    obj@formula = get("formula", parent.frame())
    obj@x = as.data.frame(get("md", parent.frame())$X)
    obj@std_curve_alpha = get("std_curve_alpha", parent.frame())
    obj@std_curve_beta = get("std_curve_beta", parent.frame())
    obj@upper_Cq = get("upper_Cq", parent.frame())
    
    obj
}
