
##' @rdname eDNA_lmer
##' @export
eDNA_lm = function(formula, data, 
                   std_curve_alpha, std_curve_beta,
                   upper_Cq = 40,
                   n_chain = 4L, iters = 1000L, verbose = FALSE,
                   sink_file = tempfile(),
                   betas_prior_mu = numeric(), betas_prior_sd = numeric(), ...)
{
    
    ml = gen_model_list_lm(formula, data)
       
    # This works because Stan ignores extra input data
    md = prep_data(ml, std_curve_alpha, std_curve_beta,
                   Cq_upper = upper_Cq, type = "model",
                   b_prior_mu = betas_prior_mu,
                   b_prior_sd = betas_prior_sd)
    md$y = ml$y
    fit = run_model(data = md, n_chain = n_chain, iters = iters,
                    verbose = verbose, sink_file = sink_file, ...)
    fit = load_slots_model(fit)
    fit = as(fit, "eDNA_model_lm")
    
    return(fit)
}

##' Fit eDNA model
##'
##' These functions fit a Bayesian latent variable model to data
##' collected from a eDNA sampling experiment. These data have a few
##' particular characteristics that justify using a specialized
##' model. More details on these characteristics and the model
##' structure, please refer to the "Getting Started" vignette for the
##' artemis package.
##'
##' There are two different modeling functions in the artemis package,
##' \code{eDNA_lm} and \code{eDNA_lmer}. \code{eDNA_lm} is for fitting
##' a fixed effects model, while \code{eDNA_lmer} is for fitting a
##' mixed or random effects model. Both models are fit using the
##' \code{rstan::sampling} function, which uses a Hamiltonian Monte
##' Carlo algorithm to estimate parameters for the model. Users are
##' encouraged to refer to the documentation for Stan and RStan at
##' \url{https://mc-stan.org/users/documentation/} for details about
##' how models are fit.
##'
##' @section Diagnosing warning and error messages:
##' 
##' The models have been written in Stan with key focus on robustness
##' and speed. However, it is possible that users might encounter
##' issues. Typically, these issues will be highlighted by warning
##' messages, which will be displayed regardless of the value of
##' \code{verbose}. Often times, these warnings can be resolved by
##' increasing the number of iterations that the HMC algorithm runs by
##' specifying \code{iters} to be a larger value. This should be the
##' first action attempted, as increasing the \code{iters} increases
##' both the warm-up and sampling iterations. If users continue to
##' have issues, additional control arguments can be passed to
##' \code{rstan::sampling} via the \code{...} argument.
##' 
##' @title Fit eDNA Model
##'
##' @aliases eDNA_lm
##' 
##' @param formula a formula, specifying the relationship between the
##'     predictors and the latent variable eDNA concentration.
##' @param data data.frame, with the response and predictors
##' @param std_curve_alpha the alpha value for the formula for
##'     converting between log(eDNA concentration) and CQ value
##' @param std_curve_beta the beta value for the formula for
##'     converting between log(eDNA concentration) and CQ value
##' @param upper_Cq numeric, the upper limit on CQ detection. Any
##'     value of log(concentration) which would result in a value
##'     greater than this limit is instead recorded as the limit.
##' @param n_chain integer, the number of chains to use. Please note
##'     that less than two is not recommended.
##' @param iters integer, the number of iterations per chain
##' @param verbose logical, when TRUE output from
##'     \code{rstan::sampling} is written to the console.
##' @param sink_file character, a file to write the console output to
##'     if \code{verbose = FALSE}, by default writes to
##'     \code{tempfile()}
##' @param ... additional arguments passed to
##'     \code{\link[rstan]{sampling}}
##' @return S4 object, with the following slots: \describe{
##'     \item{ln_conc}{matrix, the posterior samples for the latent
##'     variable, eDNA concentration} \item{Cq_star}{matrix, the
##'     posterior prediction for the observed response}
##'     \item{betas}{array, the posterior estimates for the betas for
##'     the linear model} \item{sigma_Cq}{array, the posterior
##'     estimates for the measurement error of CQ}
##'     \item{formula}{formula, the original formula used in the
##'     model} \item{x}{data.frame, the model matrix used in the
##'     model} \item{std_curve_alpha}{numeric, the std. curve value
##'     used} \item{std_curve_beta}{numeric, the std. curve value
##'     used} \item{upper_Cq}{numeric, the upper limit for observed CQ
##'     used} \item{stanfit}{stanfit, the original results from
##'     \code{rstan::sampling}} }
##' @author Matt Espe
##'
##' @examples
##'
##' ## Fixed effect model
##' ans = eDNA_lm(Cq ~ Distance, eDNA_data,
##'               std_curve_alpha = 21.2, std_curve_beta = -1.5)
##'
##' ## Mixed-effect model
##' ## This takes a while to run
##' \dontrun{
##' ans2 = eDNA_lmer(Cq ~ Distance + (1|SampleID), eDNA_data,
##'                  std_curve_alpha = 21.2, std_curve_beta = -1.5)
##'
##' }
##'
##' @export
eDNA_lmer = function(formula, data, 
                     std_curve_alpha, std_curve_beta,
                     upper_Cq = 40,
                     n_chain = 4L, iters = 500L,
                     verbose = FALSE,
                     sink_file = tempfile(), ...)
{
    
    ml = gen_model_list_lmer(formula, data)
    
    md = prep_data(ml, std_curve_alpha, std_curve_beta,
                   Cq_upper = upper_Cq, type = "model")
    md$y = ml$y

    fit = run_model(model = stanmodels$eDNA_lmer,
                    data = md, n_chain = n_chain, iters = iters,
                    verbose = verbose, sink_file = sink_file, ...)
    fit = load_slots_model(fit)
    
    fit = as(fit, "eDNA_model_lmer")
    fit@random_x = as.data.frame(md$rand_x)
    fit@random_sd = extract(fit@stanfit, "rand_sigma")$rand_sigma
    
    return(fit)
}

run_model = function(model = if(length(data$prior_mu)) stanmodels$eDNA_lm_prior else stanmodels$eDNA_lm,
                     data, n_chain,
                  iters, verbose, sink_file, ...)

{
    if(!verbose) sink(sink_file)
    
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
