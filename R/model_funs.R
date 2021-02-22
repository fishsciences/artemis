
##' @rdname eDNA_lmer
##' @export
eDNA_lm = function(formula, data, 
                   std_curve_alpha, std_curve_beta,
                   upper_Cq = 40, 
                   prior_intercept = normal(location = -15, scale = 10),
                   priors = normal(), Cq_error_type = "fixed",
                   ...)
{
    eDNA_lm_shared(model_type = "lm",
                          formula, data, 
                          std_curve_alpha, std_curve_beta,
                          upper_Cq,
                          probability_zero = 0,
                          prior_intercept,
                          priors, Cq_error_type,
                   ...)
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
##' @param QR logical, use the QR decomposition. Typically more
##'     efficient when informative priors are not used. When TRUE, the
##'     priors are ignored and a standard normal prior is used on
##'     theta, the parameters estimates on the QR decomposed model
##'     matrix.
##' @param probability_zero numeric, between 0 and 1. The probability
##'     of a non-detection from a source other than low concentration
##'     of eDNA, e.g. a filter failure. Defaults to 8% (0.08), which
##'     was the estimated p(zero) from a daily sampling experiment.
##' @param prior_intercept named list such as created by
##'     \code{rstanarm::normal}. The list must contain elements named
##'     "location" and "scale", which are the location and scale for a
##'     normal prior over the intercept. Ignored when the intercept is
##'     omitted in the model formula.
##' @param priors named list such as created by
##'     \code{rstanarm::normal}. The list must contain elements named
##'     "location" and "scale", which are the location and scale for a
##'     normal prior over the betas, and "autoscale". If a single
##'     value is provided, this value will be repeated for each
##'     beta. If \code{autoscale = TRUE}, the scale of the priors is
##'     scaled by the sd of the predictors similar to rstanarm handles
##'     them.
##' @param Cq_error_type either "fixed" or "varying", specifying if
##'     measurement error is assumed to be the same for all values of
##'     CQ, or if it increases or decreases as CQ increases.
##' @param ... additional arguments passed to
##'     \code{\link[rstan]{sampling}}
##' @return S4 object, with the following slots:
##' \describe{
##'   \item{ln_conc}{matrix, the posterior samples for the latent
##'     variable, eDNA concentration}
##'   \item{Cq_star}{matrix, the
##'     posterior prediction for the observed response}
##'   \item{betas}{array, the posterior estimates for the betas for
##'     the linear model}
##'   \item{sigma_Cq}{array, the posterior
##'     estimates for the measurement error of CQ}
##'   \item{formula}{formula, the original formula used in the
##'     model}
##'   \item{x}{data.frame, the model matrix used in the
##'     model}
##'   \item{std_curve_alpha}{numeric, the std. curve value
##'     used}
##'   \item{std_curve_beta}{numeric, the std. curve value
##'     used}
##'   \item{upper_Cq}{numeric, the upper limit for observed CQ
##'     used}
##'   \item{stanfit}{stanfit, the original results from
##'     \code{rstan::sampling}} }
##' @author Matt Espe
##'
##' @examples
##' \dontrun{
##' 
##' ## Fixed effect model
##' ans = eDNA_lm(Cq ~ Distance_m, eDNA_data,
##'               std_curve_alpha = 21.2, std_curve_beta = -1.5)
##'
##' ## Mixed-effect model
##' ## This takes a while to run
##' ans2 = eDNA_lmer(Cq ~ Distance_m + (1|SampleID), eDNA_data,
##'                  std_curve_alpha = 21.2, std_curve_beta = -1.5)
##'
##' }
##'
##' @export
eDNA_lmer = function(formula, data, 
                     std_curve_alpha, std_curve_beta,
                     upper_Cq = 40, 
                     prior_intercept = normal(location = -15, scale = 10),
                     priors = normal(), Cq_error_type = "fixed", 
                     ...)
{
    # from lm
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]

    mf[[1]] = quote(artemis:::gen_model_list_lmer)
    ml = eval(mf, parent.frame(1L))
    
    # ml = gen_model_list_lmer(formula, data)
    
    md = prep_data.lmer(ml, std_curve_alpha, std_curve_beta,
                   Cq_upper = upper_Cq, 
                   prior_int = prior_intercept,
                   prior_b = priors, error_type = Cq_error_type)

    # md$y = ml$y

    fit = run_model(model_name = "eDNA_lmer.stan", data = md, ...)
    fit = load_slots_model(fit)
    
    fit = as(fit, "eDNA_model_lmer")
    fit@random_x = as.data.frame(md$rand_x)
    fit@random_sd = extract(fit@stanfit, "rand_sigma")$rand_sigma
    
    return(fit)
}

run_model = function(model_name = "eDNA_lm.stan",
                     model = cmdstan_model(model_file),
                     model_file = system.file("stan_files",
                                              model_name,
                                              package = "artemis"),
                     data, ...)
    
{
    
    m = model$sample(data = data, ...)
    fit = read_stan_csv(m$output_files())
    
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

##' @rdname eDNA_lmer
##' @export
eDNA_zinf_lm = function(formula, data, 
                        std_curve_alpha, std_curve_beta,
                        upper_Cq = 40,
                        probability_zero = 0.08,
                        prior_intercept = normal(location = -15, scale = 10),
                        priors = normal(), Cq_error_type = "fixed",
                        model_type,
                        ...)
{
    # from lm
    eDNA_lm_shared(model_type = "zero_inf_lm",
                          formula, data, 
                          std_curve_alpha, std_curve_beta,
                          upper_Cq,
                          probability_zero,
                          prior_intercept,
                          priors, Cq_error_type,
                          ...)
}
##' @rdname eDNA_lmer
##' @export
eDNA_zinf_lmer = function(formula, data, 
                          std_curve_alpha, std_curve_beta,
                          upper_Cq = 40,
                          probability_zero = 0.08,
                          prior_intercept = normal(location = -15, scale = 10),
                          priors = normal(), Cq_error_type = "fixed",
                          model_type,
                          ...)
{
    # from lm
    eDNA_lm_shared(model_type = "zero_inf_lmer",
                   formula, data, 
                   std_curve_alpha, std_curve_beta,
                   upper_Cq,
                   probability_zero,
                   prior_intercept,
                   priors, Cq_error_type,
                   ...)
}

# houses most lm() code, which is similar between the lm and zero-inflated
eDNA_lm_shared = function(model_type,
                          formula, data, 
                          std_curve_alpha, std_curve_beta,
                          upper_Cq = 40,
                          probability_zero = 0.08,
                          prior_intercept = normal(location = -15, scale = 10),
                          priors = normal(), Cq_error_type = "fixed",
                          ...)
{
    mn = switch(model_type,
                lm = "eDNA_lm.stan",
                lmer = "eDNA_lmer.stan",
                zero_inf_lm = "eDNA_lm_zinf.stan",
                zero_inf_lmer = "eDNA_lmer_zinf.stan",
                stop("Unknown model type"))
    # from lm
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]

    mf[[1]] = quote(artemis:::gen_model_list_lm)

    ml = eval(mf, parent.frame(1L))
       
    # This works because Stan ignores extra input data
    md = prep_data.lm(ml, std_curve_alpha, std_curve_beta,
                   Cq_upper = upper_Cq, prob_zero = probability_zero,
                   prior_int = prior_intercept,
                   prior_b = priors, error_type = Cq_error_type)
    # md$y = ml$y
    fit = run_model(model_name = mn,
                    data = md, ...)
    fit = load_slots_model(fit)
    fit = as(fit, "eDNA_model_lm")
    
    return(fit)
}
