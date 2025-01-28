
##' @rdname eDNA_lmer
##' @export
eDNA_lm = function(formula, data, 
                   std_curve_alpha, std_curve_beta,
                   upper_Cq = 40, 
                   prior_intercept = normal(location = -15, scale = 10),
                   priors = normal(),
                   cache_dir = getOption("artemis_cache_dir", R_user_dir("artemis", "cache")),
                   ...)
{
    eDNA_lm_shared(model_type = "lm",
                   formula, data, 
                   std_curve_alpha, std_curve_beta,
                   upper_Cq,
                   prior_intercept,
                   priors,
                   cache_dir,
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
##' There are four different modeling functions in the artemis
##' package, \code{eDNA_lm}, \code{eDNA_lmer}, \code{eDNA_zinf_lm},
##' \code{eDNA_zinf_lmer}.  \code{eDNA_lm} is for fitting a fixed
##' effects model, while \code{eDNA_lmer} is for fitting a mixed or
##' random effects model. The *_zinf versions implement a
##' zero-inflated version of their respective lm function. All models
##' are fit using the \code{rstan::sampling} function, which uses a
##' Hamiltonian Monte Carlo algorithm to estimate parameters for the
##' model. Users are encouraged to refer to the documentation for Stan
##' and RStan at \url{https://mc-stan.org/users/documentation/} for
##' details about how models are fit.
##'
##' @section Diagnosing warning and error messages:
##' 
##' The models have been written in Stan with key focus on robustness
##' and speed. However, it is possible that users might encounter
##' issues. Typically, these issues will be highlighted by warning
##' messages coming from \code{rstan::sampling}. Often times, these
##' warnings can be resolved by increasing the number of iterations
##' that the HMC algorithm runs by specifying \code{iters} to be a
##' larger value. This should be the first action attempted, as
##' increasing the \code{iters} increases both the warm-up and
##' sampling iterations. If users continue to have issues, additional
##' control arguments can be passed to \code{rstan::sampling} via the
##' \code{...} argument.
##' 
##' @title Fit eDNA Model
##'
##' @aliases eDNA_lm
##' 
##' @param formula a formula, specifying the relationship between the
##'     predictors and the latent variable eDNA concentration.
##' @param data data.frame, with the response and predictors
##' @param std_curve_alpha the alpha (intercept) value for the formula
##'     for converting between log(eDNA concentration) and CQ value
##' @param std_curve_beta the beta (slope) value for the formula for
##'     converting between log(eDNA concentration) and CQ value
##' @param upper_Cq numeric, the upper limit on CQ detection. Any
##'     value of log(concentration) which would result in a value
##'     greater than this limit is instead recorded as the limit.
##' @param prior_intercept named list such as created by
##'     \code{rstanarm::normal}. The list must contain elements named
##'     "location" and "scale", which are the location and scale for a
##'     normal prior over the intercept. Ignored when the intercept is
##'     omitted in the model formula.
##' @param prior_random_variance the prior on variance of the random
##'     effects. Defaults to exponential distribution with rate 1.
##' @param priors named list such as created by
##'     \code{rstanarm::normal}. The list must contain elements named
##'     "location" and "scale", which are the location and scale for a
##'     normal prior over the betas, and "autoscale". If a single
##'     value is provided, this value will be repeated for each
##'     beta. If \code{autoscale = TRUE}, the scale of the priors is
##'     scaled by the sd of the predictors similar to rstanarm handles
##'     them.
##' @param cache_dir a directory which stores the Stan model code and
##'     compiled models. This directory can be set up and the models
##'     compiled via \code{compile_models()}.  Can be set with the
##'     option "artemis_cache_dir", otherwise defaults to the
##'     directory returned by \code{R_user_dir()}.
##' @param ... additional arguments passed to
##'     \code{\link[cmdstanr]{sample}}
##' @return S4 object, with the following slots:
##' \describe{
##'   \item{ln_conc}{matrix, the posterior samples for the latent
##'     variable, eDNA concentration}
##' \item{Cq_star}{matrix, the
##'     posterior prediction for the observed response}
##' \item{betas}{array, the posterior estimates for the betas for
##'     the linear model}
##' \item{sigma_ln_eDNA}{array, the posterior
##'     estimates for the measurement error of ln_eDNA}
##' \item{formula}{formula, the original formula used in the
##'     model}
##' \item{x}{data.frame, the model matrix used in the
##'     model}
##' \item{std_curve_alpha}{numeric, the std. curve intercept value
##'     used}
##' \item{std_curve_beta}{numeric, the std. curve slope value
##'     used}
##' \item{upper_Cq}{numeric, the upper limit for observed CQ
##'     used}
##' \item{fit}{fit, the original results from
##'     \code{cmdstanr::sample}} }
##' @author Matt Espe
##'
##' @examples
##' \donttest{
##' ## Fixed effect model
##' ans = eDNA_lm(Cq ~ Distance_m, eDNA_data,
##'               std_curve_alpha = 21.2, std_curve_beta = -1.5)
##'
##' ## Mixed-effect model
##' ## This takes a while to run
##' ans2 = eDNA_lmer(Cq ~ Distance_m + (1|FilterID), eDNA_data,
##'                  std_curve_alpha = 21.2, std_curve_beta = -1.5)
##'
##' }
##'
##' @export
eDNA_lmer = function(formula, data, 
                     std_curve_alpha, std_curve_beta,
                     upper_Cq = 40, 
                     prior_intercept = normal(location = -15, scale = 10),
                     priors = normal(),
                     prior_random_variance = exponential(),
                     cache_dir = getOption("artemis_cache_dir", R_user_dir("artemis", "cache")),
                     ...)
{
    eDNA_lm_shared(model_type = "lmer",
                   formula, data, 
                   std_curve_alpha, std_curve_beta,
                   upper_Cq,
                   prior_intercept,
                   priors,
                   prior_random_variance,
                   cache_dir,
                   ...)
}

run_model = function(model,
                     data, ...)
    
{
    m = model$sample(data = data, ...)
    fit = as_cmdstan_fit(m$output_files())
    fit = as(fit, "eDNA_model")
    return(fit)
}

load_slots_model = function(obj, model_type)
  # This grabs data from 'md', the input data to the model and adds it
  # to the final fit returned object, thus preserving the inputs.
{
  obj@formula = get("formula", parent.frame())
  obj@x = as.data.frame(get("md", parent.frame())$X)
  obj@std_curve_alpha = get("std_curve_alpha", parent.frame())
  obj@std_curve_beta = get("std_curve_beta", parent.frame())
  obj@upper_Cq = get("upper_Cq", parent.frame())
  
  if(model_type == "eDNA_model_lmer"){
    #obj@random_x = as.data.frame(get("md", parent.frame())$X)
    #obj@random_sd = as.data.frame(get("md", parent.frame())$X)
  }
  if(model_type %in% c("zero_inf_lm", "zero_inf_lmer")){
    obj@xz = as.data.frame(get("md", parent.frame())$Xz)
  }
  
  obj
}

##' @rdname eDNA_lmer
##' @export
eDNA_zinf_lm = function(formula, data, 
                        std_curve_alpha, std_curve_beta,
                        upper_Cq = 40,
                        prior_intercept = normal(location = -15, scale = 10),
                        priors = normal(), 
                        cache_dir = tools::R_user_dir("artemis", "cache"),
                        ...)
{
    # from lm
  eDNA_lm_shared(model_type = "zero_inf_lm",
                 formula, data, 
                 std_curve_alpha, std_curve_beta,
                 upper_Cq,
                 prior_intercept,
                 priors,
                 cache_dir,
                   ...)
}
##' @rdname eDNA_lmer
##' @export
eDNA_zinf_lmer = function(formula, data, 
                          std_curve_alpha, std_curve_beta,
                          upper_Cq = 40,
                          prior_intercept = normal(location = -15, scale = 10),
                          priors = normal(),
                          prior_random_variance = exponential(),
                          cache_dir = tools::R_user_dir("artemis", "cache"),
                          ...)
{
  
  # from lm
  eDNA_lm_shared(model_type = "zero_inf_lmer",
                 formula, data, 
                 std_curve_alpha, std_curve_beta,
                 upper_Cq,
                 prior_intercept,
                 priors,
                 prior_random_variance,
                 cache_dir,
                 ...)
}

# houses most lm() code, which is similar between the lm, lmer, and
# zero-inflated
eDNA_lm_shared = function(model_type, 
                          formula, data, 
                          std_curve_alpha, std_curve_beta,
                          upper_Cq = 40,
                          prior_intercept = normal(location = -15, scale = 10),
                          priors = normal(), 
                          prior_random_variance = exponential(),
                          cache_dir = tools::R_user_dir("artemis", "cache"),
                          ...)
{
    md_pars = get_mod_funs(model_type, cache_dir)
    
    # from lm
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]

    mf[[1]] = md_pars$gen_fun

    ml = eval(mf, parent.frame(1L))
       
    # This works because Stan ignores extra input data
    md = md_pars$prep_fun(ml, std_curve_alpha, std_curve_beta,
                          Cq_upper = upper_Cq, 
                          prior_int = prior_intercept,
                          prior_b = priors, rand_sd = prior_random_variance)
    # md$y = ml$y
    fit = run_model(model = md_pars$mod,
                    data = md, ...)
    fit = load_slots_model(fit, model_type)
    fit = as(fit, md_pars$model_class)
    
    return(fit)
}

# Put this here to allow us to easily swap out pieces later
get_mod_funs = function(model_type, cache_dir)
{
    mn = switch(model_type,
                lm = "eDNA_lm.stan",
                lmer = "eDNA_lmer.stan",
                zero_inf_lm = "eDNA_lm_zinf.stan",
                zero_inf_lmer = "eDNA_lmer_zinf.stan",
                stop("Unknown model type"))
    
    mod_file = file.path(cache_dir, mn)

    # Not sure about this - will recompile the models
    # but gets around issue with checks and installs
    compile_models(model_names = mn,
                   cache_dir,
                   rewrite = FALSE,
                   verbose = FALSE)
    
    mod = cmdstan_model(mod_file)
    
    list(mn = mn,
         mod = mod,
         gen_fun = switch(model_type,
                          lm = quote(gen_model_list_lm),
                          lmer = quote(gen_model_list_lmer),
                          zero_inf_lm = quote(gen_model_list_lm_zip),
                          zero_inf_lmer = quote(gen_model_list_lmer_zip)),
         
         prep_fun = switch(model_type,
                           lm = prep_data.lm,
                           lmer = prep_data.lmer,
                           zero_inf_lm = prep_data.zip, 
                           zero_inf_lmer = prep_data.zipr),
         
         
         model_class = switch(model_type,
                              lm = "eDNA_model_lm",
                              lmer = "eDNA_model_lmer",
                              zero_inf_lm = "eDNA_model_zip", 
                              zero_inf_lmer = "eDNA_model_zipr")
         )
}
