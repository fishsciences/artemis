## These functions are for generating simulated data to be used
## to calculate power for the eDNA experiments
## M. Espe
## March 2019

##' @rdname sim_eDNA_lmer
##' @export
sim_eDNA_lm = function(formula, variable_list,
                       betas, sigma_ln_eDNA,
                       std_curve_alpha, std_curve_beta,
                       n_sim = 1L,
                       upper_Cq = 40,
                       prob_zero = 0.08,
                       X = expand.grid(variable_list),
                       verbose = FALSE)
{
    if(!has_response(formula))
        stop("Please provide a dummy response variable for simulations")
    ml = gen_model_list_lm(formula, X)
    
    if(ncol(ml$x) != length(betas))
        stop("Please provide one beta per model term \n",
             "Provided: ", length(betas), "\n",
             "Required: ", ncol(ml$x), "\n")

    md = prep_data.sim(ml, std_curve_alpha, std_curve_beta, sigma_ln_eDNA, betas,
                       prob_zero, prior_int = normal(), prior_b = normal())

    sims = sampling(stanmodels$eDNA_sim_omni, data = md, chains = 1L,
                    algorithm = "Fixed_param", iter = n_sim, warmup = 0L,
                    refresh = ifelse(verbose, 100, -1), show_messages = verbose,
                    open_progress = FALSE)
    # hacky
    sims = as(sims, "eDNA_simulation_lm")
    sims = load_slots(sims)
    return(sims)
}

##' Simulate eDNA data
##'
##' These functions allow for computationally efficient simulation of
##' Cq values from a hypothetical eDNA sampling experiment via a
##' series of effect sizes (\code{betas}) on a number of predictor or
##' variable levels (\code{variable_levels}). The mechanism for this
##' model is described in detail in the artemis "Getting Started"
##' vignette.
##'
##' The simulation functions call to specialized functions which are
##' written in Stan and are compiled to provide speed. This also
##' allows the simulation functions and the modeling functions to
##' reflect the same process at the code level.
##'
##' @section Diagnosing "unrealistic" simulations:
##'
##' Users will find that sometimes the simulationed response (i.e. Cq
##' values) produced by this function are not similar to expected data
##' collected from a sampling experiment. This circumstance suggests
##' that there is a mismatch between the assumptions of the model and
##' the data generating process in the field. For these circumstances,
##' we suggest:
##'
##' \enumerate{
##' \item Check that the \code{betas} provided are the
##' effect sizes on the predictor on the log[eDNA concentration], and
##' not the Cq values.
##'
##' \item Check that the variable levels provided are representative
##' of real-world circumstances. For example, a sample volume of 0 ml
##' is not possible.
##'
##' \item Verify the values for the standard curve alpha and
##' beta. These are specific to each calibration for the lab, so it is
##' important that you use the same conversion between Cq values and
##' log[eDNA concentration] as the comparison data.
##'
##' }
##' 
##' @title Simulate eDNA data
##'
##' @aliases sim_eDNA_lm
##' 
##' @param formula a model formula, e.g. \code{y ~ x1 + x2}. For
##'     \code{sim_eDNA_lmer}, random intercepts can also be provided,
##'     e.g. \code{ ( 1 | rep ) }.
##' @param variable_list a named list, with the levels that each
##'     variable can take.  Please note that the variables listed in
##'     the formula, including the response variable, must be present
##'     in the variable_list or in the X design matrix.  Extra
##'     variables, i.e. variables which do not occur in the formula,
##'     are ignored.
##' @param betas numeric vector, the beta for each variable in the
##'     design matrix
##' @param sigma_ln_eDNA numeric, the measurement error on ln[eDNA].
##' @param sigma_rand numeric vector, the stdev for the random
##'     effects. There must be one sigma per random effect specified
##' @param std_curve_alpha the alpha value for the formula for
##'     converting between log(eDNA concentration) and CQ value
##' @param std_curve_beta the beta value for the formula for
##'     converting between log(eDNA concentration) and CQ value
##' @param n_sim integer, the number of cases to simulate
##' @param upper_Cq numeric, the upper limit on CQ detection. Any
##'     value of log(concentration) which would result in a value
##'     greater than this limit is instead recorded as the limit.
##' @param prob_zero numeric, between 0 and 1. The probability of
##'     seeing a non-detection (i.e., a "zero") via the zero-inflated
##'     mechanism. Defaults to 0.08.
##' @param X optional, a design matrix. By default, this is created
##'     from the variable_list using \code{expand.grid()}, which
##'     creates a balanced design matrix. However, the user can
##'     provide their own \code{X} as well, in which case the
##'     variable_list is ignored. This allows users to provide an
##'     unbalanced design matrix.
##' @param verbose logical, when TRUE output from
##'     \code{rstan::sampling} is written to the console.
##' @return S4 object of class "eDNA_simulation_{lm/lmer}" with the
##'     following slots:
##' \describe{
##' \item{ln_conc matrix}{the
##'     simulated log(concentration)}
##' \item{Cq_star matrix}{the
##'     simulated CQ values, including the measurement error}
##' \item{formula}{the formula for the simulation}
##' \item{variable_levels}{named list, the variable levels used
##'     for the simulation}
##' \item{betas}{numeric vector, the betas for
##'     the simulation}
##' \item{x}{data.frame, the design matrix}
##' \item{std_curve_alpha numeric}{the alpha for the std curve
##'     conversion}
##' \item{std_curve_beta numeric}{the alpha for the
##'     std curve conversion}
##' \item{upper_Cq}{the upper limit for CQ}
##'     }
##' 
##' @author Matt Espe
##' @examples
##' 
##' ## Includes extra variables
##' vars = list(Intercept = -10.6,
##'             distance = c(0, 15, 50),
##'             volume = c(25, 50),
##'             biomass = 100,
##'             alive = 1,
##'             tech_rep = 1:10,
##'             rep = 1:3, Cq = 1)
##'
##' ## Intercept only
##' ans = sim_eDNA_lm(Cq ~ 1, vars,
##'                       betas = c(intercept = -15),
##'                       sigma_ln_eDNA = 1e-5,
##'                       std_curve_alpha = 21.2, std_curve_beta = -1.5)
##' 
##' print(ans)
##'
##' ans = sim_eDNA_lm(Cq ~ distance + volume, vars,
##'                   betas = c(intercept = -10.6, distance = -0.05, volume = 0.1),
##'                   sigma_ln_eDNA = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5)
##' @export
sim_eDNA_lmer = function(formula, variable_list,
                         betas, sigma_ln_eDNA,
                         sigma_rand,
                         std_curve_alpha, std_curve_beta,
                         n_sim = 1L,
                         upper_Cq = 40,
                         prob_zero = 0.08,
                         X = expand.grid(variable_list),
                         verbose = FALSE)
{
    if(!has_response(formula))
        stop("Please provide a dummy response variable for simulations")
    
    if(is_lme4(formula) && length(sigma_rand) == 0)
        stop("You must provide the sd of the random effects")

    ml = gen_model_list_lmer(formula, X)

    if(ncol(ml$x) != length(betas))
        stop("Please provide one beta per model term \n",
             "Provided: ", length(betas), "\n",
             "Required: ", ncol(ml$x), "\n")

    if(length(sigma_rand) != ml$n_grp)
        stop("You must provide one sd for each random effect\n",
             "Provided: ", length(sigma_rand), "\n",
             "Required: ", ncol(ml$groups), "\n",
             "Random effects: ", colnames(ml$groups))

    
    md = prep_data.sim(ml, std_curve_alpha, std_curve_beta, sigma_ln_eDNA,
                   betas = betas, prob_zero, rand_sd = sigma_rand,
                   prior_int = normal(), prior_b = normal())

    sims = sampling(stanmodels$eDNA_sim_omni, data = md, chains = 1L,
                    algorithm = "Fixed_param", iter = n_sim, warmup = 0L,
                    refresh = ifelse(verbose, 100, -1), show_messages = verbose,
                    open_progress = FALSE)        
    # hacky
    sims = as(sims, "eDNA_simulation_lmer")
    sims = load_slots(sims)
    
    return(sims)
}


has_response = function(formula)
{
    attr(terms(formula), "response") > 0
}

load_slots = function(obj)
    ## Potentially really dangerous - think about this

    ## Want to avoid copying this code 2x for each sim, but
    ## grabbing from the parent frame might bite us
{
    obj@formula = get("formula", parent.frame())
    obj@betas = get("betas", parent.frame())
    obj@variable_levels = get("variable_list", parent.frame())
    obj@x = as.data.frame(get("md", parent.frame())$X)
    obj@std_curve_alpha = get("std_curve_alpha", parent.frame())
    obj@std_curve_beta = get("std_curve_beta", parent.frame())
    obj@upper_Cq = get("upper_Cq", parent.frame())
    
    obj
}

