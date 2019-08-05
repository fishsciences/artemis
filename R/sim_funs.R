## These functions are for generating simulated data to be used
## to calculate power for the eDNA experiments
## M. Espe
## March 2019

##' @rdname sim_eDNA_lmer 
sim_eDNA_lm = function(formula, variable_list,
                       betas, sigma_Cq,
                       std_curve_alpha, std_curve_beta,
                       n_sim = 1L,
                       upper_Cq = 40,
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

    md = prep_data(ml, std_curve_alpha, std_curve_beta, sigma_Cq, betas, type = "sim")

    if(!verbose) sink(tempfile())

    sims = sampling(stanmodels$eDNA_sim, data = md, chains = 1L,
                    algorithm = "Fixed_param", iter = n_sim,
                    refresh = ifelse(verbose, 100, -1), show_messages = verbose)

    if(!verbose) sink()
    
    # hacky
    sims = as(sims, "eDNA_simulation_lm")
    sims = load_slots(sims)
    return(sims)
}

##' Simulate eDNA data
##'
##' Simulate eDNA data
##' 
##' @title Simulate eDNA data
##'
##' @aliases sim_eDNA_lm
##' 
##' @param formula a model formula, e.g. \code{y ~ x1 + x2}. For \code{sim_eDNA_lmer},
##'     random intercepts can also be provided, e.g. \code{ ( 1 | rep ) }. 
##' @param variable_list a named list, with the levels that each variable can take.
##'     Please note that the variables listed in the formula, including the
##'     response variable, must be present in the variable_list or in the X design matrix.
##'     Extra variables, i.e. variables which do not occur in the formula, are ignored.
##' @param betas numeric vector, the beta for each variable in the design matrix
##' @param sigma_Cq numeric, the measurement error on CQ.
##' @param sigma_rand numeric vector, the stdev for the random effects. There must be
##'     one sigma per random effect specified
##' @param std_curve_alpha the alpha value for the formula for converting between
##'     log(eDNA concentration) and CQ value
##' @param std_curve_beta the beta value for the formula for converting between
##'     log(eDNA concentration) and CQ value
##' @param n_sim integer, the number of cases to simulate
##' @param upper_Cq numeric, the upper limit on CQ detection. Any value of
##'     log(concentration) which would result in a value greater than this limit is
##'     instead recorded as the limit.
##' @param X optional, a design matrix. By default, this is created
##'     from the variable_list using \code{expand.grid()}, which creates
##'     a balanced design matrix. However, the user can provide their own \code{X}
##'     as well, in which case the variable_list is ignored. This allows users to
##'     provide an unbalanced design matrix. 
##' @param verbose logical, when TRUE output from \code{rstan::sampling} is written
##'     to the console. 
##' @return S4 object of class "eDNA_simulation_{lm/lmer}" with the following slots:
##' \item {ln_conc matrix} { the simulated log(concentration)}
##' \item {Cq_star matrix} {the simulated CQ values, including the measurement error}
##' \item {formula} { the formula for the simulation}
##' \item {variable_levels} { named list, the variable levels used for the simulation}
##' \item {betas}{ numeric vector, the betas for the simulation}
##' \item {x}{ data.frame, the design matrix}
##' \item {std_curve_alpha numeric}{ the alpha for the std curve conversion}
##' \item {std_curve_beta numeric}{ the alpha for the std curve conversion}
##' \item {upper_Cq}{ the upper limit for CQ}
##' 
##' 
##' @author Matt Espe
##' @examples
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
##'                       sigma_Cq = 1e-5,
##'                       std_curve_alpha = 21.2, std_curve_beta = -1.5)
##' 
##' print(ans)
##'
##' ans = sim_eDNA_lm(Cq ~ distance + volume, vars,
##'                   betas = c(intercept = -10.6, distance = -0.05, volume = 0.1),
##'                   sigma_Cq = 1, std_curve_alpha = 21.2, std_curve_beta = -1.5)
##' 

sim_eDNA_lmer = function(formula, variable_list,
                         betas, sigma_Cq,
                         sigma_rand,
                         std_curve_alpha, std_curve_beta,
                         n_sim = 1L,
                         upper_Cq = 40,
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

    if(length(sigma_rand) != ncol(ml$groups))
        stop("You must provide one sd for each random effect\n",
             "Provided: ", length(sigma_rand), "\n",
             "Required: ", ncol(ml$groups), "\n",
             "Random effects: ", colnames(ml$groups))

    
    md = prep_data(ml, std_curve_alpha, std_curve_beta, sigma_Cq,
                  betas = betas, rand_sd = sigma_rand,
                  type = "sim")

    if(!verbose) sink(tempfile())
    
    sims = sampling(stanmodels$eDNA_sim, data = md, chains = 1L,
                    algorithm = "Fixed_param", iter = n_sim, warmup = 0L,
                    refresh = ifelse(verbose, 100, -1), show_messages = verbose)
    
    if(!verbose) sink()

    # hacky
    sims = as(sims, "eDNA_simulation_lmer")
    sims = load_slots(sims)
    
    return(sims)
}


gen_rand_eff = function(X, sigma_rand_eff)
{
    sapply(seq_along(X), function(i){
        n_levels = length(unique(X[[i]]))
        rnorm(n_levels, mean = 0, sd = sigma_rand_eff)[X[[i]]]
    })
}


gen_eDNA_conc = function(X, betas, rands)
    ## Generate the concentration from a
    ## model matrix X and a vector betas
{
    y_hat = as.matrix(X) %*% betas + rowSums(rands)
    return(y_hat)
}

gen_Cq = function(ln_conc, sigma,
                  std_alpha, std_beta,
                  thresh)
    ## Generate Cq from a conc, and then truncate to thresh
{
    Cq_hat = ln_std_curve(exp(ln_conc), std_alpha, std_beta)
    Cq_star = rnorm(length(ln_conc), Cq_hat, sigma)
    Cq_star[Cq_star > thresh] = thresh
    return(Cq_star)
}

ln_std_curve = function(x, alpha = 21.167769, beta = -1.52868305)
    ## equation, alpha, and beta from 20180927_VS_standard_curve.xlsx
{
    beta * log(x) + alpha
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

