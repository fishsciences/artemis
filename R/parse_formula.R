## Functions for parsing a user's formula
## and outputing data in the format that the Stan models
## expect


is_lme4 = function(formula)
{
    any(grepl("\\(.*\\|.*\\)", as.character(formula)))
}

gen_model_list_lm = function(formula, data)
{
    mf = model.frame(formula, data)
    y = model.response(mf)
    x = model.matrix(attr(mf, "terms"), data)

    return(list(x = x, y = y))    
}

gen_model_list_lmer = function(formula, data)
{
    mf = lFormula(formula, data)
    check_formula(mf)

    x = mf$X
    
    y = mf$fr[, as.character(mf$formula[2L])]
    if (is.matrix(y) && ncol(y) == 1L) 
        y = as.vector(y)

    rand_x = t(as.matrix(mf$reTrms$Zt))
    groups = mf$reTrms$Lind
    
    n_rand = ncol(rand_x)
    n_grp = length(unique(groups))

    return(list(y = y, x = x,
                rand_x = rand_x,
                groups = groups,
                n_rand = n_rand,
                n_grp = n_grp))
}

check_formula = function(formula)
{
    if(!all(unlist(formula$reTrms$cnms) == "(Intercept)"))
        stop("Sorry, only random intercepts are supported at this time.")
    return(NULL)
}
