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
    rt = mk_rand_mat(mf$reTrms)
    rand_x = rt$random_effects
    groups = rt$groups
    
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
    ## if(!all(unlist(formula$reTrms$cnms) == "(Intercept)"))
        ## stop("Sorry, only random intercepts are supported at this time.")
    return(NULL)
}

mk_rand_mat = function(rt)
{
    mats = lapply(rt$Ztlist, as.matrix)
    vars = mapply(function(x, var, nm)
        paste(rep(x, times = nlevels(var)), nm, sep = ":"),
        x = rt$cnms, var = rt$flist, nm = names(rt$flist),
        SIMPLIFY = FALSE)

    groups = as.integer(factor(unlist(vars), levels = unique(unlist(vars))))
    var_names = unlist(mapply(function(x, y) paste(x, rownames(y), sep = ""),
                              x = vars, y = mats,
                              SIMPLIFY = FALSE))
    rand_mat = do.call(cbind, lapply(mats, t))
    colnames(rand_mat) = var_names

    return(list(random_effects = rand_mat,
                groups = groups))
}

