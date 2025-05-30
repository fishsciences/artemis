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
    if(nrow(mf) != nrow(data))
        stop("NA values present in response or predictor cols. Please remove NA values prior to modeling.")
    
    y = model.response(mf)
    x = model.matrix(attr(mf, "terms"), data)

    return(list(x = x, y = y))    
}

gen_model_list_lmer = function(formula, data)
{
    mf = lFormula(formula, data)

    x = mf$X
    if(nrow(x) != nrow(data))
        stop("NA values present in response or predictor cols. Please remove NA values prior to modeling.")
    
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

gen_model_list_lm_zip = function(formula, data)
{
  # modeled from pscl package
  # Unpacking the mixed formula into two components
  if (!identical(formula[[3]][[1]], as.name("|")))
    stop("Please specify a zero-inflated model using `|` to separate model components.")
  
  ffc = . ~ .
  ffz = ~ .
  ffc[[2]] = formula[[2]]
  ffc[[3]] = formula[[3]][[2]]
  ffz[[3]] = formula[[3]][[3]]
  ffz[[2]] = formula[[2]]
  
  mf = model.frame(ffc, data)
  mfz = model.frame(ffz, data)
  
  if(nrow(mf) != nrow(data) || nrow(mfz) != nrow(data))
    stop("NA values present in response or predictor cols. Please remove NA values prior to modeling.")

  y = model.response(mf)
  x = model.matrix(attr(mf, "terms"), data)
  xz = model.matrix(attr(mfz, "terms"), data)
  
  return(list(x = x, xz = xz, y = y))    
}

gen_model_list_lmer_zip = function(formula, data)
{
  ffr = . ~ .
  ffz = ~ .
  ffr[[2]] = formula[[2]]
  
  ffr[[3]] = formula[[3]][[2]]
  ffz[[2]] = formula[[3]][[3]]

  mfz = model.frame(ffz, data)
  
  ## use other gen functions to parse random effects
  ans = gen_model_list_lmer(ffr, data)
  ans$xz = model.matrix(attr(mfz, "terms"), data)

  ans
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

