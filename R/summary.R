# S3 methods
##' Summary of eDNA simulations
##'
##' Summary of eDNA simulations
##' 
##' @title Summary of eDNA simulations
##'
##' @param object an object of class eDNA_simulation_*
##' @param var the simulated variable to summarize, either "Cq_star"
##'     or "ln_conc"
##' @param probs probabilities for the summary of the posterior
##'     samples
##' @param ... currently ignored
##' @return data.frame
##' @author Matt Espe
##'
##' @method summary eDNA_simulation
##' @export
summary.eDNA_simulation = function(object, var = "Cq_star",
                                   probs = c(0.025, 0.5, 0.975), ...)
{
    if(nrow(slot(object, var)) > 1){
        warning("Multiple sims present - computing summary on rep means")
        y = colMeans(slot(object, var))
    } else {
        y = slot(object, var)
    }

    ## Super ugly - clean up later
    qtl = get_marginals(y, object@x, quantile, probs = probs)
    qtl = lapply(seq_along(qtl), function(i) {
        tmp = qtl[i]
        tmp = as.data.frame(tmp, stringsAsFactors = FALSE)
        colnames(tmp) = paste0(probs * 100, "%")
        tmp = cbind(variable = names(qtl)[i],
                    level = as.numeric(rownames(qtl[[i]])),
                    tmp, stringsAsFactors = FALSE)
        tmp
    })
    qtl = do.call(rbind, qtl)
    m = get_marginals(y, object@x, mean)
    
    dt = if(var == "Cq_star"){
             get_marginals(y, object@x, p_detect,
                           thresh = object@upper_Cq)
         } else NA
    
    ans = cbind(qtl, mean = unlist(m), p_detect = unlist(dt))
    rownames(ans) = NULL
    structure(ans,
              class = c("eDNA_simulation.summary", "data.frame"))
    
}

##' Summary of eDNA model
##'
##' Summary of eDNA model
##' 
##' @title Summary of an eDNA model
##' @param object an object of class eDNA_model_*
##' @param probs probabilities for the quantiles of the posterior
##'     samples
##' @param ... currently ignored
##' @return data.frame, with summary of the fixed effects from the
##'     model
##'
##' @author Matt Espe
##' @method summary eDNA_model
##' @export
summary.eDNA_model = function(object, probs = c(0.025, 0.5, 0.975), ...)
{
    res = c(mean = mean(object@sigma_ln_eDNA), quantile(object@sigma_ln_eDNA, probs))
    nms = "ln(eDNA)_sigma"
    
    if(length(object@betas)){
        res = cbind(summarize_par(object@betas, probs), res)
        nms = c(colnames(object@x), nms)
    }
    
    if(length(object@intercept)){
        res = cbind(c(mean = mean(object@intercept),
                    quantile(object@intercept, probs)), res)
        nms = c("(Intercept)", nms)
    }

    res = data.frame(t(res))
    rownames(res) = nms
    colnames(res) = c("Mean", paste0(probs * 100, "%"))
  
    structure(res,
              iter = object@stanfit@stan_args$iter,
              class = c("eDNA_model.summary", "data.frame"))
}

summarize_par = function(x, p)
{
    means = apply(x, 2, mean, simplify = TRUE)
    qs = apply(x, 2, quantile, prob = p, simplify = TRUE)
    rbind(mean = means, qs)
}

summary.eDNA_predict_lm = function(object, probs = c(0.025, 0.5, 0.975), ...)
{
    # Placeholder
    object

}

##' Summary method for eDNA p(detect)
##'
##' Summary method for eDNA p(detect)
##' @title Summary method for eDNA p(detect)
##' @param object an object of class eDNA_p_detect
##' @param probs probabilities for summary, passed to \code{quantile}
##' @param ... ignored
##' @return a data.frame, with summary statistics of the object
##' @author Matt Espe
##' @method summary eDNA_p_detect
##' @export
summary.eDNA_p_detect = function(object, probs = c(0.025, 0.5, 0.975), ...)
{
    if(!is.matrix(object))
        return(data.frame(n_rep = attr(object, "reps"),
                          p_detect = as.numeric(object)))

    mn = apply(object, 2, mean)
    
    ci = apply(object, 2, quantile, probs)

    as.data.frame(cbind(n_rep = attr(object, "reps"), mean = mn, t(ci)))
}


summary.eDNA_predict_lm = function(object, probs = c(0.025, 0.5, 0.975),
                                   FUN = quantile)
{
    ll = if(!attr(object, "interval")){
             lapply(object, function(x) {
                 if(is.null(x)) return(NULL)
                 t(apply(x, 1, FUN, probs))
             })
         } else {
             object
         }
    
    
}

