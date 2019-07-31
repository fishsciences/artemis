# S3 methods
summary.eDNA_simulation = function(object, var = "Cq_star",
                                   probs = c(0.025, 0.5, 0.975), ...)
{
    ## Super ugly - clean up later
    qtl = get_marginals(slot(object, var), object@x, quantile, probs = probs)
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
    m = get_marginals(slot(object, var), object@x, mean)
    
    dt = if(var == "Cq_star"){
             get_marginals(slot(object, var), object@x, p_detect,
                           thresh = object@upper_Cq)
         } else NA
    
    ans = cbind(qtl, mean = unlist(m), p_detect = unlist(dt))
    rownames(ans) = NULL
    structure(ans,
              class = c("eDNA_simulation.summary", "data.frame"))
    
}


summary.eDNA_model = function(object, probs = c(0.025, 0.5, 0.975), ...)
{
    res = apply(object@betas, 2, quantile, prob = probs, simplify = FALSE)
    
    if(!is.null(ncol(res))) {
        res = t(res)
        res = rbind(res, quantile(object@sigma_Cq, probs))
    } else {
        res = c(res, quantile(object@sigma_Cq, probs))
    }
    res = data.frame(res)

    colnames(res) = paste0(probs * 100, "%")
    
    rownames(res) = c(colnames(object@x), "CQ sd")
    res$mean = colMeans(cbind(object@betas, object@sigma_Cq))
    
    structure(res,
              iter = object@stanfit@stan_args$iter,
              class = c("eDNA_model.summary", "data.frame"))
}
