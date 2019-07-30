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
    structure(ans,
              class = c("eDNA_simulation.summary", "data.frame"))
    
}

summary.eDNA_model = function(object, probs = c(0.025, 0.5, 0.975), ...)
{
    res = t(apply(object@betas, 2, quantile, prob = probs))
    res = as.data.frame(res)
    rownames(res) = colnames(object@x)
    res$mean = colMeans(object@betas)
    
    structure(res,
              class = c("eDNA_model.summary", "data.frame"))
}
