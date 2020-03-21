ranef.eDNA_model = function(object, FUN = quantile,
                            probs = c(0.025, 0.5, 0.975), ...)
{
    rands = extract(object@stanfit, pars = "rand_betas")$rand_betas
    ans = apply(rands, 2, FUN, probs)
    colnames(ans) = colnames(object@random_x)
    t(ans)
}
