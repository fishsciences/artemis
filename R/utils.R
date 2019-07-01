len = function(x)
{
    length(unique(x))
}

get_shared_rand = function(X)
{
    idx = lapply(X, unique)
    rep(seq(length(idx)), times = sapply(idx, length))
}

relevel_rands = function(x)
{
    m = c(0, cumsum(sapply(x, max)))
    as.data.frame(sapply(seq(ncol(x)), function(i){
        x[[i]] + m[i]
    }))
}

excludes_zero = function(ci)
{
    ci[2] < 0 | ci[1] > 0
}

in_interval = function(ci, beta)
{
    beta > ci[1] & beta < ci[2]
}

power_calc = function(interval_excludes_mat)
{
    rowSums(interval_excludes_mat) / ncol(interval_excludes_mat)
}

Cq2conc = function(Cq, alpha = 21.167769, beta = -1.528683)
    # Cq = -1.528683 * log(conc) + 21.167769
{
   exp((Cq - alpha) / beta)
}
