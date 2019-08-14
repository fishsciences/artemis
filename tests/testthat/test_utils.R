
context("Utilities")
test_that("Random effect generation", {

    X = expand.grid(1:4, 1:2)
    ans = get_shared_rand(X)
    expect_equal(ncol(X), len(ans))

    X = expand.grid(1:4, 1:2, 1:10)
    ans = get_shared_rand(X)
    expect_equal(ncol(X), len(ans))

})


test_that("Re-level random effects", {
    X = expand.grid(1:4, 1:2, 1:10)
    ans = relevel_rands(X)
    ans2 = get_shared_rand(ans)
    expect_equal(max(unlist(ans)), length(ans2))
    
    ii = structure(c(0.198978687755089, -0.352891790970438, -5.70524972005077, 
                     0.200586796657588, -0.33659423124757, -4.94471656961896),
                   .Dim = 3:2,
                   .Dimnames = list(
                       c("betas[1]", "betas[2]", "betas[3]"), c("2.5%", "97.5%")))
    expect_true(all(apply(ii, 1, excludes_zero)))
})

test_that("Conf 2 prob", {
    expect_true(length(conf_to_probs(0.90)) == 2)
    expect_equal(conf_to_probs(0.9), c(0.05, 0.95))
})
