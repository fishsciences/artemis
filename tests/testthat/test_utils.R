
context("Utilities")


test_that("Conf 2 prob", {
    expect_true(length(conf_to_probs(0.90)) == 2)
    expect_equal(conf_to_probs(0.9), c(0.05, 0.95))
})
