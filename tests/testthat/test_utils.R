
context("Utilities")


test_that("Conf 2 prob", {
    expect_true(length(conf_to_probs(0.90)) == 2)
    expect_equal(conf_to_probs(0.9), c(0.05, 0.95))
})

test_that("Interval funs", {
    expect_true(excludes_zero(c(0.5, 1)))
    expect_true(!excludes_zero(c(-0.5, 1)))
    expect_true(all(excludes_zero(data.frame(c(0.5, 1), c(0.85, 1.5)))))

    expect_true(within_accuracy(c(0.9, 1.1), 1))
    expect_true(!within_accuracy(c(0.5, 1.5), 1))

    expect_true(in_interval(c(0.5, 1.5), 1))
    expect_true(!in_interval(c(0.5, 1.5), 1.75))


})
