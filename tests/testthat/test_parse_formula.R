context("Formula parsing")

test_that("Choosing correct parser", {
    expect_true(is_lme4(mpg ~ cyl + (1|cyl)))
    expect_true(!is_lme4(mpg~cyl))

})

test_that("lm model.list", {
    d = gen_model_list_lm(mpg ~ cyl + wt, mtcars)

    expect_is(d, "list")
    expect_true(length(d) == 2)

    expect_true(all(names(d) %in% c("x", "y")))
    expect_true(ncol(d$x) == 3)

    # Test changes in formulas
    d = gen_model_list_lm(mpg ~ as.factor(cyl) + I(sqrt(wt)), mtcars)

    expect_is(d, "list")
    expect_true(length(d) == 2)

    expect_true(all(names(d) %in% c("x", "y")))
    expect_true(ncol(d$x) == 4)

})

test_that("lmer model.list", {
    d = gen_model_list_lmer(mpg ~ cyl + (1|cyl), mtcars)

    expect_is(d, "list")
    expect_true(length(d) == 5)
    expect_true(all(names(d) %in% c("x","y","groups", "n_grp", "n_levels")))
    expect_true(ncol(d$groups) == 1)
    expect_true(ncol(d$x) == 2)

    d = gen_model_list_lmer(mpg ~ cyl + am + (1|cyl/am), mtcars)

    expect_is(d, "list")
    expect_true(length(d) == 5)
    expect_true(all(names(d) %in% c("x","y","groups", "n_grp", "n_levels")))
    expect_true(ncol(d$groups) == 2)
    expect_true(ncol(d$x) == 3)

    d = gen_model_list_lmer(mpg ~ factor(cyl) + factor(am) + (1|cyl/am), mtcars)

    expect_is(d, "list")
    expect_true(length(d) == 5)
    expect_true(all(names(d) %in% c("x","y","groups", "n_grp", "n_levels")))
    expect_true(ncol(d$groups) == 2)
    expect_true(ncol(d$x) == 4)

})
