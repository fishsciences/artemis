context("Formula parsing")

test_that("Choosing correct parser", {
    expect_true(is_lme4(mpg ~ cyl + (1|cyl)))
    expect_true(!is_lme4(mpg~cyl))
    expect_true(!is_lme4(mpg~factor(cyl)))

    expect_true(is_lme4(mpg ~ factor(cyl) + (1|cyl)))
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

    d = gen_model_list_lm( ~ cyl + wt, mtcars)
    
})

test_that("lmer model.list", {
    d = gen_model_list_lmer(mpg ~ cyl + (1|cyl), mtcars)

    expect_is(d, "list")
    expect_true(length(d) == 6)
    expect_true(all(names(d) %in% c("y", "x", "rand_x", "groups", "n_rand", "n_grp")))
    expect_true(length(d$groups) == ncol(d$rand_x))
    expect_true(ncol(d$x) == 2)

    d = gen_model_list_lmer(mpg ~ cyl + am + (1|cyl/am), mtcars)

    expect_is(d, "list")
    expect_true(length(d) == 6)
    expect_true(all(names(d) %in% c("y", "x", "rand_x", "groups", "n_rand", "n_grp")))
    expect_true(length(d$groups) == ncol(d$rand_x))
    expect_true(ncol(d$x) == 3)

    d = gen_model_list_lmer(mpg ~ factor(cyl) + factor(am) + (1|cyl/am), mtcars)

    expect_is(d, "list")
    expect_true(length(d) == 6)
    expect_true(all(names(d) %in% c("y", "x", "rand_x", "groups", "n_rand", "n_grp")))
    expect_true(length(d$groups) == ncol(d$rand_x))
    expect_true(ncol(d$x) == 4)

})

test_that("Error on NA", {
    d = mtcars
    d$mpg[3] = NA
    expect_error(gen_model_list_lmer(mpg ~ cyl + (1|cyl), d))
    expect_error(gen_model_list_lmer(mpg ~ cyl, d))
})
