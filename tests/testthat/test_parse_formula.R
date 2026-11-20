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
    expect_true(length(d) == 7)
    expect_true(all(names(d) %in% c("y", "x", "rand_x", "groups", "group_names", "n_rand", "n_grp")))
    expect_true(length(d$groups) == ncol(d$rand_x))
    expect_true(ncol(d$x) == 2)

    expect_equal(d$group_names, "(Intercept):cyl")
    expect_true(colnames(d$rand_x)[1] == "(Intercept):cyl[4]")

    d = gen_model_list_lmer(mpg ~ cyl + am + (1|cyl/am), mtcars)

    expect_is(d, "list")
    expect_true(length(d) == 7)
    expect_true(all(names(d) %in% c("y", "x", "rand_x", "groups", "group_names", "n_rand", "n_grp")))
    expect_true(length(d$groups) == ncol(d$rand_x))
    expect_true(ncol(d$x) == 3)

    d = gen_model_list_lmer(mpg ~ factor(cyl) + factor(am) + (1|cyl/am), mtcars)

    expect_is(d, "list")
    expect_true(length(d) == 7)
    expect_true(all(names(d) %in% c("y", "x", "rand_x", "groups", "group_names", "n_rand", "n_grp")))
    expect_true(length(d$groups) == ncol(d$rand_x))
    expect_true(ncol(d$x) == 4)

    # observation-level random effect (one level per row)
    d = mtcars
    d$obs = seq_len(nrow(d))
    d = gen_model_list_lmer(mpg ~ cyl + (1|obs), d)
    expect_true(ncol(d$rand_x) == nrow(mtcars))
})

test_that("Random effect columns are contiguous by group and match their names", {
    set.seed(1)
    d = data.frame(y = rnorm(24), x = rnorm(24),
                   site = gl(3, 8), filt = gl(4, 1, 24))
    ml = gen_model_list_lmer(y ~ x + (1 + x | site) + (1 | filt), d)
    nm = unname(colnames(ml$rand_x))

    # Stan's segmented sum-to-zero needs each group in one block
    expect_false(is.unsorted(ml$groups))
    expect_equal(sub("\\[.*", "", nm), ml$group_names[ml$groups])
    expect_equal(sort(ml$group_names), c("(Intercept):filt", "(Intercept):site", "x:site"))
    sizes = setNames(tabulate(ml$groups), ml$group_names)
    expect_equal(unname(sizes[c("(Intercept):filt", "(Intercept):site", "x:site")]),
                 c(4L, 3L, 3L))

    # each column holds the values its name describes
    term = sub(":.*", "", nm)
    fac = sub(".*:(.*)\\[.*", "\\1", nm)
    lev = sub(".*\\[(.*)\\]", "\\1", nm)
    for(j in seq_along(nm)){
        val = if(term[j] == "(Intercept)") 1 else d[[term[j]]]
        expect_equal(unname(ml$rand_x[, j]), val * (d[[fac[j]]] == lev[j]),
                     info = nm[j])
    }

    # separate terms on one factor are already contiguous
    ml2 = gen_model_list_lmer(y ~ x + (1 | site) + (0 + x | site), d)
    expect_false(is.unsorted(ml2$groups))
    expect_equal(sub("\\[.*", "", unname(colnames(ml2$rand_x))), ml2$group_names[ml2$groups])
})

test_that("Error on NA", {
    d = mtcars
    d$mpg[3] = NA
    expect_error(gen_model_list_lmer(mpg ~ cyl + (1|cyl), d))
    expect_error(gen_model_list_lmer(mpg ~ cyl, d))
})

test_that("Zero-inflated mixed formula requires a bar", {
    expect_error(gen_model_list_lmer_zip(Cq ~ Distance_m + (1|FilterID), eDNA_data),
                 "zero-inflated")
    d = eDNA_data
    d$Volume_mL[1] = NA
    expect_error(gen_model_list_lmer_zip(Cq ~ Distance_m + (1|FilterID) | Volume_mL, d),
                 "NA values")
    ans = gen_model_list_lmer_zip(Cq ~ Distance_m + (1|FilterID) | Volume_mL, eDNA_data)
    expect_true(all(c("x", "xz", "rand_x") %in% names(ans)))
})
