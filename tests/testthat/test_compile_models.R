context("Compile models")

test_that("Compile models", {
    if(is.null(cmdstan_version(error_on_NA = FALSE)))
        skip("Skipping model compilation tests - cmdstan installation not found")

    td = tempdir()
    compile_models(cache_dir = td)
    expect_true(compiled_models_ok(cache_dir = td))

    # an outdated cached model is replaced
    mn = "eDNA_sim_omni.stan"
    cat("// outdated\n", file = file.path(td, mn), append = TRUE)
    compile_models(mn, cache_dir = td, verbose = FALSE)
    expect_equal(unname(tools::md5sum(file.path(td, mn))),
                 unname(tools::md5sum(system.file("stan", mn, package = "artemis"))))
})
