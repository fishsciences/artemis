context("Compile models")

test_that("Compile models", {
    if(is.null(cmdstan_version(error_on_NA = FALSE)))
        skip("Skipping model compilation tests - cmdstan installation not found")

    td = tempdir()
    compile_models(cache_dir = td)
    expect_true(compiled_models_ok(cache_dir = td))

})
