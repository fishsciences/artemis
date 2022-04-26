library(testthat)
library(artemis)

# Checks for cmdstan installed - if not installed or path not set,
# does not run tests
if(!is.null(cmdstan_version(error_on_NA = FALSE))){
    if(compiled_models_ok())
        test_check("artemis")
}
