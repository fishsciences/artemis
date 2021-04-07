library(testthat)
library(artemis)

if(length(cmdstan_path) == 0)
    install_cmdstan()

test_check("artemis")
