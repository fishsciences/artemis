# Compile the artemis Stan models into the installed package during R CMD INSTALL
# Matt Espe
# November 2026
#
# Adapted from the install.libs.R template in the instantiate package by
# William Michael Landau (copyright Eli Lilly and Company, MIT license),
# https://github.com/wlandau/instantiate. As in instantiate, the .stan files
# are copied to <package>/bin/stan and compiled there with cmdstanr. Unlike
# instantiate, a skipped or failed compile never stops the install: artemis
# then compiles the models into the user cache on first use.

min_cmdstan = "2.37.0"
bin_stan = file.path(R_PACKAGE_DIR, "bin", "stan")
dir.create(bin_stan, recursive = TRUE, showWarnings = FALSE)
stan_files = list.files(file.path("..", "inst", "stan"), pattern = "\\.stan$",
                        full.names = TRUE)
file.copy(stan_files, bin_stan, overwrite = TRUE)

v = if(requireNamespace("cmdstanr", quietly = TRUE))
        cmdstanr::cmdstan_version(error_on_NA = FALSE)
# pkgbuild::compile_dll() (load_all, document, test) reinstalls into a
# devtools_install_* temp dir on every call; skip the slow compile there
if(grepl("devtools_install_", R_PACKAGE_DIR, fixed = TRUE)) {
    NULL
} else if(is.null(v) || package_version(v) < min_cmdstan) {
    message("CmdStan ", min_cmdstan, " or newer not found: ",
            "artemis models will be compiled on first use")
} else {
    for(f in file.path(bin_stan, basename(stan_files))) {
        message("Compiling ", basename(f))
        tryCatch(cmdstanr::cmdstan_model(f, force_recompile = TRUE),
                 error = function(e)
                     message("Compiling ", basename(f), " failed; ",
                             "it will be compiled on first use:\n",
                             conditionMessage(e)))
    }
}
