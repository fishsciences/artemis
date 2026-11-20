##' This compiles the artemis models into a cache directory. Most users
##' do not need to call it. The models are compiled when the package is
##' installed, as long as CmdStan is already installed at that time.
##' When those compiled models are missing, the modeling and simulation
##' functions call this function themselves on first use. This function
##' is called for its side-effects.
##'
##' A cached \code{.stan} file that differs from the copy in the
##' installed package is replaced and recompiled. This keeps the cache
##' current after the package is updated. Use \code{rewrite = TRUE} to
##' force a rebuild of the cached models. To rebuild the models compiled
##' at install (e.g. after CmdStan is updated), reinstall artemis.
##'
##' @title Compile artemis models
##' @param model_names the names of the models to compile
##' @param models full file path to the location of models
##' @param cache_dir a directory to store the Stan model code and
##'     compiled models. Can be set with the option
##'     "artemis_cache_dir", otherwise defaults to the directory
##'     returned by \code{R_user_dir()}. If the directory does not
##'     exist, it will be created.
##' @param rewrite logical, whether to re-compile all models,
##'     overwriting previous versions even when they are current.
##' @param verbose logical, if TRUE messages will report steps being
##'     taken
##' @param ... additional arguments passed to
##'     \code{cmdstanr::cmdstan_model()}
##' @return NULL
##' @author Matt Espe
##' @export
compile_models = function(model_names = c("eDNA_lm.stan",
                                          "eDNA_lmer.stan",
                                          "eDNA_lm_zinf.stan",
                                          "eDNA_lmer_zinf.stan",
                                          "eDNA_pois_lm.stan",
                                          "eDNA_pois_lmer.stan",
                                          "eDNA_sim_omni.stan"),
                          models = system.file("stan", model_names,
                                               package = "artemis"),
                          cache_dir = getOption("artemis_cache_dir", R_user_dir("artemis", "cache")),
                          rewrite = FALSE,
                          verbose = TRUE, ...)
{
    if(package_version(cmdstan_version()) < "2.37.0")
        stop("artemis requires CmdStan 2.37.0 or newer; found ", cmdstan_version(),
             ". Update with cmdstanr::install_cmdstan()")

    if(!dir.exists(cache_dir)) {
        if(verbose) message("Creating cache directory: ", cache_dir)
        dir.create(cache_dir, recursive = TRUE)
    }
    
    model_files = file.path(cache_dir, model_names)
    is_windows = .Platform$OS.type == "windows"
    out_files = gsub("\\.stan$", ifelse(is_windows, ".exe", ""), model_files)

    source_md5 = unname(md5sum(models))
    if(anyNA(source_md5))
        stop("Stan model files not found: ", paste(models[is.na(source_md5)], collapse = ", "))

    # refresh cached code that differs from the installed package
    cached = unname(md5sum(model_files))
    stale = rewrite | is.na(cached) | cached != source_md5
    if(any(stale)){
        if(verbose) message("Copying .stan files to cache")
        file.copy(models[stale], model_files[stale], overwrite = TRUE)
    }

    to_compile = stale | !file.exists(out_files)
    if(any(to_compile)){
        if(verbose) message("Compiling executables")
        m = lapply(model_files[to_compile], cmdstan_model, force_recompile = TRUE, ...)
    }

    if(verbose) message("Models compiled and ready to use!")

    return(invisible(NULL))
}

##' Checks if the pre-compiled models for the artemis package exist,
##' either from the package installation or in the cache directory.
##'
##' @title Check Compiled Models
##' @param model_names Character string of the base-names of the
##'     compiled models. On Windows platforms, the extension "exe" is
##'     added
##' @param cache_dir directory where the cached compiled models should
##'     be found.
##' @param issue_error logical, if TRUE an error will be raised if the
##'     models are not all found
##' @return logical, TRUE if the compiled models were found and FALSE
##'     otherwise
##' @author Matt Espe
##' @export
##' @examples{
##' compiled_models_ok()
##' }
compiled_models_ok = function(model_names = c("eDNA_lm",
                                              "eDNA_lmer",
                                              "eDNA_sim_omni"),
                              cache_dir = getOption("artemis_cache_dir", R_user_dir("artemis", "cache")),
                              issue_error = FALSE)
{
    exe_names = exe_name(model_names)
    models_ok = all(file.exists(installed_exe(model_names)) |
                    file.exists(file.path(cache_dir, exe_names)))
    if(issue_error && !models_ok){
        stop("Pre-compiled model file not found. Please check:\n",
             "1. CmdStan was installed before artemis, or\n",
             "2. models have been compiled using compile_models()\n",
             "For more help, see ?compile_models")
    }

    models_ok
}

# model file name ("eDNA_lm" or "eDNA_lm.stan") to executable name
exe_name = function(model_name, is_windows = .Platform$OS.type == "windows")
{
    paste0(sub("\\.stan$", "", model_name), if(is_windows) ".exe" else "")
}

# path to the executable compiled at install time (may not exist)
installed_exe = function(model_name)
{
    file.path(system.file(package = "artemis"), "bin", "stan", exe_name(model_name))
}

# CmdStanModel for a model file, preferring the install-time executable
# and falling back to compiling into the cache
get_model = function(model_name, cache_dir, exe = installed_exe(model_name))
{
    if(file.exists(exe))
        return(cmdstan_model(exe_file = exe))
    compile_models(model_name, cache_dir = cache_dir, verbose = FALSE)
    cmdstan_model(file.path(cache_dir, model_name))
}

                          
