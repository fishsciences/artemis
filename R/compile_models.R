##' This compiles the models needed for the artemis package. This must
##' be done before using the \code{eDNA_lm*} functions. This function
##' is called for its side-effects.
##'
##' @title Compile artemis models
##' @param model_names the names of the models to compile
##' @param models full file path to the location of models
##' @param cache_dir a directory to store the Stan model code and
##'     compiled models. Can be set with the option
##'     "artemis_cache_dir", otherwise defaults to the directory
##'     returned by \code{R_user_dir()}. If the directory does not
##'     exist, it will be created.
##' @param rewrite logical, whether to re-compile the models,
##'     overwriting previous versions.
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
                          rewrite = TRUE,
                          verbose = TRUE, ...)
{
    if(!dir.exists(cache_dir)) {
        if(verbose) message("Creating cache directory: ", cache_dir)
        dir.create(cache_dir, recursive = TRUE)
    }
    
    model_files = file.path(cache_dir, model_names)
    is_windows = .Platform$OS.type == "windows"
    out_files = gsub("\\.stan$", ifelse(is_windows, ".exe", ""), model_files)

    if(!all(file.exists(model_files)) || rewrite){
        if(verbose) message("Copying .stan files to cache")
        file.copy(models, model_files, overwrite = rewrite, ...)
    }

    if(!all(file.exists(out_files)) || rewrite){
        if(verbose) message("Compiling executables")
        m = lapply(model_files, cmdstan_model, force_recompile = rewrite, ... )
    }

    if(verbose) message("Models compiled and ready to use!")

    return(invisible(NULL))
}

##' Checks if the pre-compiled models for the artemis package exist.
##'
##' @title Check Compiled Models
##' @param model_names Character string of the base-names of the
##'     compiled models. On Windows platforms, the extension "exe" is
##'     added
##' @param cache_dir directory where the compiled models should be
##'     found.
##' @param issue_error logical, if TRUE an error will be raised if the
##'     models are not all found in the specified cache_dir
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
    is_windows = .Platform$OS.type == "windows"
    out_files = gsub("\\.stan$", ifelse(is_windows, ".exe", ""), model_names)
    out = file.path(cache_dir, out_files)
    models_ok = all(file.exists(out))
    if(issue_error && !models_ok){
        stop("Pre-compiled model file not found. Please check:\n",
             "1. cache_dir is set to proper location\n",
             "2. models have been compiled using compile_models()\n",
             "For more help, see ?compile_models")
    }

    models_ok
}

                          
