##' This compiles the models needed for the package. 
##'
##' Compile models.
##' @title Compile models
##' @param model_names the names of the models to compile
##' @param models full file path to the location of models
##' @param rewrite logical, whether to re-compile the models
##' @param ... additional arguments passed to \code{cmdstanr::cmdstan_model()}
##' @return NULL
##' @author Matt Espe
##' @export
compile_models = function(model_names = c("eDNA_lm.stan",
                                          "eDNA_lmer.stan",
                                          "eDNA_sim_omni.stan"),
                          models = system.file("stan", model_names,
                                               package = "artemis"),
                          cache_dir = tools::R_user_dir("artemis", "cache"),
                          rewrite = TRUE,
                          verbose = TRUE, ...)
{
    if(!dir.exists(cache_dir)) {
        if(verbose) message("Creating cache directory: ", cache_dir)
        dir.create(cache_dir, recursive = TRUE)
    }
    
    out_files = file.path(cache_dir, model_names)
    if(verbose) message("Copying .stan files to cache")
    file.copy(models, out_files, overwrite = rewrite, ...)
    if(verbose) message("Compiling executables")
    m = lapply(out_files, cmdstan_model, force_recompile = rewrite, ... )
    if(verbose) message("All done! Models ready for sampling!")
    return(invisible(NULL))
}

