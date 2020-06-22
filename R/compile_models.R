##' This compiles the models needed for the package. 
##'
##' Compile models.
##' @title Compile models
##' @param model_names the names of the models to compile
##' @param models full file path to the location of models
##' @param rewrite logical, whether to re-compile the models
##' @param ... additional arguments passed to \code{cmdstan_model}
##' @return NULL
##' @author Matt Espe
##' @export
compile_models = function(model_names = c("eDNA_omni.stan", "eDNA_sim_omni.stan"),
                          models = system.file("stan_files", model_names,
                                               package = "artemis"),
                          rewrite = TRUE, ...)
{
    m = lapply(models, cmdstan_model, force_recompile = rewrite, ... )
    return(invisible(NULL))
    
}

