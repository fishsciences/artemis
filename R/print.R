##' Print eDNA simulation
##'
##' Print method.
##' 
##' @title Print eDNA simulation
##' @param x object of class eDNA_simulation
##' @param FUN a function to use to summarize the results, default is \code{summary}
##' @param digits number of digits to show 
##' @param show_variables logical, should the variable levels used in the simulation
##'     be included
##' @param ... additional arguments passed to \code{print}
##' @return NULL
##' @author Matt Espe
##' @method print eDNA_simulation
##' @export
print.eDNA_simulation = function(x, FUN = summary, digits = getOption("digits"),
                                 show_variables = FALSE, ...) {
    cat("\nformula: "); print(x@formula)

    if(show_variables){
    cat("\nvariable levels:\n")
    
    vars = sapply(seq_along(x@variable_levels), function(i)
        paste0(names(x@variable_levels)[i],  " :",
               paste(x@variable_levels[[i]], collapse = " ")))
    
    cat(vars, sep = "\n")
    }
    cat("\nStandard curve parameters: Cq = alpha + beta * log(concentration)\n")
    cat("\tStandard curve alpha = ", x@std_curve_alpha, "\n")
    cat("\tStandard curve beta = ", x@std_curve_beta, "\n")
    cat("\n ln concentration: \n")
    print(summary(x, var = "ln_conc"), digits = digits, row.names = FALSE, ...)

    cat("\n simulated Cq: \n")
    print(summary(x), digits = digits, row.names = FALSE, ...) 
    invisible(x)
}

##' Print eDNA model results
##'
##' Print method.
##'
##' @title Print eDNA model results
##' @param x object of class "eDNA_model_*"
##' @param digits number of digits to show 
##' @param ... additional arguments passed to \code{print}
##' @return NULL
##'
##' @author Matt Espe
##' @method print eDNA_model 
##' @export
print.eDNA_model = function(x, digits = getOption("digits"), ...)
{
    cat("\nformula: "); print(x@formula)
    cat("\nStandard curve parameters: Cq = alpha + beta * log(concentration)\n")
    cat("\tStandard curve alpha = ", x@std_curve_alpha, "\n")
    cat("\tStandard curve beta = ", x@std_curve_beta, "\n\n")

    cat("\nParameter estimates:\n")
    print(summary(x, ...), digits = digits, ...)
    invisible(x)
}

##' Print method for p(detect)
##'
##' Print method for p(detect)
##' @title Print eDNA p(detect)
##' @param x object of class "eDNA_p_detect"
##' @param digits number of digits to show 
##' @param ... additional arguments passed to \code{print}
##' @return x 
##' @author Matt Espe
##' @method print eDNA_p_detect
##' @export
print.eDNA_p_detect = function(x, digits = getOption("digits"), ...)
{
    cat("Variable levels: \n")
    print(attr(x, "variable_levels"), digits, ...)
    cat("\n")
    tmp = if(is.matrix(x)){
              summary(x)
          } else {
              data.frame(n_reps = attr(x, "reps"),
                         p_detect = as.numeric(x))
    }
    
    print.default(tmp, digits, ...)

    invisible(x)          
}
