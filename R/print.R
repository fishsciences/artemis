print.eDNA_simulation = function(x, FUN = summary, digits = getOption("digits"),
                                 show_variables = FALSE) {
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
    print(summary(x, var = "ln_conc", digits = digits))

    cat("\n simulated Cq: \n")
    print(summary(x, digits = digits))    
    return(invisible(NULL))
}

print.eDNA_simulation.summary = function(x, digits = getOption("digits"), ...)
{
    x = as.data.frame(sapply(x,
                             function(y, digs)
                                 if(is.numeric(y)) round(y, digs) else y
                           , digs = digits),
                      stringsAsFactors = FALSE)
    NextMethod()
}

print.eDNA_model = function(x, digits = getOptions("digits"), ...)
{
    cat("\nformula: "); print(x@formula)
    cat("\nStandard curve parameters: Cq = alpha + beta * log(concentration)\n")
    cat("\tStandard curve alpha = ", x@std_curve_alpha, "\n")
    cat("\tStandard curve beta = ", x@std_curve_beta, "\n")
    invisible(x)
}



print.eDNA_model.summary = function(x, digits = getOptions("digits"), ...)
{
    NextMethod()
    invisible(x)
}
