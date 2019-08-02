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
