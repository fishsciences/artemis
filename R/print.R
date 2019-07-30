print.eDNA_simulation = function(x, FUN = summary, digits = getOption("digits")) {
    cat("\nformula: "); print(x@formula)
    cat("\nvariable levels:\n")
    
    vars = sapply(seq_along(x@variable_levels), function(i)
        paste0(names(x@variable_levels)[i],  " :",
               paste(x@variable_levels[[i]], collapse = " ")))
    
    cat(vars, sep = "\n")
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
