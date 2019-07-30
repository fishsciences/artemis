setMethod("print", "eDNA_simulation",
          function(x, FUN = summary, digits = 3) {
              cat("\nformula: "); print(x@formula)
              cat("\nvariable levels:\n")
              
              vars = sapply(seq_along(x@variable_levels), function(i)
                  paste0(names(x@variable_levels)[i],  " :",
                         paste(x@variable_levels[[i]], collapse = " ")))

              cat(vars, sep = "\n")
              cat("\n ln concentration: \n")
              print(get_marginals(x@ln_conc, x@x, FUN, digits))

              cat("\n simulated Cq: \n")
              print(get_marginals(x@Cq_star, x@x, FUN, digits))
              
              return(invisible(NULL))
          })
