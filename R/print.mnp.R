print.mnp <- function (x, digits = max(3, getOption("digits") - 3), ...)
  {
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    param <- apply(x$param, 2, mean)
    if (length(param)) {
      cat("Parameters:\n")
      print.default(format(param, digits = digits), print.gap = 2,
                    quote = FALSE)
    }
    else cat("No parameters\n")
    cat("\n")
    invisible(x)
  }
