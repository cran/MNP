print.summary.mnp <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
      sep = "") 

  cat("\nCoefficients:\n")
  printCoefmat(x$coef.table, digits = digits, na.print = "NA", ...)
  
  cat("\nCovariances:\n")
  printCoefmat(x$cov.table, digits = digits, na.print = "NA", ...)
  
  cat("\nNumber of alternatives:", x$n.alt)
  cat("\nNumber of observations:", x$n.obs)
  cat("\n\n")
  invisible(x)
}
