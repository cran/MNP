summary.mnp <- function(object, CI=c(2.5, 97.5),...){

  p <- object$n.alt
  param <- object$param
  n.cov <- ncol(param) - p*(p-1)/2
  n.draws <- nrow(param)
  param.table <- cbind(apply(param, 2, mean), apply(param, 2, sd),
                       apply(param, 2, quantile, min(CI)/100),
                       apply(param, 2, quantile, max(CI)/100)) 
  colnames(param.table) <- c("mean", "std.dev.", paste(min(CI), "%", sep=""),
                             paste(max(CI), "%", sep=""))
  
  ans <- list(call=object$call, base = object$base, n.alt=p, n.obs=if(is.matrix(object$y))
              nrow(object$y) else length(object$y), n.draws = n.draws,
              coef.table=param.table[1:n.cov,],
              cov.table=param.table[(n.cov+1):ncol(param),])  
  class(ans) <- "summary.mnp"
  return(ans)
}
