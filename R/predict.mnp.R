predict.mnp <- function(object, newdata = NULL, newdraw = NULL, moredraw = 1,
                        type = c("prob", "choice", "order", "latent"),
                        verbose = FALSE, ...){

  if (NA %in% match(type, c("prob", "choice", "order", "latent")))
    stop("Invalid input for `type'.")
      
  p <- object$n.alt
  if (is.null(newdraw))
    param <- object$param
  else
    param <- newdraw
  coef <- coef(object)
  n.cov <- ncol(coef)
  n.draws <- nrow(param)
  cov <- param[,(n.cov+1):ncol(param)]
  
  ## get X matrix ready
  if (is.null(newdata)) 
    x <- object$x
  else {
    call <- object$call
    x <- xmatrix.mnp(as.formula(call$formula), data = newdata,
                     choiceX = call$choiceX,
                     cXnames = eval(call$cXnames),
                     base = object$base, n.dim = p-1,
                     lev = object$alt, MoP = is.matrix(object$y),
                     verbose = FALSE, extra = FALSE)
    if (nrow(x) > 1) 
      x <- as.matrix(x[apply(is.na(x), 1, sum)==0,])
    else if (sum(is.na(x))>0)
      stop("Invalid input for `newdata'.")
  }
  
  n.obs <- nrow(x)
  if (verbose) {
    if (n.obs == 1)
      cat("There is one observation to predict. Please wait...\n")
    else
      cat("There are", n.obs, "observations to predict. Please wait...\n")
  }
  
  alt <- object$alt
  if (object$base != alt[1]) 
    alt <- c(object$base, alt[1:(length(alt)-1)])
  
  ## computing W
  W <- array(NA, dim=c(p-1, n.obs, n.draws, moredraw),
             dimnames=c(alt[2:p], NULL, 1:n.draws, 1:moredraw))
  tmp <- floor(n.draws/10)
  inc <- 1
  Sigma <- cov.mnp(object)
  for (i in 1:n.draws) {
    for (j in 1:n.obs) 
      W[,j,i,] <- matrix(x[j,], ncol=n.cov) %*% matrix(coef[i,]) +
        mvrnorm(moredraw, mu = rep(0, p-1), Sigma = Sigma[,,i])
    if (i == inc*tmp & verbose) {
      cat("", inc*10, "percent done.\n")
      inc <- inc + 1
    }
  }
  ans <- list()
  if ("latent" %in% type)
    ans$w <- W
  else
    ans$w <- NULL

  ## computing Y
  Y <- array(NA, dim = c(n.obs, n.draws, moredraw),
             dimnames=list(NULL, 1:n.draws, 1:moredraw))
  O <- array(NA, dim=c(p, n.obs, n.draws, 1:moredraw),
             dimnames=list(alt, NULL, 1:n.draws, 1:moredraw))
  for (i in 1:n.obs) 
    for (j in 1:n.draws) {
      for (k in 1:moredraw) {
        Y[i,j,k] <- alt[match(max(c(0, W[,i,j,k])), c(0,W[,i,j,k]))]
        O[,i,j,k] <- rank(c(0, -W[,i,j,k]))
      }
    }
  if ("choice" %in% type)
    ans$y <- Y
  else
    ans$y <- NULL
  if ("order" %in% type)
    ans$o <- O
  else
    ans$o <- NULL

  ## computing P
  if ("prob" %in% type) {
    P <- array(NA, dim = c(n.obs, p, moredraw),
               dimnames = c(NULL, alt, 1:moredraw))
    for (i in 1:p)
      for (j in 1:moredraw)
        P[,i,j] <- apply(Y[,,j,drop=FALSE]==alt[i], 1, mean) 
    ans$p <- P
  }
  else
    ans$p <- NULL

  return(ans)
}
