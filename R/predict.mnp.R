predict.mnp <- function(object, newdata = NULL,
                        type = c("prob", "choice", "latent"),
                        verbose = FALSE, ...){

  if (NA %in% match(type, c("prob", "choice", "latent")))
    stop("Invalid input for `type'.")
      
  p <- object$n.alt
  param <- object$param
  n.cov <- ncol(param) - p*(p-1)/2
  n.draws <- nrow(param)
  coef <- param[,1:n.cov]
  cov <- param[,(n.cov+1):ncol(param)]
  alt <- object$alt
  
  ## get X matrix ready
  if (is.null(newdata)) 
    x <- object$x
  else {
    call <- object$call
    x <- xmatrix.mnp(as.formula(call$formula), data = newdata,
                     choiceX = call$choiceX,
                     cXnames = call$cXnames,
                     base = object$base, n.dim = p-1,
                     lev = object$alt, MoP = is.matrix(object$y),
                     verbose = FALSE, extra = FALSE)    
  }

  n.obs <- nrow(x)
  if (verbose)
    cat("There are", n.obs, "observations to predict. Please wait...\n")
  ## computing W
  W <- array(NA, dim=c(p-1, n.obs, n.draws))
  tmp <- floor(n.draws/10)
  inc <- 1
  for (i in 1:n.draws) {
    Sigma <- matrix(0, p-1, p-1)
    count <- 1
    for (j in 1:(p-1)) {
      Sigma[j,j:(p-1)] <- cov[i,count:(count+p-j-1)]
      count <- count + p - j
    }
    diag(Sigma) <- diag(Sigma)/2
    Sigma <- Sigma + t(Sigma)
    for (j in 1:n.obs) 
      W[,j,i] <- matrix(x[j,], ncol=n.cov) %*% matrix(coef[i,]) +
        mvrnorm(1, mu = rep(0, p-1), Sigma = Sigma)
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
  Y <- matrix(NA, nrow = n.obs, ncol = n.draws)
  for (i in 1:n.obs) 
    for (j in 1:n.draws)
      Y[i,j] <- alt[match(max(c(0, W[,i,j])), c(0,W[,i,j]))]
  if ("choice" %in% type)
    ans$y <- Y
  else
    ans$y <- NULL

  ## computing P
  P <- matrix(NA, nrow = n.obs, ncol = p)
  colnames(P) <- alt
  for (i in 1:p)
    P[,i] <- apply(Y==alt[i], 1, mean) 
  if ("prob" %in% type)
    ans$p <- P
  else
    ans$p <- NULL

  return(ans)
}
