predict.mnp <- function(object, newdata = NULL, newdraw = NULL,
                        type = c("prob", "choice", "order", "latent"),
                        verbose = FALSE, ...){

  if (NA %in% match(type, c("prob", "choice", "order", "latent")))
    stop("Invalid input for `type'.")
      
  p <- object$n.alt
  if (is.null(newdraw))
    param <- object$param
  else
    param <- newdraw
  n.cov <- ncol(param) - p*(p-1)/2
  n.draws <- nrow(param)
  coef <- param[,1:n.cov]
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
  }

  n.obs <- nrow(x)
  if (verbose)
    cat("There are", n.obs, "observations to predict. Please wait...\n")

  alt <- object$alt
  if (object$base != alt[1]) 
    alt <- c(object$base, alt[1:(length(alt)-1)])
  
  ## computing W
  W <- array(NA, dim=c(p-1, n.obs, n.draws), dimnames=c(alt[2:p],
                                               NULL, 1:n.draws))
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
  Y <- matrix(NA, nrow = n.obs, ncol = n.draws, dimnames=list(NULL, 1:n.draws))
  O <- array(NA, dim=c(p, n.obs, n.draws), dimnames=list(alt, NULL, 1:n.draws))
  for (i in 1:n.obs) 
    for (j in 1:n.draws) {
      Y[i,j] <- alt[match(max(c(0, W[,i,j])), c(0,W[,i,j]))]
      O[,i,j] <- rank(c(0, -W[,i,j]))
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
    P <- matrix(NA, nrow = n.obs, ncol = p)
    colnames(P) <- alt
    for (i in 1:p)
      P[,i] <- apply(Y==alt[i], 1, mean) 
    ans$p <- P
  }
  else
    ans$p <- NULL

  return(ans)
}