mprobit <- function(formula, data = parent.frame(), choiceX = NULL,
                    cXnames = NULL, base = NULL, n.draws = 5000,
                    p.var = "Inf", p.df = n.dim+1, p.scale = 1,
                    p.alpha0 = 1, coef.start = 0, cov.start = 1,
                    burnin = 0, thin = 0, verbose = FALSE) {  
  call <- match.call()
  mf <- match.call(expand = FALSE)
  mf$choiceX <- mf$cXnames <- mf$base <- mf$n.draws <- mf$p.var <-
    mf$p.df <- mf$p.scale <- mf$p.alpha0 <- mf$coef.start <-
      mf$cov.start <- mf$verbose <- mf$burnin <- mf$thin <- NULL  
  mf[[1]] <- as.name("model.frame.default")
  mf$na.action <- 'na.pass'
  mf <- eval.parent(mf)

  ## obtaining Y
  tmp <- ymatrix.mnp(mf, base=base, extra=TRUE, verbose=verbose)
  Y <- tmp$Y
  MoP <- tmp$MoP
  lev <- tmp$lev
  base <- tmp$base
  p <- tmp$p
  n.dim <- p - 1
  if(verbose)
    cat("\nThe base category is `", base, "'.\n\n", sep="") 
  if (p < 3)
    stop(paste("Error: The number of alternatives should be at least 3."))
  if(verbose) 
    cat("The total number of alternatives is ", p, ".\n\n", sep="") 
  
  ### obtaining X
  tmp <- xmatrix.mnp(formula, data=eval.parent(data),
                     choiceX=call$choiceX, cXnames=cXnames, 
                     base=base, n.dim=n.dim, lev=lev, MoP=MoP,
                     verbose=verbose, extra=TRUE)
  X <- tmp$X
  coefnames <- tmp$coefnames
  n.cov <- ncol(X) / n.dim
  n.obs <- nrow(X)
  if (verbose)
    cat("The dimension of beta is ", n.cov, ".\n\n", sep="")

  ## checking the prior for beta
  p.imp <- FALSE 
  if (p.var == Inf) {
    p.imp <- TRUE
    p.prec <- diag(0, n.cov)
    if (verbose)
      cat("Improper prior will be used for beta.\n\n")
  }
  else if (is.matrix(p.var)) {
    if (ncol(p.var) != n.cov || nrow(p.var) != n.cov)
      stop("Error: The dimension of `p.var' should be ", n.cov, " x ", n.cov, sep="")
    if (sum(sign(eigen(p.var)$values) < 1) > 0)
      stop("Error: `p.var' must be positive definite.")
    p.prec <- solve(p.var)
  }
  else {
    p.var <- diag(p.var, n.cov)
    p.prec <- solve(p.var)
  }
  p.mean <- rep(0, n.cov)

  ## checking prior for Sigma
  p.df <- eval(p.df)
  if (length(p.df) > 1)
    stop(paste("Error: `p.df' must be a positive integer."))
  if (p.df < n.dim)
    stop(paste("Error: `p.df' must be at least ", n.dim, ".", sep=""))
  if (abs(as.integer(p.df) - p.df) > 0)
    stop(paste("Error: `p.df' must be a positive integer."))
  if (!is.matrix(p.scale))  
    p.scale <- diag(p.scale, n.dim)
  if (ncol(p.scale) != n.dim || nrow(p.scale) != n.dim)
    stop("Error: `p.scale' must be ", n.dim, " x ", n.dim, sep="")
  if (sum(sign(eigen(p.scale)$values) < 1) > 0)
    stop("Error: `p.scale' must be positive definite.")
  else if (p.scale[1,1] != 1) {
    p.scale[1,1] <- 1
    warning("p.scale[1,1] will be set to 1.")
  }
  if (p.alpha0 <= 0)
    stop("Error: `p.alpha0' must be positive scalar")
  Signames <- NULL
  for(j in 1:n.dim)
    for(k in 1:n.dim)
      if (j<=k)
        Signames <- c(Signames, paste(if(MoP) lev[j] else lev[j+1],
                                      ":", if(MoP) lev[k] else lev[k+1], sep="")) 

  ## checking starting values
  if (length(coef.start) == 1)
    coef.start <- rep(coef.start, n.cov)
  else if (length(coef.start) != n.cov)
    stop(paste("Error: The dimenstion of `coef.start' must be  ",
               n.cov, ".", sep=""))
  if (!is.matrix(cov.start)) {
    cov.start <- diag(n.dim)*cov.start
    cov.start[1,1] <- 1
  }
  else if (ncol(cov.start) != n.dim || nrow(cov.start) != n.dim)
    stop("Error: The dimension of `cov.start' must be ", n.dim, " x ", n.dim, sep="")
  else if (sum(sign(eigen(cov.start)$values) < 1) > 0)
    stop("Error: `cov.start' must be a positive definite matrix.")
  else if (cov.start[1,1] != 1) {
    cov.start[1,1] <- 1
    warning("cov.start[1,1] will be set to 1.")
  }
  
  ## checking thinnig and burnin intervals
  if (burnin < 0)
    stop("Error: `burnin' should be a non-negative integer.") 
  if (thin < 0)
    stop("Error: `thin' should be a non-negative integer.")
  keep <- thin + 1
  
  ## running the algorithm
  n.par <- n.cov + n.dim*(n.dim+1)/2
  if(verbose)
    cat("Starting Gibbs sampler...\n")
  # recoding NA into -1
  Y[is.na(Y)] <- -1 
  param <- .C("cMNPgibbs", as.integer(n.dim),
              as.integer(n.cov), as.integer(n.obs), as.integer(n.draws),
              as.double(p.mean), as.double(p.prec), as.integer(p.df),
              as.double(p.scale*p.alpha0), as.double(X), as.integer(Y), 
              as.double(coef.start), as.double(cov.start), 
              as.integer(p.imp), as.integer(burnin), as.integer(keep), 
              as.integer(verbose), as.integer(MoP),
              pdStore = double(n.par*(ceiling((n.draws-burnin)/keep)+1)),
              PACKAGE="MNP")$pdStore 
  param <- matrix(param, ncol=n.par,
                  nrow=(ceiling((n.draws-burnin)/keep)+1), byrow=TRUE)
  colnames(param) <- c(coefnames, Signames)

  ##recoding -1 back into NA
  Y[Y==-1] <- NA
  ## returning the object
  res <- list(param =param, x = X, y = Y, call = call, n.alt = p,
              p.mean = if(p.imp) NULL else p.mean, p.var = p.var,
              p.df = p.df, p.scale = p.scale, p.alpha0 = p.alpha0,
              burnin = burnin, thin = thin, seed = .Random.seed)
  class(res) <- "mnp"
  return(res)
}
  


