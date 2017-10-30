# [ Iteratively reweighted least squares for canonical links - read design matrix one row at a time ]

# `f` is formula; `b`, `b1`, `b1inv` and `b2` are the cumulant function, its
# first derivative, inverse of b', and second derivative; `init`
# is a function that provides an initial 'mu' as a function of 'y'.
# For simplicity, assumes unit dispersion and canonical link

irls <- function (f, b, b1, b1inv, b2, init, tol = 1e-6) {
  # initialize
  y <- model.response(model.frame(f))
  X <- model.matrix(f) # design matrix
  mu <- init(y)
  eta <- b1inv(mu) # g(mu)
  lhood <- sum(y * eta - b(eta))
  
  # iterate
  repeat {
    W <- sapply(eta, b2) # b''(theta) = V(mu)
    z <- eta + (y - mu) / W # working response
    
    ### come up with new beta estimates by reading design matrix one row at a time
    
    XtWX = array(0, dim = c(ncol(X), ncol(X))) # initialize empty p x p matrix - t(X)WX
    XtWz = array(0, dim = c(ncol(X), 1))       # initialize empty p x 1 matrix - t(X)WY
    for (i in 1:nrow(X)) {
      XtWX = XtWX + outer(X[i,] * W[i], X[i,])
      XtWz = XtWz + X[i,] * W[i] * z[i]
    }
    # replace the function below
    
    beta = solve(XtWX) %*% XtWz
    
    eta = array(0, dim = c(nrow(X), 1))
    for (i in 1:nrow(X)) {
      eta[i] = X[i,] %*% beta
    }
    lhood.new <- sum(y * eta - b(eta))
    if (abs((lhood.new - lhood) / lhood) < tol) break # converged?
    lhood <- lhood.new
    mu <- sapply(eta, b1) # b'(theta)
  }
  
  # report
  list(coef = beta, var = solve(XtWX))
}


# Faster IRLS function where both loops through the design matrix are implemented in C ++
library(Rcpp)
sourceCpp("suff_stats.cpp")
sourceCpp("eta_update.cpp")
irls_cpp <- function (f, b, b1, b1inv, b2, init, tol = 1e-6) {
  # initialize
  y <- model.response(model.frame(f))
  X <- model.matrix(f) # design matrix
  mu <- init(y)
  eta <- b1inv(mu) # g(mu)
  lhood <- sum(y * eta - b(eta))
  
  # iterate
  repeat {
    W <- sapply(eta, b2) # b''(theta) = V(mu)
    z <- eta + (y - mu) / W # working response
    
    ### come up with new beta estimates by reading design matrix one row at a time
    
    XtWX = array(0, dim = c(ncol(X), ncol(X))) # initialize empty p x p matrix - t(X)WX
    XtWz = array(0, dim = c(ncol(X), 1))       # initialize empty p x 1 matrix - t(X)WY
    
    ss = v(XtWX, XtWz, X, W, z) # C++ implementation of the loop below
    
    # for (i in 1:nrow(X)) {
    #   XtWX = XtWX + outer(X[i,] * W[i], X[i,])
    #   XtWz = XtWz + X[i,] * W[i] * z[i]
    # }
    
    beta = solve(ss$XtWX) %*% ss$XtWz
    
    eta = eta(beta, X) # C ++ implementation of the loop below

    # eta = array(0, dim = c(nrow(X), 1))
    # for (i in 1:nrow(X)) {
    #   eta[i] = X[i,] %*% beta
    # }
    
    lhood.new <- sum(y * eta - b(eta))
    if (abs((lhood.new - lhood) / lhood) < tol) break # converged?
    lhood <- lhood.new
    mu <- sapply(eta, b1) # b'(theta)
  }
  
  # report
  list(coef = beta, var = solve(ss$XtWX))
}

# e.g. Poisson with canonical link:
poisson.irls <- function (f, tol = 1e-6)
  irls(f, exp, exp, log, exp, function (y) y + 0.5, tol)

poisson.irls.cpp <- function (f, tol = 1e-6)
  irls_cpp(f, exp, exp, log, exp, function (y) y + 0.5, tol)

#################### TEST ####################
# generate poisson response, two predictors
n = 2000
X = cbind(rep(1,n), rnorm(n, mean=0, sd=4))
beta_true = c(0.5,1)
Y = rpois(n=n, lambda = exp(X %*% beta_true))

library(microbenchmark)
poisson.irls(Y~X-1)
poisson.irls.cpp(Y~X-1)
summary(glm(Y~X-1, family = poisson))
microbenchmark(poisson.irls(Y~X-1), poisson.irls.cpp(Y~X-1))
