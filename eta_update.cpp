// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// C ++ implementation of eta update

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
colvec eta(colvec beta, mat X) {
  int n = size(X, 0);
  vec eta = ones<vec>(n);
  for(int i=0; i<n; ++i) {
    eta.row(i) = X.row(i) * beta;
  }
  return eta;
}


/*** R
n=20
X = cbind(rep(1,n), rnorm(n, mean=0, sd=4))
beta = c(1,1)
eta(beta, X)
*/
