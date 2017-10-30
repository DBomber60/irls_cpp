// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List v(mat XtWX, mat XtWz, arma::mat X, colvec W, colvec Z) {
  int n = size(X, 0);
  for(int i=0; i<n; ++i) {
    double w = as_scalar(W.row(i));
    double z = as_scalar(Z.row(i));
    XtWX += (X.row(i)*w).t()*(X.row(i));
    XtWz += (X.row(i)*w*z).t();
  }
  return List::create(Named("XtWX")=XtWX, Named("XtWz")=XtWz);
  }


/*** R
n=20
X = cbind(rep(1,n), rnorm(n, mean=0, sd=4))
XtWX = array(0, dim = c(ncol(X), ncol(X)))
XtWz = array(0, dim = c(ncol(X), 1))
W = Z = rep(1,n)
ss = v(XtWX, XtWz, X, W, Z)

solve(ss$XtWX) %*% ss$XtWz
solve(t(X) %*% diag(W) %*% X) %*% (t(X) %*% diag(W) %*% Z)

*/

