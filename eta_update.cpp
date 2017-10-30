// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


/*** R
timesTwo(42)
*/
