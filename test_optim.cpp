#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector test_optim() {
  NumericVector testvec(10);
  return(testvec);
}

