#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

//' C++ (fast) implementation of correlation for large sample sets
// [[Rcpp::export]]
NumericMatrix fast_corr(NumericMatrix x) {
  int ntaxa = x.nrow(), nsample = x.ncol();
  NumericMatrix corr(nsample, nsample);
  for(int i = 0; i < nsample; i++) {
    for(int j = 0; j < nsample; j++) {
      if(i == j) {
        corr(i,j) = 1;
      } else if(j > i) {
        double dp1 = 0;
        double dp2 = 0;
        double dp3 = 0;
        for(int k = 0; k < ntaxa; k++) {
          dp1 += x(k,i)*x(k,j);
          dp2 += x(k,i)*x(k,i);
          dp3 += x(k,j)*x(k,j);
        }
        corr(i,j) = dp1/(sqrt(dp2)*sqrt(dp3));
        corr(j,i) = corr(i,j);
      }
    }
  }
  return(corr);
}