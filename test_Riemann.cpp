#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

using Eigen::MatrixXd;
using Eigen::VectorXd;

//' Calculate distances between positive-definite covariance matrices
// [[Rcpp::export]]
double mat_dist(Eigen::MatrixXd A, Eigen::MatrixXd B, bool use_Riemann) {
  double d = -1;
  if(use_Riemann) {
    Eigen::LLT<MatrixXd> LLT_A;
    LLT_A.compute(A);
    MatrixXd L = LLT_A.matrixL();
    MatrixXd invL = L.inverse();
    MatrixXd temp = invL*B*(invL.transpose());
    Eigen::EigenSolver<MatrixXd> es(temp);
    d = sqrt(es.eigenvalues().array().log().abs2().sum());
  } else {
    MatrixXd C = A - B;
    d = sqrt(C.array().abs2().sum());
  }
  return(d);
}

