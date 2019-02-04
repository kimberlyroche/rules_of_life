#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXXd;
using Eigen::VectorXd;

// [[Rcpp::export]]
double logd_matrixt(double s1, double s2, double s3, double s4, double s5, double s6, double s7,
                    MatrixXd ilr_data, MatrixXd week_kernel, MatrixXd season_kernel, MatrixXd group_kernel,
                    MatrixXd age_kernel, MatrixXd indiv_kernel, MatrixXd plate_kernel, MatrixXd conc_kernel) {
  double N = ilr_data.cols();
  double P = ilr_data.rows();
  double upsilon = P + 10;
  MatrixXd K = MatrixXd::Identity(P,P);
  // constrainst imposed in optim()
  //  MatrixXd A = exp(s1)*week_kernel + exp(s2)*season_kernel + exp(s3)*group_kernel +
  //               exp(s4)*age_kernel + exp(s5)*indiv_kernel + exp(s6)*plate_kernel; + exp(s7)*conc_kernel;
  MatrixXd A = s1*week_kernel + s2*season_kernel + s3*group_kernel +
               s4*age_kernel + s5*indiv_kernel + s6*plate_kernel + s7*conc_kernel;
  // returns lower triangular L of A
  MatrixXd cholA = A.llt().matrixL();
  VectorXd logDiagCholA = cholA.diagonal().array().log();
  double logDetA = 0.0;
  for(int i = 0; i < logDiagCholA.size(); i++) {
    logDetA += logDiagCholA[i];
  }
  logDetA *= 2;
  MatrixXd K_inv = K.inverse();
  MatrixXd cholAInv = cholA.inverse();
  MatrixXd A_inv = (cholAInv.transpose())*cholAInv;
  MatrixXd S = MatrixXd::Identity(P,P) + (1/upsilon)*K_inv*ilr_data*A_inv*(ilr_data.transpose());
  MatrixXd cholS = S.llt().matrixL();
  VectorXd logDiagCholS = cholS.diagonal().array().log();
  double logDetS = 0.0;
  for(int i = 0; i < logDiagCholS.size(); i++) {
    logDetS += logDiagCholS[i];
  }
  logDetS *= 2;
  double d = (P/2)*logDetA + ((upsilon+N+P-1)/2)*logDetS;
  return(d);
}
