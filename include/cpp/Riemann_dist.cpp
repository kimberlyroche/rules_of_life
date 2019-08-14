#include <Rcpp.h>
#include <RcppEigen.h>
#include <omp.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

using Eigen::MatrixXd;
using Eigen::VectorXd;

//' Calculate distances between a matrix of positive-definite covariance matrices concatenated
//' horizontally (i.e. D rows x DN columns where D=features, N=samples)
// [[Rcpp::export]]
Eigen::MatrixXd Riemann_dist_samples(Eigen::MatrixXd samples, int n_indiv, int n_samples_per) {
  // Eigen parallelization doesn't appear to fight with this
  int P = samples.rows();
  MatrixXd distances((n_indiv*n_samples_per), (n_indiv*n_samples_per));
  #pragma omp parallel shared(samples, P, distances) num_threads(omp_get_max_threads())
  {
  #pragma omp for collapse(2)
  for(int i=0; i < n_indiv*n_samples_per; i++) {
    for (int j=0; j < n_indiv*n_samples_per; j++) {
      if(j >= i) {
        // get samples A and B
        int a_idx = i*P;
        MatrixXd A = samples.middleCols(a_idx, P);
        int b_idx = j*P;
        MatrixXd B = samples.middleCols(b_idx, P);
        // calculate distance
        Eigen::LLT<MatrixXd> LLT_A;
        LLT_A.compute(A);
        MatrixXd L = LLT_A.matrixL();
        MatrixXd invL = L.inverse();
        MatrixXd temp = invL*B*(invL.transpose());
        Eigen::EigenSolver<MatrixXd> es(temp);
        double d = sqrt(es.eigenvalues().array().log().abs2().sum());
        distances(i,j) = d;
      }
    }
  }
  }
  // copy the lower triangular from the upper
  // these distances are demonstrably symmetric
  for(int i=0; i < n_indiv*n_samples_per; i++) {
    for (int j=0; j < n_indiv*n_samples_per; j++) {
      if(j < i) {
        distances(i,j) = distances(j,i);
      }
    }
  }
  return(distances);
}

//' Calculate distances between positive-definite covariance matrices A and B
// [[Rcpp::export]]
double Riemann_dist_pair(Eigen::MatrixXd A, Eigen::MatrixXd B) {
  Eigen::LLT<MatrixXd> LLT_A;
  LLT_A.compute(A);
  MatrixXd L = LLT_A.matrixL();
  MatrixXd invL = L.inverse();
  MatrixXd temp = invL*B*(invL.transpose());
  Eigen::EigenSolver<MatrixXd> es(temp);
  double d = sqrt(es.eigenvalues().array().log().abs2().sum());
  // Frobenius difference alternative
  // MatrixXd C = A - B;
  // double d = sqrt(C.array().abs2().sum());
  return(d);
}

