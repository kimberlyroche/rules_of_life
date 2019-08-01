#include <Rcpp.h>
#include <RcppEigen.h>
#include <omp.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::export]]
Eigen::MatrixXd mat_dist(Eigen::MatrixXd samples, int n_indiv, int n_samples_per) {
  Rcout << "OpenMP max threads: " << omp_get_max_threads() << std::endl;
  // Eigen parallelization might be a problem here, we'll see
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

