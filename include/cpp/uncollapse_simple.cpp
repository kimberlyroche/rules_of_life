#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/random/mersenne_twister.hpp>  // included in BH  
#include "MatDist.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::ArrayXXd;
using Eigen::Map;

// these should be passed as references

// [[Rcpp::export]]
List uncollapse_simple(const Eigen::Map<Eigen::VectorXd> eta,
                       const Eigen::Map<Eigen::MatrixXd> X, 
                       const Eigen::Map<Eigen::MatrixXd> Theta,
                       const Eigen::Map<Eigen::MatrixXd> Gamma, 
                       const Eigen::Map<Eigen::MatrixXd> Xi, 
                       const double upsilon){
  List out(2);
  out.names() = CharacterVector::create("Lambda", "Sigma");
  int Q = Gamma.rows();
  int D = Xi.rows() + 1;
  int N = X.cols();
  int iter = eta.size()/(N*(D-1)); // assumes result is an integer !!!
  double upsilonN = upsilon + N;
  const MatrixXd GammaInv(Gamma.lu().inverse());
  const MatrixXd GammaInvN(GammaInv + X*X.transpose());
  const MatrixXd GammaN(GammaInvN.lu().inverse());
  const MatrixXd LGammaN(GammaN.llt().matrixL());
  //const Map<const MatrixXd> Eta(NULL);
  const MatrixXd ThetaGammaInvGammaN(Theta*GammaInv*GammaN);
  const MatrixXd XTGammaN(X.transpose()*GammaN);
  
  // Storage for output
  MatrixXd LambdaDraw0((D-1)*Q, iter);
  MatrixXd SigmaDraw0((D-1)*(D-1), iter);
  
  // iterate over all draws of eta - we're forgoing parallelization rn
  // boost::random::mt19937 rng(1);
  // storage for computation
  MatrixXd LambdaN(D-1, Q);
  MatrixXd XiN(D-1, D-1);
  MatrixXd LSigmaDraw(D-1, D-1);
  MatrixXd ELambda(D-1, Q);
  MatrixXd EEta(D-1, N);
  for (int i=0; i < iter; i++){
    //R_CheckUserInterrupt();
    const Map<const MatrixXd> Eta(&eta(i*N*(D-1)),D-1, N);
    LambdaN.noalias() = Eta*XTGammaN+ThetaGammaInvGammaN;
    ELambda = LambdaN-Theta;
    EEta.noalias() = Eta-LambdaN*X;
    XiN.noalias() = Xi+ EEta*EEta.transpose() + ELambda*GammaInv*ELambda.transpose();
    
    // Draw Random Component
    LSigmaDraw = rInvWishRevCholesky(upsilonN, XiN);
    // Note: Below is valid even though LSigmaDraw is reverse cholesky factor
    Eigen::Ref<VectorXd> LambdaDraw_tmp = LambdaDraw0.col(i);
    Eigen::Map<MatrixXd> LambdaDraw(LambdaDraw_tmp.data(), D-1, Q);
    LambdaDraw = rMatNormalCholesky(LambdaN, LSigmaDraw, LGammaN.matrix());
    
    Eigen::Ref<VectorXd> SigmaDraw_tmp = SigmaDraw0.col(i);
    Eigen::Map<MatrixXd> SigmaDraw_tosquare(SigmaDraw_tmp.data(), D-1, D-1);
    SigmaDraw_tosquare.noalias() = LSigmaDraw*LSigmaDraw.transpose();
  }
  
  IntegerVector dLambda = IntegerVector::create(D-1, Q, iter);
  IntegerVector dSigma = IntegerVector::create(D-1, D-1, iter);
  NumericVector nvLambda = wrap(LambdaDraw0);
  NumericVector nvSigma = wrap(SigmaDraw0);
  nvLambda.attr("dim") = dLambda;
  nvSigma.attr("dim") = dSigma;
  out[0] = nvLambda;
  out[1] = nvSigma;
  return out;
}












