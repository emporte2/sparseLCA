#ifndef distributions_mvn
#define distributions_mvn

# include <RcppArmadillo.h>
# include <math.h>
# include <Rmath.h>
# include <R.h>
# include <Rcpp.h>
# include <iostream>

// FUNCTION DECLARATIONS
// ---------------------
arma::vec log_ddirichlet(const arma::vec eta, const arma::vec alpha);

arma::mat sample_covariance(const arma::mat Ym);


const double log2pi = std::log(2.0 * M_PI);

arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov);


arma::vec dmvnrm_arma(arma::mat x,
                      arma::rowvec mean,
                      arma::mat sigma,
                      bool logd = false,
                      bool sumd=false);


arma::vec dmvnorm_arma_prec(arma::mat x, arma::rowvec mean, arma::mat precision, bool log = false);


arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);


double log_dwishart_sample(arma::cube Sigma_sample, int df, arma::mat Winv);

double log_ptheta(arma::mat Theta, arma::vec mu, double kappa,arma::cube Sigma_sample);


double log_dirichlet(arma::vec pvec, arma::vec alphavec);

arma::mat wishartArma_int(int nu, arma::mat V);


arma::mat iwishartArma_int(int nu, arma::mat S);
#endif
