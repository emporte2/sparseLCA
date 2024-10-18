#ifndef mvn_likelihood
#define mvn_likelihood

# include <RcppArmadillo.h>
# include <math.h>
# include <Rmath.h>
# include <R.h>
# include <Rcpp.h>
# include <iostream>

// FUNCTION DECLARATIONS
// ---------------------
arma::mat log_conditional_likelihood_mvn_cp(const arma::mat Y,
                                            const int K,
                                            const arma::mat Thetam,
                                            const arma::cube Sigmaa);

arma::mat expected_log_likelihood_mvn_cp(const arma::mat Y,
                                         const arma::mat T,
                                         const int K,
                                         const arma::mat Theta,
                                         const arma::cube Sigma);

arma::vec expected_log_p_allocations(const arma::mat Tmat,const int K,const arma::vec etav);

arma::vec log_likelihood_mvn_cp(const arma::mat Y,
                                const int K,
                                const arma::mat Thetam,
                                const arma::cube Sigmaa,
                                const arma::vec eta);

#endif
