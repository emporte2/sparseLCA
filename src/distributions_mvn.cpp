# include <RcppArmadillo.h>
# include <math.h>
# include <Rmath.h>
# include <R.h>
//# include <RcppDist.h>

/* mvnMix.cpp: functions need for a MVN MFM model  */
# include "helpers.h"




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec log_ddirichlet(const arma::vec eta, const arma::vec alpha)
{
  arma::vec OutM(1); OutM(0) = 0;
  double Bnum = 1;
  double Bconst = 1;
  double alpha_sum = 0;
  int k = eta.n_elem;
  for( int j = 0; j < k ; j ++)
  {
    OutM += alpha(j)*log(eta(j));
    Bnum *= tgamma(alpha(j));
    alpha_sum += alpha(j);
  }
  Bconst = Bnum/tgamma(alpha_sum);
  OutM += -log(Bconst);
  return(OutM);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat sample_covariance(const arma::mat Ym)
{
  // Scaled by N, not N-1!!
  int n = Ym.n_rows;
  int p = Ym.n_cols;
  arma::mat Ymi(p,1);  Ymi.fill(-999.-999);
  arma::mat Ymean(p,1);  Ymean.fill(0.0);

  for( int pp = 0; pp < p ; pp ++)
  {
    for( int i=0; i < n; i ++)
    {
    Ymean(pp,0) += Ym(i,pp);
    }
  }
  Ymean = Ymean/n;

  arma::mat OutM(p,p); OutM.fill(0.0);

  for( int i = 0; i < n; i ++)
  {
    for( int pp = 0; pp < p; pp ++)
    {
      Ymi(pp,0) = Ym(i,pp) - Ymean(pp,0);
    }

    OutM += Ymi*Ymi.t();
  }
  OutM *= 1.0/(n);
  return(OutM);

}


/* Reference for the folllowing MVN density functions
  https://gallery.rcpp.org/articles/dmvnorm_arma/
 */
const double log2pi = std::log(2.0 * M_PI);


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov) {
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec dmvnrm_arma(arma::mat x,
                      arma::rowvec mean,
                      arma::mat sigma,
                      bool logd = false,
                      bool sumd=false) {
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
  }

  if ( (logd == false) & !sumd) {
    out = exp(out);
  }
  if(sumd)
  {
    out = sum(out);
  }
  return(out);
}
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec dmvnorm_arma_prec(arma::mat x, arma::rowvec mean, arma::mat precision, bool log = false) {
  int p = mean.n_elem;
  arma::mat sigma(p,p);
  sigma = inv_sympd(precision);
  arma::vec distval = Mahalanobis(x,  mean, sigma);
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;

  if (log) {
    return(logretval);
  } else {
    return(exp(logretval));
  }
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double log_dwishart_sample(arma::cube Sigma_sample, int df, arma::mat Winv)
{
  //unnormalized!
  int n = Sigma_sample.n_slices;
  int p = Winv.n_cols;
  double out =0;
  for( int i=0; i<n; i++)
  {
    arma::mat Temp(p,p); Temp = Sigma_sample.slice(i);
  out+= ((n-p-1.0)/2.0)*log_det_sympd(Temp) - trace(Winv*Temp)/2.0;
  }

  return(out);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double log_ptheta(arma::mat Theta, arma::vec mu, double kappa,arma::cube Sigma_sample)
{
  int p=Theta.n_cols;
  int K=Theta.n_rows;
  double kappainv = 1.0/kappa;
  double out = 0;
  for( int k=0; k<K; k++)
  {
    arma::mat Sigmatemp(p,p);
    Sigmatemp = kappainv*Sigma_sample.slice(k);
    out += dmvnrm_arma(Theta.row(k), mu.t(),Sigmatemp,true,true)(0);
  }
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double log_dirichlet(arma::vec pvec, arma::vec alphavec)
{
  double out; out=lgamma(sum(alphavec));
  int K = alphavec.n_elem;
  for(int k=0; k<K; k++)
  {
    out += (alphavec(k)-1)*log(pvec(k))-lgamma(alphavec(k));
  }
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat wishartArma_int(int nu, arma::mat V) {

  int p = V.n_rows;
  arma::vec zero_vec(p); zero_vec = arma::zeros(p);
  arma::mat X(nu,p);  X.fill(-999.-999);

  X = mvrnormArma(nu, zero_vec, V);

  return X.t()*X;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat iwishartArma_int(int nu, arma::mat S) {

  int p = S.n_rows;
  arma::vec zero_vec(p); zero_vec = arma::zeros(p);
  arma::mat X(nu,p);  X.fill(-999.-999);
  arma::mat X2_temp(p,p);

  X = mvrnormArma(nu, zero_vec, inv(S));

  X2_temp = inv(X.t()*X);
  return((X2_temp));
}


