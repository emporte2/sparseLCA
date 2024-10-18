# include <RcppArmadillo.h>
# include <math.h>
# include <Rmath.h>
# include <R.h>
# include <iostream>
using namespace Rcpp;
//# include "distributions_mixreg.cpp"
# include "distributions_mvn.h"
# include "helpers.h"



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat log_conditional_likelihood_mvn_cp(const arma::mat Y,
                     const int K,
                     const arma::mat Thetam,
                     const arma::cube Sigmaa)
{
  //int p = betam.n_cols;
  int n = Y.n_rows;


    arma::mat llikM(n,K); llikM.fill(0);
    for( int k=0; k < K; k++)
    {
      arma::vec out(n); out.fill(0);
      out = dmvnrm_arma(Y, Thetam.row(k), Sigmaa.slice(k), true,false);
      llikM.col(k) = out;
    }

  return(llikM);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat expected_log_likelihood_mvn_cp(const arma::mat Y,
                                   const arma::mat T,
                                   const int K,
                                   const arma::mat Theta,
                                   const arma::cube Sigma)
{
  int n=Y.n_rows;
  arma::mat L1(n,K); L1.fill(0);
  L1 = log_conditional_likelihood_mvn_cp(Y, K,Theta,Sigma);
  arma::vec out(1); out.fill(0);
  for(int i=0; i<n; i++)
  {
    arma::vec temp(K); temp=(L1.row(i)%T.row(i)).t();
    out(0) += sum(temp);
  }
  return(out);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec expected_log_p_allocations(const arma::mat Tmat,const int K,const arma::vec etav)
{
  arma::vec out(1);
  for(int j=0; j<K; j++)
  {
    out += sum(Tmat.col(j))*log(etav(j));
  }
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec log_likelihood_mvn_cp(const arma::mat Y,
                                            const int K,
                                            const arma::mat Thetam,
                                            const arma::cube Sigmaa,
                                            const arma::vec eta)
{
  //int p = betam.n_cols;
  int n = Y.n_rows;
  arma::mat condll(n,K);
  condll= log_conditional_likelihood_mvn_cp(Y,K,Thetam,Sigmaa);
  arma::vec out(n);

  arma::mat llikM(n,K); llikM.fill(0);
  for( int i=0; i < n; i++)
  {
    arma::vec temp(K); temp.fill(0);
    temp = condll.row(i).t()+log(eta);
    out(i) = logsumexp_vec(temp);
  }

  return(out);
}
/*
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec anchor_likelihood_perms(const arma::vec Y,
                                  const arma::mat X,
                                  const int K,
                                  const arma::mat betam,
                                  const double sigsq,
                                  const arma::mat AnchorsTall,
                                  const arma::mat Perms)
{
  int Q=Perms.n_rows;
  int n_anchors = AnchorsTall.n_rows;
  arma::vec out(Q); out.fill(0);
  //arma::mat out2(Q,n_anchors); out2.fill(0);
  arma::mat LikeM(n_anchors,K);
  LikeM = log_conditional_likelihood_mixreg_cp(Y,X,K,betam,sigsq);

  for(int q=0; q<Q; q++)
  {
    for(int i=0; i<n_anchors; i++)
    {
      out(q) += LikeM(i,Perms(q,AnchorsTall(i,1)-1)-1);
      //out2(q,i) =LikeM(i,Perms(q,AnchorsTall(i,1)-1)-1);
    }
  }
  return(out);
}
 */

/*
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec anchor_perms_pvec(const arma::vec Y,
                                  const arma::mat X,
                                  const int K,
                                  const arma::mat betam,
                                  const double sigsq,
                                  const arma::mat AnchorsTall,
                                  const arma::mat Perms)
{
  int Q=Perms.n_rows;
  arma::vec out(Q); out.fill(0);
  out=anchor_likelihood_perms(Y,X,K,betam,sigsq,AnchorsTall,Perms);
  return(exp(out- logsumexp_vec(out)));
}
 */
