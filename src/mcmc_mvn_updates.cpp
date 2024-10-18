# include <RcppArmadillo.h>
# include <math.h>
# include <Rmath.h>
# include <R.h>
# include <iostream>
# include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
# include "distributions_mvn.h"
# include "helpers.h"
# include "mvn_likelihood.h"


// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// arma::field<arma::cube> test_field(int K, int p, int nit)
// {
//   arma::field<arma::cube> out(3,1);
//   arma::cube eta(nit,K,1); eta.fill(0);
//   arma::cube beta(nit,K,p); beta.fill(1);
//   arma::cube ssq(nit,1,1); ssq.fill(-999);
//   out(0,0) = eta; out(1,0)=beta; out(2,0)=ssq;
//   return(out);
// }

// [[Rcpp::depends(RcppArmadillo)]]
arma::vec update_betaj_mcmc(arma::mat yj, const arma::mat Xj,
                           const int K, const double ssq,
                           const arma::vec prior_mean,
                           const arma::mat V_inv)
{
  int p = Xj.n_cols; int nj=Xj.n_rows;
  arma::vec out(p);
  arma::mat A(p,p); A.fill(0);
  arma::mat Ainv(p,p);
  arma::mat S(nj,nj); S.fill(0);
  arma::vec mu(p); mu.fill(0);
  S.diag() += 1.0/ssq;

  A=V_inv + Xj.t()*S*Xj;
  Ainv = inv(A);
  mu = Ainv*(Xj.t()*S*yj + V_inv*prior_mean);
  out = mvrnormArma(1, mu, Ainv).row(0).t();
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
arma::mat update_beta_mcmc(arma::vec Y, const arma::mat X,
                           const int K,
                           const arma::vec Sv, const double ssq,
                           const arma::vec prior_mean,
                           const arma::mat V_inv)
{
  int p = X.n_cols;
  arma::vec nvec(K); nvec = nvec_count(Sv,K);
  arma::mat out(K,p);
  for( int k=0; k<K; k++)
  {
    arma::uvec ind = find(Sv==k+1);
    arma::mat Xj(nvec(k),p); Xj=X.rows(ind);
    arma::vec yj(nvec(k)); yj=Y(ind);
    out.row(k) = update_betaj_mcmc(yj,Xj,K,ssq,prior_mean,V_inv).t();
  }
  return(out);
}



// [[Rcpp::depends(RcppArmadillo)]]
arma::mat update_Sigmaj_mcmc(const arma::mat YP, const arma::vec mu0, const double kappa0,
                                   const int nu0, const arma::mat W_inv)
{
  int n = YP.n_rows;
  int p = YP.n_cols;
  arma::mat YP_mean(p,1);  YP_mean.fill(0.0);
  arma::mat mu0M(p,1);  mu0M.fill(-999.-999);

  int nu_post = 0; nu_post = -999;
  arma::mat W_post(p,p); W_post.fill(0.0);

  if( n ==0 )
  {
    nu_post = nu0;
    W_post = inv_sympd(W_inv);
  }

  if( n >0 )
  {
    for( int pp = 0; pp < p; pp++)
    {
      for( int i = 0; i < n; i ++)
      {
        YP_mean(pp,0) += YP(i,pp);
        mu0M(pp,0) = mu0[pp];
      }

      YP_mean(pp,0) = YP_mean(pp,0)/n;

    }
    nu_post = n + nu0;
    arma::mat Cm(p,p); Cm = sample_covariance(YP)*n;
    arma::mat Dm(p,p); Dm = (kappa0*n/(kappa0 + n))*(YP_mean - mu0M)*trans(YP_mean - mu0M);
    W_post = inv(W_inv + Cm + Dm);
  }

  arma::mat Sig_inv(p,p);
  Sig_inv = wishartArma_int(nu_post, W_post);


  return(inv_sympd(Sig_inv));

}

// [[Rcpp::depends(RcppArmadillo)]]
arma::mat update_Sigmaj_ind_mcmc(const arma::mat YP, const arma::vec theta,
                             const int df0, const arma::mat prior_Scale)
{
  int n = YP.n_rows;
  int p = YP.n_cols;
  //arma::mat YP_mean(p,1);  YP_mean.fill(0.0);
  //arma::mat mu0M(p,1);  mu0M.fill(-999.-999);

  int df_post = 0; df_post = -999;
  arma::mat D_post(p,p); D_post.fill(0.0);
  arma::mat Dm(p,p); Dm.fill(0.0);
  arma::mat D_postinv(p,p); D_postinv.fill(0.0);

  if( n ==0 )
  {
    df_post = df0;
    D_post = prior_Scale;
  }

  if( n >0 )
  {

    for( int i = 0; i < n; i ++)
    {
      Dm += ((trans(YP.row(i)) - theta)*trans(trans(YP.row(i)) - theta));
    }
    df_post = n + df0;
    D_post = (prior_Scale +Dm);

  }
  D_postinv = inv_sympd(D_post);
  arma::mat Sig_inv(p,p);
  Sig_inv = wishartArma_int(df_post, D_postinv);


  //return(inv_sympd(Sig_inv));
  return(Sig_inv);

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube update_Sigma_ind_mcmc(const arma::mat YY,const arma::vec s, const int K, const arma::mat theta,
                             const double df0, const arma::mat prior_Scale)
{
  //int n = YY.n_rows;
  int p = YY.n_cols;
  arma::cube out(p,p,K);
  arma::vec nvec(K); nvec = nvec_count(s,K);

  for( int k=0; k<K; k++)
  {
    arma::uvec ind = find(s==k+1);
    arma::mat YYk(nvec(k),p); YYk=YY.rows(ind);
    arma::vec theta_temp(p); theta_temp=trans(theta.row(k));

    out.slice(k) = update_Sigmaj_ind_mcmc(YYk, theta_temp, df0, prior_Scale);
  }


  return(out);

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube update_Sigma_mcmc(const arma::mat YY,const arma::vec s, const int K, const arma::vec mu0, const double kappa0,
                                   const int df, const arma::mat Winv)
{
  //int n = YY.n_rows;
  int p = YY.n_cols;
  arma::cube out(p,p,K);
  arma::vec nvec(K); nvec = nvec_count(s,K);

  for( int k=0; k<K; k++)
  {
    arma::uvec ind = find(s==k+1);
    arma::mat YYk(nvec(k),p); YYk=YY.rows(ind);
    out.slice(k) = update_Sigmaj_mcmc(YYk, mu0, kappa0, df, Winv);
  }


  return(out);

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec update_eta(const int K, const arma::vec Sv, const double alphap)
{
  arma::vec apost(K); apost.fill(alphap);
  arma::vec nvec(K); nvec = nvec_count(Sv,K);
  apost +=nvec;
  arma::vec out(K);
  for(int k=0; k<K; k++)
  {
    out(k)= R::rgamma(apost(k),1);
  }
  return(out/sum(out));
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double rg_check(const double ap,
                  const double bp)
{
  double out =0;
  out=1.0/R::rgamma(ap, bp);
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_class_probs_mvn_cp(const arma::mat Y,
                                    const int KK,
                                    const arma::mat current_theta,
                                    const arma::cube current_Sigma,
                                    const arma::vec etav)
{
  int N = Y.n_rows; //int p=Y.n_cols;
  arma::mat fm1(N,KK); fm1.zeros();
  arma::mat fm2(N,KK); fm2.zeros();
  arma::mat oneM(N,1); oneM.ones();
  fm1 = log_conditional_likelihood_mvn_cp(Y, KK,current_theta,current_Sigma)  + oneM*log(etav.t());

  for( int i=0; i < N; i++)
  {
    double sumt=0;
    sumt = logsumexp_vec(fm1.row(i).t());
    fm2.row(i) = exp(fm1.row(i) - sumt);
  }
  return(fm2);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec sample_S_mvn_mcmc(const arma::mat Y,
                        const int KK,
                        const arma::mat current_Theta,
                        const arma::cube current_Sigma,
                        const arma::vec etav)
{
  int N = Y.n_rows;
  arma::mat Tm(N,KK); Tm.fill(0);
  Tm = update_class_probs_mvn_cp(Y, KK,current_Theta,current_Sigma,etav);
  arma::vec s_out(N);
  arma::uvec classes(KK);
  for( int k=0; k<KK; k++)
  {
    classes(k) = k+1;
  }

  for(int n=0; n<N; n++)
  {
    arma::vec pr(KK);
    pr = Tm.row(n).t();
    s_out(n) = RcppArmadillo::sample(classes,1,false, pr)(0);
  }

  return(s_out);
}


// [[Rcpp::depends(RcppArmadillo)]]
arma::vec update_thetaj_mcmc(const arma::mat YP, const arma::mat Sigmaj, const arma::vec mu0, const double kappa0)
{
  int n = YP.n_rows;
  int p = YP.n_cols;
  arma::mat YP_mean(p,1);  YP_mean.fill(0.0);

  arma::vec mu_post(p);  mu_post.fill(-999);
  double kappa_post = 0; kappa_post= -999;

  if( n ==0 )
  {
    mu_post = mu0;
    kappa_post = kappa0;
  }

  if( n >0 )
  {
    for( int pp = 0; pp < p; pp++)
    {
      for( int i = 0; i < n; i ++)
      {
        YP_mean(pp,0) += YP(i,pp);
      }

      YP_mean(pp,0) = YP_mean(pp,0)/n;
      mu_post(pp) = (kappa0*mu0(pp) + n*YP_mean(pp,0)) /( kappa0 + n);

    }
    kappa_post = n + kappa0;
  }
  arma::vec theta_samp(p); theta_samp.fill(-990);
  double invkappa = 1.0/kappa_post;
  theta_samp = (mvrnormArma(1, mu_post, invkappa*Sigmaj).row(0)).t();

  return(theta_samp);

}

// [[Rcpp::depends(RcppArmadillo)]]
arma::vec update_thetaj_ind_mcmc(const arma::mat YP, const arma::mat Sigmaj_inv, const arma::vec mu0, const arma::mat B_inv)
{
  int n = YP.n_rows;
  int p = YP.n_cols;
  arma::mat YP_mean(p,1);  YP_mean.fill(0.0);

  arma::vec mu_post(p);  mu_post.fill(-999);
  arma::mat V_post(p,p);  V_post.fill(-999);
  if( n ==0 )
  {
    mu_post = mu0;
    V_post = arma::inv(B_inv);
  }

  if( n >0 )
  {
    for( int pp = 0; pp < p; pp++)
    {
      for( int i = 0; i < n; i ++)
      {
        YP_mean(pp,0) += YP(i,pp);
      }

      YP_mean(pp,0) = YP_mean(pp,0)/n;
      //mu_post(pp) = (B_inv*mu0(pp) + n*YP_mean(pp,0)) /( kappa0 + n);

    }
    V_post = arma::inv(B_inv + n*Sigmaj_inv);
    mu_post = (V_post*(B_inv*mu0 + n*Sigmaj_inv*YP_mean));
  }
  arma::vec theta_samp(p); theta_samp.fill(-990);
  //double invkappa = 1.0/kappa_post;
  theta_samp = (mvrnormArma(1, mu_post, V_post).row(0)).t();

  return(theta_samp);

}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_theta_ind_mcmc(const arma::mat YP, const arma::vec s, const int K, const arma::cube Sigma, const arma::vec mu0, const arma::mat B_inv)
{
  int p = YP.n_cols;
  arma::mat out(K,p);
  arma::vec nvec(K); nvec = nvec_count(s,K);


  for( int k=0; k<K; k++)
  {
    arma::uvec ind = find(s==k+1);
    arma::mat YYk(nvec(k),p); YYk=YP.rows(ind);
    out.row(k) = update_thetaj_ind_mcmc(YYk, Sigma.slice(k),mu0,B_inv).t();
  }

  return(out);

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_theta_mcmc(const arma::mat YP, const arma::vec s, const int K, const arma::cube Sigma, const arma::vec mu0, const double kappa0)
{
  int p = YP.n_cols;
  arma::mat out(K,p);
  arma::vec nvec(K); nvec = nvec_count(s,K);


  for( int k=0; k<K; k++)
  {
    arma::uvec ind = find(s==k+1);
    arma::mat YYk(nvec(k),p); YYk=YP.rows(ind);
    out.row(k) = update_thetaj_mcmc(YYk, Sigma.slice(k),mu0,kappa0).t();
  }

  return(out);

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double log_e0_gamma_ratio(double e0, double e0_old, arma::vec gamma,  const int Kbig)
{
  // log p(gamma | e0); eq 10 SFM 2016
  double out=0;
  out += lgamma(Kbig*e0) -Kbig*lgamma(e0) +Kbig*lgamma(e0_old)-lgamma(Kbig*e0_old);
  double lnum = 0; double ldenom=0;
  for(int k=0; k<Kbig; k++)
  {
    lnum += (e0-1)*log(gamma(k));
    ldenom += (e0_old-1)*log(gamma(k));
  }
  return(out+lnum-ldenom);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double log_eqd1_ratio(double e0, double e0_old, arma::vec cvec,  const int Kbig)
{
  double out=0;
  int N=cvec.n_elem;
  arma::vec nvec = count_classes(cvec, Kbig);
  arma::uvec(nonzero) = find(nvec>0);
  //int kplus = nonzero.n_elem;
  //out += gamma_small(Kbig*e0)/tgamma(N+Kbig*e0)*tgamma(N+Kbig*e0_old)/gamma_small(Kbig*e0_old);
  out += lgamma(Kbig*e0) -lgamma(N+Kbig*e0) +lgamma(N+Kbig*e0_old)-lgamma(Kbig*e0_old);
  double lnum = 0; double ldenom=0;
  for(int k=0; k<Kbig; k++)
  {
    lnum += lgamma(nvec(k)+e0)-lgamma(e0);
    ldenom += lgamma(nvec(k)+e0_old)-lgamma(e0_old);
  }
  return(out+lnum-ldenom);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double log_prior_ratio(double e0, double e0_old, double ae, double be)
{
  // This is with be=scale parameterization
  double out=0;
  //out+= exp( (e0_old - e0)/be)*pow((e0/e0_old), ae-1);
  out += (e0_old - e0)/be + (ae-1)*log(e0/e0_old);
  return(out);
}
