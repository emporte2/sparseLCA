// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// log_ddirichlet
arma::vec log_ddirichlet(const arma::vec eta, const arma::vec alpha);
RcppExport SEXP _sparselca_log_ddirichlet(SEXP etaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(log_ddirichlet(eta, alpha));
    return rcpp_result_gen;
END_RCPP
}
// sample_covariance
arma::mat sample_covariance(const arma::mat Ym);
RcppExport SEXP _sparselca_sample_covariance(SEXP YmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Ym(YmSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_covariance(Ym));
    return rcpp_result_gen;
END_RCPP
}
// Mahalanobis
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov);
RcppExport SEXP _sparselca_Mahalanobis(SEXP xSEXP, SEXP centerSEXP, SEXP covSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type center(centerSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov(covSEXP);
    rcpp_result_gen = Rcpp::wrap(Mahalanobis(x, center, cov));
    return rcpp_result_gen;
END_RCPP
}
// dmvnrm_arma
arma::vec dmvnrm_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool logd, bool sumd);
RcppExport SEXP _sparselca_dmvnrm_arma(SEXP xSEXP, SEXP meanSEXP, SEXP sigmaSEXP, SEXP logdSEXP, SEXP sumdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type logd(logdSEXP);
    Rcpp::traits::input_parameter< bool >::type sumd(sumdSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnrm_arma(x, mean, sigma, logd, sumd));
    return rcpp_result_gen;
END_RCPP
}
// dmvnorm_arma_prec
arma::vec dmvnorm_arma_prec(arma::mat x, arma::rowvec mean, arma::mat precision, bool log);
RcppExport SEXP _sparselca_dmvnorm_arma_prec(SEXP xSEXP, SEXP meanSEXP, SEXP precisionSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type precision(precisionSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnorm_arma_prec(x, mean, precision, log));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
RcppExport SEXP _sparselca_mvrnormArma(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma(n, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// log_dwishart_sample
double log_dwishart_sample(arma::cube Sigma_sample, int df, arma::mat Winv);
RcppExport SEXP _sparselca_log_dwishart_sample(SEXP Sigma_sampleSEXP, SEXP dfSEXP, SEXP WinvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type Sigma_sample(Sigma_sampleSEXP);
    Rcpp::traits::input_parameter< int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Winv(WinvSEXP);
    rcpp_result_gen = Rcpp::wrap(log_dwishart_sample(Sigma_sample, df, Winv));
    return rcpp_result_gen;
END_RCPP
}
// log_ptheta
double log_ptheta(arma::mat Theta, arma::vec mu, double kappa, arma::cube Sigma_sample);
RcppExport SEXP _sparselca_log_ptheta(SEXP ThetaSEXP, SEXP muSEXP, SEXP kappaSEXP, SEXP Sigma_sampleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Sigma_sample(Sigma_sampleSEXP);
    rcpp_result_gen = Rcpp::wrap(log_ptheta(Theta, mu, kappa, Sigma_sample));
    return rcpp_result_gen;
END_RCPP
}
// log_dirichlet
double log_dirichlet(arma::vec pvec, arma::vec alphavec);
RcppExport SEXP _sparselca_log_dirichlet(SEXP pvecSEXP, SEXP alphavecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pvec(pvecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alphavec(alphavecSEXP);
    rcpp_result_gen = Rcpp::wrap(log_dirichlet(pvec, alphavec));
    return rcpp_result_gen;
END_RCPP
}
// wishartArma_int
arma::mat wishartArma_int(int nu, arma::mat V);
RcppExport SEXP _sparselca_wishartArma_int(SEXP nuSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(wishartArma_int(nu, V));
    return rcpp_result_gen;
END_RCPP
}
// iwishartArma_int
arma::mat iwishartArma_int(int nu, arma::mat S);
RcppExport SEXP _sparselca_iwishartArma_int(SEXP nuSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(iwishartArma_int(nu, S));
    return rcpp_result_gen;
END_RCPP
}
// logsumexp_vec
double logsumexp_vec(arma::vec avec);
RcppExport SEXP _sparselca_logsumexp_vec(SEXP avecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type avec(avecSEXP);
    rcpp_result_gen = Rcpp::wrap(logsumexp_vec(avec));
    return rcpp_result_gen;
END_RCPP
}
// logsumexp
double logsumexp(const double a, const double b);
RcppExport SEXP _sparselca_logsumexp(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(logsumexp(a, b));
    return rcpp_result_gen;
END_RCPP
}
// normalize_vec
arma::mat normalize_vec(arma::vec v);
RcppExport SEXP _sparselca_normalize_vec(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(normalize_vec(v));
    return rcpp_result_gen;
END_RCPP
}
// normalize_logvec
arma::mat normalize_logvec(arma::vec lv);
RcppExport SEXP _sparselca_normalize_logvec(SEXP lvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type lv(lvSEXP);
    rcpp_result_gen = Rcpp::wrap(normalize_logvec(lv));
    return rcpp_result_gen;
END_RCPP
}
// nvec_count
arma::vec nvec_count(arma::vec svec, int k);
RcppExport SEXP _sparselca_nvec_count(SEXP svecSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type svec(svecSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(nvec_count(svec, k));
    return rcpp_result_gen;
END_RCPP
}
// count_if
int count_if(arma::vec c_vec, const int cval);
RcppExport SEXP _sparselca_count_if(SEXP c_vecSEXP, SEXP cvalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c_vec(c_vecSEXP);
    Rcpp::traits::input_parameter< const int >::type cval(cvalSEXP);
    rcpp_result_gen = Rcpp::wrap(count_if(c_vec, cval));
    return rcpp_result_gen;
END_RCPP
}
// count_classes
arma::vec count_classes(arma::vec c_vec, const int Kbig);
RcppExport SEXP _sparselca_count_classes(SEXP c_vecSEXP, SEXP KbigSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type c_vec(c_vecSEXP);
    Rcpp::traits::input_parameter< const int >::type Kbig(KbigSEXP);
    rcpp_result_gen = Rcpp::wrap(count_classes(c_vec, Kbig));
    return rcpp_result_gen;
END_RCPP
}
// table_by_class
arma::cube table_by_class(arma::mat YY, arma::vec cvec, const int Kbig);
RcppExport SEXP _sparselca_table_by_class(SEXP YYSEXP, SEXP cvecSEXP, SEXP KbigSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type YY(YYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cvec(cvecSEXP);
    Rcpp::traits::input_parameter< const int >::type Kbig(KbigSEXP);
    rcpp_result_gen = Rcpp::wrap(table_by_class(YY, cvec, Kbig));
    return rcpp_result_gen;
END_RCPP
}
// update_Sigma_ind_mcmc
arma::cube update_Sigma_ind_mcmc(const arma::mat YY, const arma::vec s, const int K, const arma::mat theta, const double df0, const arma::mat prior_Scale);
RcppExport SEXP _sparselca_update_Sigma_ind_mcmc(SEXP YYSEXP, SEXP sSEXP, SEXP KSEXP, SEXP thetaSEXP, SEXP df0SEXP, SEXP prior_ScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type YY(YYSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type df0(df0SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type prior_Scale(prior_ScaleSEXP);
    rcpp_result_gen = Rcpp::wrap(update_Sigma_ind_mcmc(YY, s, K, theta, df0, prior_Scale));
    return rcpp_result_gen;
END_RCPP
}
// update_Sigma_mcmc
arma::cube update_Sigma_mcmc(const arma::mat YY, const arma::vec s, const int K, const arma::vec mu0, const double kappa0, const int df, const arma::mat Winv);
RcppExport SEXP _sparselca_update_Sigma_mcmc(SEXP YYSEXP, SEXP sSEXP, SEXP KSEXP, SEXP mu0SEXP, SEXP kappa0SEXP, SEXP dfSEXP, SEXP WinvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type YY(YYSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const double >::type kappa0(kappa0SEXP);
    Rcpp::traits::input_parameter< const int >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Winv(WinvSEXP);
    rcpp_result_gen = Rcpp::wrap(update_Sigma_mcmc(YY, s, K, mu0, kappa0, df, Winv));
    return rcpp_result_gen;
END_RCPP
}
// update_eta
arma::vec update_eta(const int K, const arma::vec Sv, const double alphap);
RcppExport SEXP _sparselca_update_eta(SEXP KSEXP, SEXP SvSEXP, SEXP alphapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type Sv(SvSEXP);
    Rcpp::traits::input_parameter< const double >::type alphap(alphapSEXP);
    rcpp_result_gen = Rcpp::wrap(update_eta(K, Sv, alphap));
    return rcpp_result_gen;
END_RCPP
}
// rg_check
double rg_check(const double ap, const double bp);
RcppExport SEXP _sparselca_rg_check(SEXP apSEXP, SEXP bpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type ap(apSEXP);
    Rcpp::traits::input_parameter< const double >::type bp(bpSEXP);
    rcpp_result_gen = Rcpp::wrap(rg_check(ap, bp));
    return rcpp_result_gen;
END_RCPP
}
// update_class_probs_mvn_cp
arma::mat update_class_probs_mvn_cp(const arma::mat Y, const int KK, const arma::mat current_theta, const arma::cube current_Sigma, const arma::vec etav);
RcppExport SEXP _sparselca_update_class_probs_mvn_cp(SEXP YSEXP, SEXP KKSEXP, SEXP current_thetaSEXP, SEXP current_SigmaSEXP, SEXP etavSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int >::type KK(KKSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type current_theta(current_thetaSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type current_Sigma(current_SigmaSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type etav(etavSEXP);
    rcpp_result_gen = Rcpp::wrap(update_class_probs_mvn_cp(Y, KK, current_theta, current_Sigma, etav));
    return rcpp_result_gen;
END_RCPP
}
// sample_S_mvn_mcmc
arma::vec sample_S_mvn_mcmc(const arma::mat Y, const int KK, const arma::mat current_Theta, const arma::cube current_Sigma, const arma::vec etav);
RcppExport SEXP _sparselca_sample_S_mvn_mcmc(SEXP YSEXP, SEXP KKSEXP, SEXP current_ThetaSEXP, SEXP current_SigmaSEXP, SEXP etavSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int >::type KK(KKSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type current_Theta(current_ThetaSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type current_Sigma(current_SigmaSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type etav(etavSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_S_mvn_mcmc(Y, KK, current_Theta, current_Sigma, etav));
    return rcpp_result_gen;
END_RCPP
}
// update_theta_ind_mcmc
arma::mat update_theta_ind_mcmc(const arma::mat YP, const arma::vec s, const int K, const arma::cube Sigma, const arma::vec mu0, const arma::mat B_inv);
RcppExport SEXP _sparselca_update_theta_ind_mcmc(SEXP YPSEXP, SEXP sSEXP, SEXP KSEXP, SEXP SigmaSEXP, SEXP mu0SEXP, SEXP B_invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type YP(YPSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type s(sSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type B_inv(B_invSEXP);
    rcpp_result_gen = Rcpp::wrap(update_theta_ind_mcmc(YP, s, K, Sigma, mu0, B_inv));
    return rcpp_result_gen;
END_RCPP
}
// log_e0_gamma_ratio
double log_e0_gamma_ratio(double e0, double e0_old, arma::vec gamma, const int Kbig);
RcppExport SEXP _sparselca_log_e0_gamma_ratio(SEXP e0SEXP, SEXP e0_oldSEXP, SEXP gammaSEXP, SEXP KbigSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type e0(e0SEXP);
    Rcpp::traits::input_parameter< double >::type e0_old(e0_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const int >::type Kbig(KbigSEXP);
    rcpp_result_gen = Rcpp::wrap(log_e0_gamma_ratio(e0, e0_old, gamma, Kbig));
    return rcpp_result_gen;
END_RCPP
}
// log_eqd1_ratio
double log_eqd1_ratio(double e0, double e0_old, arma::vec cvec, const int Kbig);
RcppExport SEXP _sparselca_log_eqd1_ratio(SEXP e0SEXP, SEXP e0_oldSEXP, SEXP cvecSEXP, SEXP KbigSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type e0(e0SEXP);
    Rcpp::traits::input_parameter< double >::type e0_old(e0_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cvec(cvecSEXP);
    Rcpp::traits::input_parameter< const int >::type Kbig(KbigSEXP);
    rcpp_result_gen = Rcpp::wrap(log_eqd1_ratio(e0, e0_old, cvec, Kbig));
    return rcpp_result_gen;
END_RCPP
}
// log_prior_ratio
double log_prior_ratio(double e0, double e0_old, double ae, double be);
RcppExport SEXP _sparselca_log_prior_ratio(SEXP e0SEXP, SEXP e0_oldSEXP, SEXP aeSEXP, SEXP beSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type e0(e0SEXP);
    Rcpp::traits::input_parameter< double >::type e0_old(e0_oldSEXP);
    Rcpp::traits::input_parameter< double >::type ae(aeSEXP);
    Rcpp::traits::input_parameter< double >::type be(beSEXP);
    rcpp_result_gen = Rcpp::wrap(log_prior_ratio(e0, e0_old, ae, be));
    return rcpp_result_gen;
END_RCPP
}
// log_conditional_likelihood_mvn_cp
arma::mat log_conditional_likelihood_mvn_cp(const arma::mat Y, const int K, const arma::mat Thetam, const arma::cube Sigmaa);
RcppExport SEXP _sparselca_log_conditional_likelihood_mvn_cp(SEXP YSEXP, SEXP KSEXP, SEXP ThetamSEXP, SEXP SigmaaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Thetam(ThetamSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type Sigmaa(SigmaaSEXP);
    rcpp_result_gen = Rcpp::wrap(log_conditional_likelihood_mvn_cp(Y, K, Thetam, Sigmaa));
    return rcpp_result_gen;
END_RCPP
}
// expected_log_likelihood_mvn_cp
arma::mat expected_log_likelihood_mvn_cp(const arma::mat Y, const arma::mat T, const int K, const arma::mat Theta, const arma::cube Sigma);
RcppExport SEXP _sparselca_expected_log_likelihood_mvn_cp(SEXP YSEXP, SEXP TSEXP, SEXP KSEXP, SEXP ThetaSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type T(TSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_log_likelihood_mvn_cp(Y, T, K, Theta, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// expected_log_p_allocations
arma::vec expected_log_p_allocations(const arma::mat Tmat, const int K, const arma::vec etav);
RcppExport SEXP _sparselca_expected_log_p_allocations(SEXP TmatSEXP, SEXP KSEXP, SEXP etavSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Tmat(TmatSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type etav(etavSEXP);
    rcpp_result_gen = Rcpp::wrap(expected_log_p_allocations(Tmat, K, etav));
    return rcpp_result_gen;
END_RCPP
}
// log_likelihood_mvn_cp
arma::vec log_likelihood_mvn_cp(const arma::mat Y, const int K, const arma::mat Thetam, const arma::cube Sigmaa, const arma::vec eta);
RcppExport SEXP _sparselca_log_likelihood_mvn_cp(SEXP YSEXP, SEXP KSEXP, SEXP ThetamSEXP, SEXP SigmaaSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Thetam(ThetamSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type Sigmaa(SigmaaSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(log_likelihood_mvn_cp(Y, K, Thetam, Sigmaa, eta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sparselca_log_ddirichlet", (DL_FUNC) &_sparselca_log_ddirichlet, 2},
    {"_sparselca_sample_covariance", (DL_FUNC) &_sparselca_sample_covariance, 1},
    {"_sparselca_Mahalanobis", (DL_FUNC) &_sparselca_Mahalanobis, 3},
    {"_sparselca_dmvnrm_arma", (DL_FUNC) &_sparselca_dmvnrm_arma, 5},
    {"_sparselca_dmvnorm_arma_prec", (DL_FUNC) &_sparselca_dmvnorm_arma_prec, 4},
    {"_sparselca_mvrnormArma", (DL_FUNC) &_sparselca_mvrnormArma, 3},
    {"_sparselca_log_dwishart_sample", (DL_FUNC) &_sparselca_log_dwishart_sample, 3},
    {"_sparselca_log_ptheta", (DL_FUNC) &_sparselca_log_ptheta, 4},
    {"_sparselca_log_dirichlet", (DL_FUNC) &_sparselca_log_dirichlet, 2},
    {"_sparselca_wishartArma_int", (DL_FUNC) &_sparselca_wishartArma_int, 2},
    {"_sparselca_iwishartArma_int", (DL_FUNC) &_sparselca_iwishartArma_int, 2},
    {"_sparselca_logsumexp_vec", (DL_FUNC) &_sparselca_logsumexp_vec, 1},
    {"_sparselca_logsumexp", (DL_FUNC) &_sparselca_logsumexp, 2},
    {"_sparselca_normalize_vec", (DL_FUNC) &_sparselca_normalize_vec, 1},
    {"_sparselca_normalize_logvec", (DL_FUNC) &_sparselca_normalize_logvec, 1},
    {"_sparselca_nvec_count", (DL_FUNC) &_sparselca_nvec_count, 2},
    {"_sparselca_count_if", (DL_FUNC) &_sparselca_count_if, 2},
    {"_sparselca_count_classes", (DL_FUNC) &_sparselca_count_classes, 2},
    {"_sparselca_table_by_class", (DL_FUNC) &_sparselca_table_by_class, 3},
    {"_sparselca_update_Sigma_ind_mcmc", (DL_FUNC) &_sparselca_update_Sigma_ind_mcmc, 6},
    {"_sparselca_update_Sigma_mcmc", (DL_FUNC) &_sparselca_update_Sigma_mcmc, 7},
    {"_sparselca_update_eta", (DL_FUNC) &_sparselca_update_eta, 3},
    {"_sparselca_rg_check", (DL_FUNC) &_sparselca_rg_check, 2},
    {"_sparselca_update_class_probs_mvn_cp", (DL_FUNC) &_sparselca_update_class_probs_mvn_cp, 5},
    {"_sparselca_sample_S_mvn_mcmc", (DL_FUNC) &_sparselca_sample_S_mvn_mcmc, 5},
    {"_sparselca_update_theta_ind_mcmc", (DL_FUNC) &_sparselca_update_theta_ind_mcmc, 6},
    {"_sparselca_log_e0_gamma_ratio", (DL_FUNC) &_sparselca_log_e0_gamma_ratio, 4},
    {"_sparselca_log_eqd1_ratio", (DL_FUNC) &_sparselca_log_eqd1_ratio, 4},
    {"_sparselca_log_prior_ratio", (DL_FUNC) &_sparselca_log_prior_ratio, 4},
    {"_sparselca_log_conditional_likelihood_mvn_cp", (DL_FUNC) &_sparselca_log_conditional_likelihood_mvn_cp, 4},
    {"_sparselca_expected_log_likelihood_mvn_cp", (DL_FUNC) &_sparselca_expected_log_likelihood_mvn_cp, 5},
    {"_sparselca_expected_log_p_allocations", (DL_FUNC) &_sparselca_expected_log_p_allocations, 3},
    {"_sparselca_log_likelihood_mvn_cp", (DL_FUNC) &_sparselca_log_likelihood_mvn_cp, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_sparselca(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
