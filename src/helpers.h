#ifndef helpers
#define helpers

# include <RcppArmadillo.h>
# include <math.h>
# include <Rmath.h>
# include <R.h>
# include <Rcpp.h>
# include <iostream>

// FUNCTION DECLARATIONS
// ---------------------
double logsumexp_vec(arma::vec avec);
double logsumexp(const double a, const double b);

arma::mat normalize_A(arma::mat A);
arma::mat normalize_A_reg(arma::mat A);
arma::mat normalize_vec(arma::vec v);
arma::mat normalize_logvec(arma::vec lv);


arma::vec nvec_count(arma::vec svec, int k);
arma::vec ordered_remove(arma::vec list, int index);

arma::vec ordered_insert_next(arma::vec list);


int ordered_next(arma::vec list);


bool in_vec_int(int n, arma::vec v);


int which_vec_int(int n, arma::vec v);



arma::vec int_seq(const int k);



int count_if(arma::vec c_vec, const int cval);


arma::vec count_classes(arma::vec c_vec, const int Kbig);


arma::cube table_by_class(arma::mat YY, arma::vec cvec,
                          const int Kbig);



#endif
