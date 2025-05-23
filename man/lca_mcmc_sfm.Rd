% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gibbs_sfm.R
\name{lca_mcmc_sfm}
\alias{lca_mcmc_sfm}
\title{MCMC for binary sparse LCA}
\usage{
lca_mcmc_sfm(
  YY,
  C,
  nit,
  nthin,
  nburn,
  prior.list,
  tuning.list,
  Rvec = NULL,
  permute = TRUE
)
}
\arguments{
\item{YY}{Numeric matrix with n rows J columns.}

\item{C}{integer; number of classes}

\item{nit}{integer; number of iteractions after burnin, thin}

\item{nthin}{integer; number of samples to thin}

\item{nburn}{integer; number of burn-in samples}

\item{prior.list}{list with all prior hyperparameters}

\item{tuning.list}{list will all tuning parameters}
}
\value{
list with all posterior samples
}
\description{
MCMC for binary sparse LCA
}
