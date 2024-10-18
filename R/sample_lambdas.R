#### Sample from full conditional of shrinkage factors
#### LCA project
#### September 2024
#### Erica M. Porter

# theta.current is a CxJ matrix of current component means
# mu0 is a Jx1 vector of prior means for components
# C is the number of latent classes
# Rvec is a Jx1 vector of ranges for the items (fixed at 1 for now)
# nu1 and nu2 are fixed hyperparameters for the shrinkage prior

#' @importFrom rmutil rginvgauss
#'

update_lambda <- function(theta.current, mu0, C, Rvec, nu1, nu2){
  J <- length(Rvec)
  avec <- rep(2*nu2, J)
  bvec <- apply((theta.current-mu0)^2/Rvec^2, 2, sum)
  Pc <- nu1 - C/2

  m.arg <- sqrt(avec/bvec)
  s.arg <- 1/bvec
  f.arg <- Pc

  lambda.sample <- rmutil::rginvgauss(n=J, m=m.arg, s=s.arg, f=f.arg)

  return(lambda.sample)
}

