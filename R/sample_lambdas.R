#### Sample from full conditional of shrinkage factors
#### LCA project
#### December 2024
#### Erica M. Porter

# theta.current is a CxJ matrix of current component means
# mu0 is a Jx1 vector of prior means for components
# C is the number of latent classes
# Rvec is a Jx1 vector of ranges for the items (fixed at 1 for now)
# nu1 and nu2 are fixed hyperparameters for the shrinkage prior

#' @importFrom GeneralizedHyperbolic rgig
#'

update_lambda <- function(theta.current, mu0, C, Rvec, nu1, nu2){
  J <- length(Rvec)

lambda.sample <- rep(NA,J)

  avec <- rep(2*nu2, J)
  Pc <- nu1 - C/2

  bvec <- rep(NA, J)
  for(j in 1:J){
    bvec[j] <- sum(((theta.current[,j] - mu0[j])^2)/(Rvec[j]^2))
  }

  for(j in 1:J){
    param <- c(avec[j], bvec[j], Pc)
    lambda.sample[j] <- GeneralizedHyperbolic::rgig(1, param=param)
  }

  return(lambda.sample)
}

