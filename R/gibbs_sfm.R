#' @title MCMC for binary sparse LCA
#'
#'
#' @param YY Numeric matrix with n rows J columns.
#' @param C integer; number of classes
#' @param nit integer; number of iteractions after burnin, thin
#' @param nthin integer; number of samples to thin
#' @param nburn integer; number of burn-in samples
#' @param prior.list list with all prior hyperparameters
#' @param tuning.list list will all tuning parameters
#' @return list with all posterior samples
#'
#'
#' @importFrom tmvtnorm rtmvnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib sparselca
#'
#' @export
lca_mcmc_sfm <- function(YY,  C, nit, nthin, nburn,
                         prior.list, tuning.list, Rvec=NULL, permute=TRUE, order.gamma=FALSE,
                         update.theta=TRUE,update.Sigma=TRUE,update.lambda=TRUE,
                         update.gamma=TRUE, update.e0=TRUE)
{
  N <- nrow(YY)
  J <- ncol(YY)

  if(permute & order.gamma)
  {
    stop("Only one of permute and order.gamma can be TRUE.")
  }

  # create truncation set for data YY
  A <- matrix(-Inf,nrow=N,ncol=J)
  B <- rep(Inf,nrow=N,ncol=J)
  A <- t(apply(YY,1,function(x){out <- rep(-Inf,J); out[which(x==1)] <- 0; out}))
  B <- t(apply(YY,1,function(x){out <- rep(Inf,J); out[which(x==0)] <- 0; out}))

  # objects for prior hyperparameters
  mu.current <- prior.list$mu0 # will eventually be updated
  #B0 <- prior.list$B0
  #e0 <- prior.list$e
  ae.0 <- prior.list$ae
  be.0 <- ae.0*C^2 # originally ae.0*C: testing a higher b
  nu1 <- prior.list$nu1
  nu2 <- prior.list$nu2
  delta.e0 <- tuning.list$delta
  d0 <- 2.5 + (J-1)/2
  D0 <- prior.list$D0

  # create objects to store results
  ystar.samples <- array(NA,dim=c(nit,N,J))
  c.samples <- matrix(NA, nrow=nit,ncol=N) # latent allocations
  theta.samples <- array(NA, dim=c(nit,C,J))
  Sigma.inv.samples <- array(NA, dim=c(nit,J,J,C))
  Sigma.samples <- array(NA, dim=c(nit,J,J,C))
  cplus.samples <- matrix(NA, nrow=nit) # latent allocations
  e0.samples <- matrix(NA,nrow=nit,ncol=2)
  lambda.samples <- matrix(NA, nrow=nit, ncol=J)

  gamma.samples <- matrix(NA, nrow=nit,ncol=C) # mixture weights
  pi.samples <- array(NA, dim=c(nit, C,J))
  kplus.samples <- rep(NA,nit) # number of non empty components

  # initialize parameters
  cv.current <- sample(1:C,size=N,replace=TRUE)
  nv.current <- count_classes(cv.current,C)
  Sigma.inv.current <- array(NA, dim=c(J,J,C))
  for(c in 1:C)
  {
    Sigma.inv.current[,,c] <- diag(rep(1,J)) # identity is default value
  }
  Sigma.current <- array(NA, dim=c(J,J,C))
  for(c in 1:C)
  {
    Sigma.current[,,c] <- Sigma.inv.current[,,c]
  }

  if(is.null(prior.list$e0))
  {
  e0.current.s <- c(rgamma(1,ae.0, rate=be.0),0)
  e0.current <- e0.current.s[1]
  } else{
    e0.current.s <- c(1/C,0)
    e0.current <- e0.current.s[1]
  }

  if(!update.gamma)
  {
    gamma.current <- rep(1/C,C)
  }


  theta.current <- matrix(rnorm(J*C),nrow=C,ncol=J)
  if(!update.theta)
  {
    theta.current <- matrix(rep(qnorm(apply(YY,2,mean)),C),nrow=C,ncol=J,byrow=TRUE)
  }
  ystar.current <- matrix(rnorm(N*J),nrow=N,ncol=J)

  if(!is.null(prior.list$lambda))
  {
    lambda.current <- prior.list$lambda
    Lambda.temp <- diag(sqrt(lambda.current))
    B0.current <- Lambda.temp%*%diag(Rvec)%*%Lambda.temp
    B0.current <- diag(Rvec)
  } else{
  B0.current <- diag(Rvec)
  lambda.current <-rep(1,J)
  }

  ## add everything else
  it <- 1
  it2 <- 1

  while( it/nthin <= nit)
  {
    # update y star
    for(i in 1:N)
    {
      ci <- cv.current[i]
      ystar <- tmvtnorm::rtmvnorm(1, mean=theta.current[ci,],sigma=Sigma.current[,,ci],
                        lower=A[i,],upper=B[i,],algorithm = "gibbs",burn.in=10,thinning=5)

      ystar.current[i,] <- ystar
    }

    if(update.theta)
    {
    # update theta
    theta.current <- update_theta_ind_mcmc(ystar.current,cv.current, C, Sigma.inv.current,
                                           mu.current, B0.current) # takes precision matrix

    if(permute)
    {
      theta.current <- theta.current[sample(1:C,C),]
    }
    }
    # could write a function here to transform thetas into pi's for FYI

    if(update.gamma)
    {
    # update gamma
    gamma.current <- update_eta(C, cv.current, e0.current[1])
    if(order.gamma)
    {
      gamma.current <- sort(gamma.current)
    }
    }


    if(sum(is.na(gamma.current))>0){stop('na in gamma.current')}
    if(sum(is.na(theta.current))>0){stop('na in theta.current')}
    if(sum(is.na(ystar.current))>0){stop('na in ystar.current')}

    # update C - latent allocations
    cv.current <- sample_S_mvn_mcmc(ystar.current,C, theta.current,Sigma.current, gamma.current)
    nv.current <- count_classes(cv.current,C)
    cplus.current <- sum(nv.current>0)
    if(sum(is.na(cv.current))>0){stop('na in cv.current')}

    # update B matrix - coming soon
    if(update.lambda)
    {
    lambda.current <- update_lambda(theta.current, mu.current, C, Rvec, nu1, nu2)
    Lambda <- diag(sqrt(lambda.current))
    B0.current <- Lambda%*%diag(Rvec)%*%Lambda
    }

    if(update.Sigma)
    {
    # update Sigma
    Sigma.inv.current <- update_Sigma_ind_mcmc(ystar.current, cv.current,C, theta.current, d0, D0)
    for(c in 1:C)
    {
      Sigma.current[,,c] <- solve(Sigma.inv.current[,,c])
    }
    }

    if(update.e0)
    {
    # update e0
    e0.current.s <- update_e0(e0.current, gamma.current, C, delta=delta.e0, ae.0, be.0)
    e0.current <- e0.current.s[1]
    }

    if(it2 > nburn)
    {
      if(it%%nthin==0)
      {

        cat("iteration:",it2,"\n")

        c.samples[it/nthin,] <- cv.current
        gamma.samples[it/nthin,] <- gamma.current
        pi.samples[it/nthin,,] <- pnorm(theta.current)
        cplus.samples[it/nthin] <- cplus.current
        ystar.samples[it/nthin,,] <- ystar.current
        theta.samples[it/nthin,,] <- theta.current
        Sigma.inv.samples[it/nthin,,,] <- Sigma.inv.current
        Sigma.samples[it/nthin,,,] <- Sigma.current
        e0.samples[it/nthin,] <- e0.current.s
        lambda.samples[it/nthin,] <- lambda.current


      }
      it <- it+1
    }
    it2 <- it2+1
  }
  return(list(c=c.samples, gamma=gamma.samples, pi=pi.samples,
              cplus=cplus.samples, ystar=ystar.samples, theta=theta.samples,
              Sigma=Sigma.samples, e0=e0.samples))
}
