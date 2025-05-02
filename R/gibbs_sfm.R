#' @title MCMC for binary sparse LCA
#'
#'
#' @param YY Numeric matrix with n rows J columns.
#' @param C integer; number of classes
#' @param nit integer; number of iterations across chains after burnin, thin
#' @param nthin integer; number of samples to thin
#' @param nburn integer; number of burn-in samples
#' @param nchain integer; number of chains
#' @param prior.list list with all prior hyperparameters
#' @param tuning.list list will all tuning parameters
#' @return list with all posterior samples
#'
#'
#' @importFrom tmvtnorm rtmvnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom Rcpp evalCpp
#' @importFrom GeneralizedHyperbolic rgig
#' @useDynLib sparselca
#'
#' @export
lca_mcmc_sfm <- function(YY,  C, nit, nthin=1, nburn, nchain=1,
                         prior.list, tuning.list, permute=TRUE, order.gamma=FALSE,
                         update.theta=TRUE,update.Sigma=FALSE,update.lambda=TRUE, update.mu=TRUE,
                         update.gamma=TRUE, update.e0=TRUE, update.D0=FALSE)
{
  N <- nrow(YY)
  J <- ncol(YY)

  nit.per.chain <- floor(nit/nchain)
  nit <- nit.per.chain*nchain
  cat("Total iterations split among ",nchain, " chains:", nit,
      "\n iterations per chain: ",nit.per.chain, "\n")

  if(permute & order.gamma)
  {
    stop("Only one of permute and order.gamma can be TRUE.")
  }

  # create truncation set for data YY
  A <- matrix(-Inf,nrow=N,ncol=J)
  B <- rep(Inf,nrow=N,ncol=J)
  A <- t(apply(YY,1,function(x){out <- rep(-Inf,J); out[which(x==1)] <- 0; out}))
  B <- t(apply(YY,1,function(x){out <- rep(Inf,J); out[which(x==0)] <- 0; out}))

  delta.e0 <- tuning.list$delta


  # create objects to store results
  ystar.samples <- array(NA,dim=c(nit,N,J))
  c.samples <- matrix(NA, nrow=nit,ncol=N) # latent allocations
  theta.samples <- array(NA, dim=c(nit,C,J))
  mu.samples <- matrix(NA, nrow=nit, ncol=J)
  Sigma.inv.samples <- array(NA, dim=c(nit,J,J,C))
  Sigma.samples <- array(NA, dim=c(nit,J,J,C))
  cplus.samples <- matrix(NA, nrow=nit) # latent allocations
  e0.samples <- matrix(NA,nrow=nit,ncol=2)
  lambda.samples <- matrix(NA, nrow=nit, ncol=J)

  gamma.samples <- matrix(NA, nrow=nit,ncol=C) # mixture weights
  pi.samples <- array(NA, dim=c(nit, C,J))
  kplus.samples <- rep(NA,nit) # number of non empty components

  for( ch in 1:nchain)
  {
  # initialize parameters
  cv.current <- sample(1:C,size=N,replace=TRUE)
  nv.current <- count_classes(cv.current,C)
  ystar.current <- matrix(rnorm(n*J), nrow=n,ncol=J)


  hyp <- default.hyp(prior.list, YY,C)
  be.0 <- hyp$be0
  ae.0 <- hyp$ae0
  e0.current <- hyp$e0
  d0 <- hyp$d0
  D0.current <- hyp$D0
  g0 <- hyp$g0
  G0 <- hyp$G0
  lambda.current <- hyp$lambda
  r0.vec <- hyp$r0.vec
  nu1 <- hyp$nu1
  nu2 <- hyp$nu2
  Sigma.current <- hyp$Sigma
  Sigma.inv.current <- hyp$Sigma.inv
  theta.current <- hyp$theta
  gamma.current <- hyp$gamma
  mu.current <- hyp$mu
  B0.current <- hyp$B
  B0.inv.current <- solve(hyp$B)



  ## add everything else
  it <- 1
  it2 <- 1

  while( it/nthin <= nit.per.chain)
  {
    # update y star
    ystar.current <- update_ystar(theta.current, Sigma.current, A,B, cv.current)

    if(update.theta)
    {
    # update theta
    theta.current <- update_theta_ind_mcmc(ystar.current,cv.current, C, Sigma.inv.current,
                                           mu.current, B0.inv.current) # takes precision matrix

    if(permute)
    {
      theta.current <- theta.current[sample(1:C,C),]
    }
    }

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
    lambda.current <- update_lambda(theta.current, mu.current, C, r0.vec, nu1, nu2)
    Lambda <- diag(sqrt(lambda.current))
    B0.current <- Lambda%*%diag(r0.vec)%*%Lambda
    Lambda.i <- diag(sqrt(1/lambda.current))
    B0.inv.current <- Lambda.i%*%diag(1/r0.vec)%*%Lambda.i
    }

    if(update.mu)
    {
    mu.current <- mvtnorm::rmvnorm(1,apply(theta.current,2,mean), (1/C)*B0.current)
    }

    if(update.Sigma)
    {
    # update Sigma
    Sigma.inv.current <- update_Sigma_ind_mcmc(ystar.current, cv.current,C, theta.current, d0, D0)
    for(c in 1:C)
    {
      Sigma.current[,,c] <- solve(Sigma.inv.current[,,c])
    }
    if(update.D0)
    {
      D0 <- rWishart(1, df=g0+C*d0, Sigma=G0+apply(Sigma.inv.current,2:3,sum))[,,1]
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



        c.samples[nit.per.chain*(ch-1) + it/nthin,] <- cv.current
        gamma.samples[nit.per.chain*(ch-1) + it/nthin,] <- gamma.current
        pi.samples[nit.per.chain*(ch-1) + it/nthin,,] <- pnorm(theta.current)
        cplus.samples[nit.per.chain*(ch-1) + it/nthin] <- cplus.current
        ystar.samples[nit.per.chain*(ch-1) + it/nthin,,] <- ystar.current
        theta.samples[nit.per.chain*(ch-1) + it/nthin,,] <- theta.current
        mu.samples[nit.per.chain*(ch-1) + it/nthin,] <- mu.current
        Sigma.inv.samples[nit.per.chain*(ch-1) + it/nthin,,,] <- Sigma.inv.current
        Sigma.samples[nit.per.chain*(ch-1) + it/nthin,,,] <- Sigma.current
        e0.samples[nit.per.chain*(ch-1) + it/nthin,] <- e0.current.s
        lambda.samples[nit.per.chain*(ch-1) + it/nthin,] <- lambda.current


      }
      it <- it+1
    }
    it2 <- it2+1
  }
  cat("chain:",ch," of ", nchain," complete. \n")

  }
  return(list(c=c.samples, gamma=gamma.samples, pi=pi.samples,
              cplus=cplus.samples, ystar=ystar.samples, theta=theta.samples,
              mu = mu.samples,
              Sigma=Sigma.samples, e0=e0.samples))
}
