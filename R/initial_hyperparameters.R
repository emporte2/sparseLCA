default.hyp <- function(prior.list, Y, C=NULL,
                        update.theta=TRUE, update.lambda=TRUE,
                        update.mu=TRUE, update.Sigma=TRUE,
                        update.D0=TRUE, update.gamma=TRUE,
                        update.e0=TRUE)
{
 J <- ncol(Y)

# theta updates hyperparameters
if(update.theta)
{
  if(is.null(prior.list$r0.vec))
  {
    r0.vec <- rep(0.5, J)
    cat("r0.vec initialized by default. \n")
  } else
  {
    r0.vec <- prior.list$r0.vec
  }

  if(update.lambda)
  {
    if(is.null(prior.list$nu1))
    {
      cat("nu1 initialized by default. \n")
      nu1 <- 0.5
    } else {
      nu1 <- prior.list$nu1

    }
    if( is.null(prior.list$nu2))
    {
      cat("nu2 initialized by default. \n")
      nu2 <- 0.5
    } else {
      nu2 <- prior.list$nu2
    }

  }
  # initialize or set default lambda
  lambda <- rep(1, J)

  if(update.mu)
  {
  if( is.null(prior.list$m0))
  {
    m0 <-  qnorm(apply(Y,2,mean)+1e-12) # data dependent mu value if unspecified
    cat("m0 initialized by default. \n")

  }else {
    m0 <- prior.list$m0
  }
  }
  # initialize mu
  if( is.null(prior.list$mu))
  {
    mu <-  qnorm(apply(Y,2,mean)+1e-12) # data dependent mu value if unspecified
    cat("mu initialized by default. \n")

  } else
  {
    mu <- prior.list$mu
  }

  # initialize B
  B <- diag(sqrt(lambda))%*%diag(r0.vec)%*%diag(sqrt(lambda))


}
  # initialize theta from data
  theta <- matrix(rep(qnorm(apply(Y,2,mean)),C),nrow=C,ncol=J,byrow=TRUE)
# set hyperparameters for Sigma updates
if(update.Sigma)
{
  if(is.null(prior.list$d0))
  {
  d0 <- 2.5 + (J-1)/2
  cat("d0 initialized by default. \n")

  } else
  {
  d0 <- prior.list$d0
  }
  if(is.null(prior.list$D0))
  {
    D0 <- diag(1/(r0.vec^2), J)*100
    cat("D0 initialized by default. \n")

  } else {
    D0 <- prior.list$D0
  }

  if(update.D0)
  {
    if(is.null(prior.list$g0))
    {
      g0 <- 0.5 + (J-1)/2
      cat("g0 initialized by default. \n")

    } else
    {
      g0 <- prior.list$g0
    }
    if(is.null(prior.list$G0))
    {
      G0 <- diag(1/(r0.vec)^2, J)*100*g0/d0
      cat("G0 initialized by default. \n")

    } else
    {
      G0 <- prior.list$G0
    }
  } else
  {
    g0=NULL
    G0=NULL
  }
} else
{
  d0=NULL
  D0=NULL
  g0=NULL
  G0=NULL
}

# initialize Sigma regardless of its update
Sigma.inv.current <- array(NA, dim=c(J,J,C))

for(c in 1:C)
  {
    Sigma.inv.current[,,c] <- diag(rep(1,J)) # identity is default value
  }
Sigma.current <- Sigma.inv.current


# set hyperparameters for gamma updates
if(update.gamma)
{
if(update.e0)
{
  if( is.null(prior.list$ae.0))
  {
    ae.0 <- 10
    cat("ae.0 initialized by default. \n")

  } else{
    ae.0 <- prior.list$ae.0

  }
  if(is.null(prior.list$be.0))
  {
    if(!is.null(C))
    {
      be.0 <- ae.0*C^2
      cat("be.0 initialized by default. \n")

    }else {stop("C must be provided if be0 is null.")}
  } else
  {
    be.0 <- prior.list$be.0
  }
  e0 <- (rgamma(1,ae.0, rate=be.0))
} else
{
  e0 <- rep(1/C, C)
}


} else {
  e0=NULL
  ae.0=NULL
  be.0=NULL
}

# initialize gamma regardless of update
gamma <- rep(1/C,C)



return(list(theta=theta, Sigma=Sigma.current, Sigma.inv = Sigma.inv.current,gamma=gamma,
            lambda=lambda, r0.vec=r0.vec, B=B, nu1=nu1, nu2=nu2, mu=mu,
            e0=e0, ae0=ae.0, be0=be.0, d0=d0, D0=D0, g0=g0, G0=G0))
}
