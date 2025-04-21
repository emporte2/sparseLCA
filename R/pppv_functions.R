library(mvtnorm)

gen.ppy <- function(output, C, J,n, nit)
{
  ysim <- array(NA, dim=c(nit, n, J))
  for( i in 1:n)
  {
    for(it in 1:nit)
    {
      c.samp <- sample(1:C,1,prob=output$gamma[it,])
      ysim[it,i,] <- rmvnorm(1, output$theta[it,c.samp,], output$Sigma[it,,,c.samp])
    }
  }

  ysim2 <- apply(ysim, 1:3,function(x){as.numeric(x>0)})
  return(ysim2)
}

enumerate.binary.response.patterns <- function(r,J)
{
  outcomes <- list()
  for( j in 1:J)
  {
    outcomes[[j]] <- seq(0,r-1)
  }
  expand.grid(outcomes)
}

log.p.resp.c <- function(resp, item.probs, gamma)
{
  sum(resp*log(item.probs) + (1-resp)*log(1-item.probs)) + log(gamma)

}

count.resp.patterns <- function(YY, resp.patterns)
{
  J <- ncol(YY)
  observed.patterns <- sapply(1:nrow(YY),function(i){ which( apply(rps,1,function(x){(sum(x==YY[i,]))})==J) } )
  observed.counts <- sapply(1:nrow(resp.patterns),function(i){sum(observed.patterns==i)})
  observed.counts
}

ppy.resp.patterns <- function(Ytilde.array, resp.patterns)
{
  out <-t(sapply(1:dim(Ytilde.array)[1],function(it) { count.resp.patterns(Ytilde.array[it,,], resp.patterns) } ))
  out
}


log.cell.probs <- function(param, C.fit, response.patterns)
{
  nit <- nrow(param$gamma)
  log.cell.probM <- matrix(NA, nrow=nit, ncol=nrow(response.patterns))

  for( it in 1:nit)
  {
    log.cell.prob <- rep(NA, nrow(rps))

    for(resp in 1:nrow(rps))
    {
      logp.resp.vec <- rep(0,C.fit)
      for( c.temp in 1:C.fit)
      {
        logp.resp.vec[c.temp] <- log.p.resp.c(response.patterns[resp,], param$pi[it,c.temp,], param$gamma[it,c.temp])

      }
      log.cell.prob[resp] <- logsumexp_vec(logp.resp.vec)
    }
    log.cell.probM[it, ] <- log.cell.prob
  }
  apply(log.cell.probM,2,function(x){mean(exp(x))})
}

g.squared <- function(cell.counts,cell.probs,  n=NULL)
{
  if(is.null(n)) { n <- sum(cell.counts)}
  G2 <- 2*sum( cell.counts*(log(cell.counts/(cell.probs*n )+ 1e-12)) )
  G2
}

g.squared.distribution <- function(Ytilde.rps, cell.probs, n=NULL)
{
  g2 <- rep(NA, nrow(Ytilde.rps))
  for(i in 1:nrow(Ytilde.rps))
  {
  g2[i] <- g.squared(Ytilde.rps[i,],cell.probs,  n=NULL)
  }
  g2
}

posterior.predictive.p <- function(YY, output, C, J,n, nit)
{
  responses <- enumerate.binary.response.patterns(2,J)

  # count frequency of patterns in observed data
  obs.response.patterns <- count.resp.patterns(YY, responses)

  # generate posterior predicted samples
  Y.tilde <- gen.ppy(output, C, J,n, nit)

  # frequency table of response patterns in Y.tilde
  pp.response.patterns <-ppy.resp.patterns(Y.tilde, responses)

  # calculated expected cell counts
  cell.probs <- log.cell.probs(output, C, responses)

  # g2 for y-tilde
  gsquared.tilde <- g.squared.distribution(pp.response.patterns,cell.probs)
  gsquared.obs <- g.squared(obs.response.patterns, cell.probs, n=n)


  # pppv
  mean(gsquared.obs < gsquared.tilde)


}

