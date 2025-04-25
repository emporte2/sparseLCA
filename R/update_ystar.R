update_ystar <- function(theta.current, Sigma.current,A,B,c.current)
{
  N <- nrow(A)
  ystar.current <- matrix(NA, nrow=N,ncol=ncol(A))

  for(i in 1:N)
  {
    ci <- c.current[i]
    ystar <- tmvtnorm::rtmvnorm(1, mean=theta.current[ci,],sigma=Sigma.current[,,ci],
                                lower=A[i,],upper=B[i,],algorithm = "gibbs")
    ystar.current[i,] <- ystar
  }
  ystar.current
}
