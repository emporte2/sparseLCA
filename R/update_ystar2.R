update_ystar <- function(theta.current, Sigma.current,A,B,c.current,
                         ystar.previous)
{
  N <- nrow(A)
  ystar.current <- matrix(NA, nrow=N,ncol=ncol(A))
  C <- nrow(theta.current)
  J <- ncol(theta.current)
  mult.mean <- array(NA, dim=c(C,J,J-1))
  conditional.sig <- matrix(NA, nrow=C, ncol=J)
  for(c in 1:C)
  {
    for(j in 1:J)
    {
    mult.mean[c,j,] <- Sigma.current[j,-j,c]%*%solve(Sigma.current[-j,-j,c])
    conditional.sig[c,j] <- Sigma.current[j,j,c] - crossprod(mult.mean[c,j,], Sigma.current[-j,j,c])
    }
  }


   for(i in 1:N)
   {
     ci <- c.current[i]
     for(j in 1:J)
     {
       mean1 <- theta.current[ci,j] + mult.mean[c,j,]%*%(ystar.previous[i,-j] - theta.current[ci,-j])
       ystar.current[i,j] <- rtruncnorm(1, a=A[i,j], b=B[i,j], mean = mean1, sd = sqrt(conditional.sig[c,j]))
     }

   }
  ystar.current
}
