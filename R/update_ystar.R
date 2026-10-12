update_ystar_old <- function(theta.current, Sigma.current,A,B,c.current,
                         ystar.previous)
{
  N <- nrow(A)
  ystar.current <- matrix(NA, nrow=N,ncol=ncol(A))


  #ystar.current <- update_Ystar_mcmc(c.current, theta.current,Sigma.current,  A ,B)
   for(i in 1:N)
   {
     ci <- c.current[i]
  #   #ystar <- tmvtnorm::rtmvnorm(1, mean=theta.current[ci,],sigma=Sigma.current[,,ci],
  #   #                            lower=A[i,],upper=B[i,],algorithm = "gibbs")
     ystar <- update_ystar_mcmc(theta.current[ci,],Sigma.current[,,ci],
                                A[i,] ,B[i,])
  #   #ystar <- tmvmixnorm::rtmvn(1, Mean=theta.current[ci,],Sigma=Sigma.current[,,ci],
  #   #                            lower=A[i,],upper=B[i,],int=ystar.previous[i,])
  #
   ystar.current[i,] <- ystar
   }
  ystar.current
}
