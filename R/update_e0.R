#' @param e0.current real; current value of the e0.
#' @param gamma.current vector; current value of the gamma (mixture weights)
#' @param delta real; standard deviation of normal random walk
#' @param ae.h real; prior hyperparameter with e0~gamma(shape=ae.h, rate=be.h)
#' @param be.h real; prior hyperparameter with e0~gamma(shape=ae.h, rate=be.h)
#' @return updated value of e0
#'
#' @export
update_e0 <- function(e0.current, gamma.current, KK, delta=0.01, ae.h, be.h)
{
  # generate proposal
  e0.proposed <- rnorm(1,e0.current,delta)
  accept <- 0
  if( e0.proposed > 0)
  {
    # note: log prior ratio uses scale parameterization
    logr <- log_e0_gamma_ratio(e0.proposed, e0.current, gamma.current,  KK) + log_prior_ratio(e0.proposed, e0.current,ae.h, 1/be.h)
    accept <- as.numeric(log(runif(1))<min(logr,0))
  }
  return(c( c(e0.current,e0.proposed)[accept+1], accept))

}
