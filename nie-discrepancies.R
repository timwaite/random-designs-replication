
RestrictedGaussianBump1D <- function(x, x_star=x_star, bandwidth=bandwidth) {
  num <-   dnorm(x, x_star, sd= bandwidth)
  denom <-  pnorm((1-x_star)/bandwidth) - pnorm((-1-x_star)/bandwidth)
  num/denom *(x<=1) * (x>=-1)
}

InitializeNieDiscrepancy1D <- function(f,p,U,tau2, lower=-1,upper=1, bandwidth=1e-2) {
  # NB training and test density both U[-1,1]
  #
  # find x_star
  # - NB using optim as optimize seems to find the minimum not maximum 
  afun <- function(x) { t(f(x)) %*% solve(U) %*% f(x) * (x>=lower) * (x<=upper) }
  x_star <- optim(runif(1), Vectorize(afun), lower=lower, upper=upper, method="L-BFGS-B", control=list(fnscale=-1))$par
  
  #
  # compute L2 projection coefficients
  #
  alpha <- vector(length=p)
  for (j in 1:p) {
    tmpfun = function(x) f(x)[j] * sqrt(2*RestrictedGaussianBump1D(x, x_star, bandwidth))*1/2
    alpha[j] <- integrate( Vectorize(tmpfun), lower=lower, upper=upper )$value
  }
  coefs <- solve(U)%*% alpha
  
  l2_proj <- function(x) c(t(f(x))%*%coefs)
  
  discrep <- function(x) tau*(sqrt(2*RestrictedGaussianBump1D(x, x_star, bandwidth)) - l2_proj(x))/(c(sqrt(1-t(alpha)%*%solve(U)%*%alpha)))/sqrt(2)
  # NB sqrt(2) factor relates to difference between Nie's notation and ours - their L2 norm is wrt a prob density (in this case uniform on [-1,1]), ours is wrt Lebesgue measure 
  return(list(discrepancy=discrep, alpha=alpha, coefs=coefs))
}
 

#nie_discrep_linear = Vectorize(InitializeNieDiscrepancy1D(f=function(x) c(1,x), p=2, U=diag(c(1,1/3)))$discrepancy)
#Uquad = diag(c(1,1/3,1/5))
#Uquad[1,3] <- 1/3 -> Uquad[3,1]
#nie_discrep_quad = Vectorize(InitializeNieDiscrepancy1D(f=function(x) c(1,x,x^2), p=3, U=Uquad)$discrepancy)

