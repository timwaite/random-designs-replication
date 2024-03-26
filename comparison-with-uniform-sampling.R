# Comparison with uniform random sampling under specific discrepancies 
#  discrepancies come from Appendix of Nie et al (2017), Can. J. Stat 46, 104-122.

# Example: simple linear model, one factor, X=[-1,1]

source("nie-discrepancies.R")

f <- function(x) c(1,x)
U = diag(c(1,1/3))
tau=1

# initialize discrepancy, use tau2 =1 as can get results for other tau2 by transformation
nie_discrep_linear = Vectorize(InitializeNieDiscrepancy1D(f=f, p=2, U=U, tau2=1, bandwidth=0.01)$discrepancy)

## Verify orthogonality and L2-norm 
#integrate(function(x) nie_discrep_linear(x), lower=-1,upper=1)
#integrate(function(x) nie_discrep_linear(x)*x, lower=-1, upper=1)
#integrate(function(x) nie_discrep_linear(x)**2, lower=-1, upper=1)

M <- function(xi) {
  Fmat <- apply(xi,1,f)
  t(Fmat)%*%Fmat
}

ISE = function(xi, discrep, f, U) {
  psi = discrep(xi)
  Fmat = t(sapply(xi,f))
  Mxi = t(Fmat) %*% Fmat
  Kmat <- 2*Fmat%*%solve(Mxi) %*% U %*% solve(Mxi) %*% t(Fmat)
  ISB_quad_form = c(t(psi) %*% Kmat %*% psi)
  IV = c(2*sum(diag(U %*% solve(Mxi))))
  return(c(IV=IV, ISB_quad_form=ISB_quad_form))
}

#
# sample uniform design realisations 
#
Nmc = 100000
xi_unif <- list()
for (i in 1:Nmc) {
  xi_unif[[i]] <- runif(8,min=-1,max=1)
}

#
# compute the mean IV and mean ISB quadratic form by averaging over design realisations 
#
unif_results <- t(sapply(xi_unif, ISE, nie_discrep_linear, f, U))
(unif_results_summary <- apply(unif_results,2,mean))

#
# compare performance of uniform random design to random translation designs
#
library(tidyverse)

# recall best RT design for different tau2
heuristics_results = read_csv("heuristics-results.csv")
best_Psis <- heuristics_results %>% 
  filter(model=="Linear (1 factor)") %>% 
  group_by(tau2) %>% 
  filter(minPsi==min(minPsi)) %>%
  ungroup()

# compute efficiency bounds for different tau2
best_Psis %>% transmute(tau2=tau2,
                        MIV_RT=MIV, 
                        MIV_unif = unif_results_summary["IV"],
                        MISB_RT=MISB,
                        MISB_unif = tau2 + tau2 * unif_results_summary["ISB_quad_form"],
                        Psi_RT=minPsi,
                      IMSP_unif = unifMIV + unifMISB,
                     Efficiency_unif_UB = minPsi/unifIMSPE) -> best_Psis

View(best_Psis)

xtable( best_Psis %>% slice( seq(from=2,to=90,by=11) ))

#
# table for paper, Section 4.3 
#
library(xtable)
xtable(t(best_Psis %>% slice( seq(from=2,to=90,by=11) ) %>% select(tau2, minPsi, unifImspe, unifEfficiencyUB)))

#
# plot of design realisations
#

delta_split <- heuristics_results %>% filter(model=="Linear (1 factor)", tau2==1, heuristic=="Split") %>% select(optDelta) %>% as.numeric()
(xi_split <-  rep(c(-1,1),each=4) + delta_split/2*c(1,3,5,7,-7,-5,-3,-1) + delta_split/2* runif(1,min=-1,max=1))

plot(xi_unif[[1]], rep(0,8), pch=4, xlim=c(-1,1))
points(xi_split, rep(1,8), pch=2)

