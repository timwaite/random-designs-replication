#
#  Supporting material for "Replication in random translation designs for model-robust prediction"  
# by Waite, T.W. (2024+)

#
# Additional example - 10 factors, 32 runs, main effects model 
#  

source("obj-funs-replication.R")
library(FrF2)
library(lhs)

ffd1 <- attr(FrF2(nruns=16, nfactors=10), "desnum")
ffd2 <- attr(FrF2(nruns=32, nfactors=10), "desnum")

#
# Define heuristic strategies to compare 
#

lin_10factor_ffd1_retain <- list(
  meandesign=function(delta) {
    rbind(ffd1)*(1-delta/2)
  },
  r=rep(2,16),
  model="Linear (10 factor)",
  heuristic="FFD1-Retain",
  deltaBound=1
)

lin_10factor_ffd1_split <- list(
  meandesign=function(delta) {
    rbind(
      ffd1*(1-delta/2), 
      ffd1%*%diag(c(1-3*delta/2, rep(1-delta/2,9)))
    )},
  r=rep(1,32),
  model="Linear (10 factor)",
  heuristic="FFD1-Split",
  deltaBound=0.5
)

lin_10factor_ffd2 <- list(
  meandesign=function(delta) {
    (1-delta/2)*ffd2
  },
  r=rep(1,32),
  model="Linear (10 factor)",
  heuristic="FFD2",
  deltaBound=1
)


#
# set up model 
#

lin_10factor_model = function(x) c(1,x)
A_lin_10factor = diag(c(2**10, rep(2**10/3,10)))

heuristics_list = list(lin_10factor_ffd1_retain, 
                       lin_10factor_ffd1_split,
                       lin_10factor_ffd2)

fullfactorial <- as.matrix( 
  expand.grid(
    rep( list(c(-1,1)), 10 ) 
    ) 
  )


Tmax_10factor <- rbind( 2*randomLHS(1000, 10) -1 , fullfactorial)
Tmc_10factor <- 2*randomLHS(1000, 10)-1

lin_10factor_stats <- NULL

tau2s = c(0.001,0.01,0.1)

for (tau2 in tau2s ) {
  for (i in 1:3 ) {
    heuristic= heuristics_list[[i]]
    cat(heuristic$heuristic, "\t tau2 = ", tau2, "\n")
    lin_10factor_stats <- rbind(lin_10factor_stats,
                                optimal_heuristic_stats(heuristic, 
                                                        tau2=tau2, 
                                                        lin_10factor_model, 
                                                        A_lin_10factor, 
                                                        Tmc_10factor, 
                                                        Tmax_10factor))
  }
}

lin_10factor_stats
# write.csv(lin_10factor_stats, file="10factor-results.csv")

library(tidyverse)
lin_10factor_stats %>% group_by(tau2) %>% mutate(Efficiency = min(minPsi)/minPsi) %>% select(tau2,heuristic,optDelta,minPsi,MIV,MISB,Efficiency)-> results10factor2

library(xtable)
xtable(results10factor2)
