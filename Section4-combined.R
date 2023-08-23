#
#  Supporting material for "Replication in random translation designs for model-robust prediction"  
# by Waite, T.W. (2023+)


# Section 4 - numerical examples 
#   Driver script to compute results for Test Problems 1 to 4
#    results are saved in heuristics-results.csv 
#   to reproduce plots from results file see Section4-combined-plots.R

source("obj-funs-replication.R")

### Linear 1 factor example set up 

# set up model
lin_1factor_model =function(x) { c(1,x) }
A_lin_1factor = diag(c(2,2/3))

# define heuristics to compare, 8 runs each 

fullsplit_lin_1factor <- list( 
  meandesign=function(delta) matrix( rep(c(-1,1),each=4) + delta/2*c(1,3,5,7,-7,-5,-3,-1), ncol=1), 
  r= rep(1,8),
  model="Linear (1 factor)",
  heuristic="Split",
  deltaBound = 1/4
)

halfsplit_lin_1factor <- list(
  meandesign=function(delta) matrix( rep(c(-1,1), each=2) +  delta/2*c(1,3,-3,-1), ncol=1),
  r=rep(2,4),
  model="Linear (1 factor)",
  heuristic="Half split",
  deltaBound = 0.5
)

retain_lin_1factor <- list(
  meandesign=function(delta) matrix( c(-1,1) + delta/2*c(1,-1), ncol=1),
  r=rep(4,2),
  model="Linear (1 factor)",
  heuristic="Retain",
  deltaBound = 1
)


### Quadratic 1 factor example set up

# set up model
quad_1factor_model = function(x) { c(1,x,x^2) }
A_quad_1factor = diag(c(2,2/3,2/5)); A_quad_1factor[1,3] <- 2/3 -> A_quad_1factor[3,1]


# define heuristics to compare, 4 runs each 
retain_quad_1factor <- list(
  meandesign=function(delta) matrix( c(-1+delta/2,0,1-delta/2) , ncol=1),
  r=c(1,2,1),
  model="Quadratic (1 factor)",
  heuristic="Retain",
  deltaBound = 2/3
)

split_quad_1factor <- list(
  meandesign = function(delta) matrix( 
    c(-1,0,0,1) + delta/2*c(1,-1,1,-1), 
    ncol=1),
  r=c(1,1,1,1),
  model="Quadratic (1 factor)",
  heuristic="Split",
  deltaBound = 0.5
)

### Linear model with 2 factors

# set up model
lin_2factor_model =function(x) { c(1,x) }
A_lin_2factor = diag(c(4,4/3,4/3))

# define heuristics to compare 

fullfactorial <- as.matrix(expand.grid(c(-1,1), c(-1,1)))
doublefactorial <- fullfactorial[rep(1:4,each=2),]

split_lin_2factor <- list( 
  meandesign=function(delta) doublefactorial + 
    delta/2*matrix(
      c(1,3,-1,-3,1,3,-1,-3,
        rep(c(1,-1),each=4)
      ),
      ncol=2), 
  r= rep(1,8),
  model="Linear (2 factor)",
  heuristic="Split",
  deltaBound = 0.5
)

retain_lin_2factor <- list( 
  meandesign=function(delta) fullfactorial + 
    delta/2*matrix(
      c(1,-1,1,-1, rep(c(1,-1),each=2)),
      ncol=2), 
  r= rep(2,4),
  model="Linear (2 factor)",
  heuristic="Retain",
  deltaBound = 1
)

### Quadratic model with 3 factors 

# set up model

quad_3factor_model = function(x) { c(1,x[1],x[2],x[3], x[2]*x[3], x[1]*x[3], x[1]*x[2], x[1]^2, x[2]^2, x[3]^2)}
A_quad_3factor = matrix(0, ncol=10, nrow=10)

indxs <- matrix( c(0,0,0,
                   1,0,0,
                   0,1,0,
                   0,0,1,
                   0,1,1,
                   1,0,1,
                   1,1,0,
                   2,0,0,
                   0,2,0,
                   0,0,2), ncol=3, nrow=10, byrow=T)

INT <- function(idx) { 2/(idx+1) * ( (idx+1) %% 2 == 1) }

for (i in 1:10) {
  for (j in 1:10) {
    A_quad_3factor[i,j] <- prod( sapply(indxs[i,] + indxs[j,], INT) )
  }
}


# set up heuristics to compare

Vopt_des = read.csv("design_3factor.csv")[,2:4]
Vopt_des_unique <- unique(Vopt_des)

clip <- function(x, delta=0.1) {
  x <- as.matrix(x)
  x<-  pmin(  x, 1-delta/2 )
  x <- pmax( x, -1+delta/2)
  return(x)
}

retain_quad_3factor <- list(
  meandesign=function(delta) {
    xi <- as.matrix(Vopt_des_unique)
    xi <- clip(xi, delta)  
    return(xi)},
  r = c(1,1,1,1,1,2,1,1,1,1,1),
  model = "Quadratic (3 factor)",
  heuristic = "Retain",
  deltaBound = 0.585 # for this delta, the support sets for points 7 and 11 (labelled rows 8 and 12 in des_unique - row 7 is missing) begin to overlap
)

split_quad_3factor <- list(
  meandesign=function(delta){
    xi <- as.matrix(Vopt_des)
    xi <- clip(xi,delta)
    xi[6,] = xi[6,] + delta/2*(c(-1,0,0))
    xi[7,] = xi[7,] + delta/2*(c(1,0,0))
    return(xi)
  },
  r=rep(1,12),
  model = "Quadratic (3 factor)",
  heuristic = "Split",
  deltaBound=0.4535  # at this delta the support sets for points 3 and 6 begin to touch 
)

heuristics_list = list(fullsplit_lin_1factor,
                       halfsplit_lin_1factor,
                       retain_lin_1factor,
                       retain_quad_1factor,
                       split_quad_1factor,
                       split_lin_2factor, 
                       retain_lin_2factor,
                       split_quad_3factor,
                       retain_quad_3factor)

models_list = c(rep(list(lin_1factor_model),3),
                rep(list(quad_1factor_model),2),
                rep(list(lin_2factor_model),2),
                rep(list(quad_3factor_model),2))

As_list = c(rep(list(A_lin_1factor),3),
            rep(list(A_quad_1factor),2),
            rep(list(A_lin_2factor),2),
            rep(list(A_quad_3factor),2))

library(lhs)
Tmc_1factor <- 2*randomLHS(n=100, k=1)-1 
Tmc_2factor <- 2*randomLHS(n=100, k=2)-1 
Tmc_3factor <- 2*randomLHS(n=100, k=3)-1 
Tmcs_list <- c( rep(list(Tmc_1factor),5), 
                rep(list(Tmc_2factor),2),
                rep(list(Tmc_3factor),2))

Tmax_1factor <- rbind( 2*randomLHS(50,1)-1, 1, -1 )
Tmax_2factor <- rbind( 2*randomLHS(50,2)-1, fullfactorial )
Tmax_3factor <- rbind( 2*randomLHS(50,3)-1, 
                       as.matrix(
                         expand.grid(c(-1,1), 
                                     c(-1,1),
                                     c(-1,1)) 
                       ) 
)
Tmaxs_list <- c( rep(list(Tmax_1factor),5), 
                 rep(list(Tmax_2factor),2),
                 rep(list(Tmax_3factor),2))
numfactors <- rep(c(1,2,3), c(5,2,2))

# pcs <- unique(c(seq(from=0.01,to=0.1,length=50),seq(from=0.0,to=0.5,length=50) ))
tau2s <- sort(unique(c(seq(from=0.0,to=0.2,length=51), seq(from=0,to=1,length=51))))

results_df <- data.frame(NULL)
for (i in 1:length(heuristics_list)) {
  heuristic = heuristics_list[[i]]
  cat("\n\n", heuristic$model, "\n")
  cat(heuristic$heuristic,"\n")
  for (tau2 in tau2s) {
    # tau2 <- 2^numfactors[i] * 1 * pc^2
    cat("tau2 = ",tau2,".......................\r")
    results_df <- rbind(results_df, 
                          optimal_heuristic_stats(heuristic, 
                                                tau2, 
                                                models_list[[i]], 
                                                As_list[[i]], 
                                                Tmcs_list[[i]], 
                                                Tmaxs_list[[i]]))
  }
}

library(tidyverse)
results_df <- mutate(results_df, pc=sqrt(tau2/(2^numfactors)))

# write.csv(results_df, "heuristics-results.csv", row.names=FALSE)

# for plots etc see Ex-combined-plots.R

# investigating leveling off phenomenon for 3 factor quadratic
numDeltas = 500
deltas <- seq(from=0.01,0.4,length=numDeltas) # grid of delta values 
lambdaMax <- rep(NA,numDeltas) # space to store max eigenvalue for each delta
worstIdx <- rep(NA,numDeltas) # space to store index of maximizing t for each delta 
lambdas <- rep(NA,nrow(Tmax_3factor))  # space for eigenvalue for each different t (fixed delta)
mivs <- rep(NA,numDeltas)

# NB this loop takes some time 
for (i in 1:numDeltas) {
  mivs[i] <- MIV.approx.repl(xi.bar = retain_quad_3factor$meandesign(deltas[i]), 
                          r = retain_quad_3factor$r, 
                          delta = deltas[i], 
                          Tmc = Tmc_3factor, 
                          f = quad_3factor_model, 
                          A = A_quad_3factor,
                          sigma2.UB = 1)
  for (j in 1:nrow(Tmax_3factor)) {
    lambdas[j] <- (biassq.approx.repl( t=Tmax_3factor[j,],
                                       xi.bar=retain_quad_3factor$meandesign(deltas[i]), 
                                       r=retain_quad_3factor$r, 
                                      delta=deltas[i], 
                                       tau2=1, 
                                       f=quad_3factor_model, 
                                       A=A_quad_3factor) -1) * deltas[i]^3
  }
  worstIdx[i] = which.max(lambdas)
  lambdaMax[i] = max(lambdas)
}

pdf("../figures/problem4-lambdaMax.pdf", width=3,height=3)
par(mar=c(4,4,2,1)+0.1, mex=0.8, cex.axis=0.8, mgp=c(2,1,0))
plot(deltas,lambdaMax,
      type="l", 
      xlab=expression(delta), 
      ylab=expression(max[t%in%T]*lambda[max](dot(K)[d(t)])))
tChanges <- abs(diff(worstIdx))>0 # at which delta does the maximizing t change?
abline(v=deltas[tChanges][1], lty=2) # add vertical lines to the plot at such values
dev.off()