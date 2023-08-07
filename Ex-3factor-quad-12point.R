#
#  Supporting material for "Replication in random translation designs for model-robust prediction" 
# by Waite, T.W (2023+)


# Section 4 - Numerical examples

# Section 4.2: Example 2 - 12 points, 3 factors, quadratic model

source("obj-funs-replication.R")
source("coord-descent.R")

f <- function(x) { c(1,x[1],x[2],x[3], x[2]*x[3], x[1]*x[3], x[1]*x[2], x[1]^2, x[2]^2, x[3]^2)}

#
# Caclulate A matrix
#
A <- matrix(0, ncol=10, nrow=10)

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
    A[i,j] <- prod( sapply(indxs[i,] + indxs[j,], INT) )
  }
}

library(MASS); fractions(A)

####
#
# V-optimal deterministic design 
#
#

library(lhs)
init <- c(2*randomLHS(n=12,k=3)-1)
ansV <- coord.descent(init, VobfunVectorWrapper, f=f, n=12, q=3, A=A)   # best value 3.9178, worth trying several repeats
des <- data.frame(matrix(ansV$xcurr,nrow=12,ncol=3))
des <- round(des,digits=3) 
des <- des[do.call(order,des),]

# write.csv(des, file="design_3factor.csv")
des = read.csv("design_3factor.csv")[,2:4]
## also tried a more thorough search incorporating d multiple random starts and consolidation of replicates via point exchange 
## gives a similar design with replication and the same objective function value to 4dp
#
# randomStarts <- 100
# bestObjective <- Inf
# des <- NULL
# for (startIdx in 1:randomStarts) { 
#   cat("\rRandom start ", startIdx, "/", randomStarts)
#   init <- c(2*randomLHS(n=12,k=3)-1)
#   desVector <- coord.descent(init, VobfunVectorWrapper, verbose=F, f=f, n=12, q=3, A=A)    
#   desMatrix <- data.frame(matrix(ansV$xcurr,nrow=12,ncol=3))
#   peResult <- pointExchange(desMatrix, desMatrix, Vobfun, f=f, A=A)
#   if (peResult$objective < bestObjective) {
#     des <- peResult$design
#   }
#   flush.console()
# }
# des <- des[do.call(order,des),]  

#
# Set up heuristic xi.bar as described in text 
# - N.B. - first check that rows 6 & 7 of des are replicates
#

des_unique <- unique(des)
r <- rep(1,11); r[6] <- 2
cbind(des_unique, r)

clip <- function(x, delta=0.1) {
  x <- as.matrix(x)
  x<-  pmin(  x, 1-delta/2 )
  x <- pmax( x, -1+delta/2)
  return(x)
}

heur.des.repl <- function(delta) {
  xi <- as.matrix(des_unique)
  xi <- clip(xi, delta)  
  return(xi)
}


heur.des.split <- function(delta) {
  xi <- as.matrix(des)
  xi <- clip(xi, delta)  
  return(xi)
}


Tmc <- 2*randomLHS(n=100, k=3)-1 
Tmax <- rbind( 2*randomLHS(50,3)-1, as.matrix( expand.grid(c(-1,1), c(-1,1),c(-1,1)) ) )
# NB inclusion of the corner points in Tmax stabilizes the results as the worst t is usually at one of these points

# tau2=2^3*1*0.01^2 # sd of psi at randomly selected point is **1**% that of the random error sd

pc <- seq(from=0.0,to=0.1,length=50) 
tau2s <- 2^3 * 1 * pc^2


minPsi.split <- rep(NA, length(tau2s))
optDelta.split <- rep(NA, length(tau2s))
minPsi.repl <- rep(NA,length(tau2s))
optDelta.repl <- rep(NA, length(tau2s))

for (i in 1:length(tau2s)) {
  #
  # Risk bound for heuristic strategy, as a function of delta
  # - optimize delta and plot the risk bound
  #
  cat("Case ",i,"/",length(tau2s),"\r")
  tau2 <- tau2s[i]
  
  PsiH.repl <- function(delta) { Psi.approx.repl( xi.bar=heur.des.repl(delta), r, delta, Tmc, f, A, sigma2.UB=1, Tmax, tau2=tau2 ) } 
  optH.repl <- optimize(Vectorize(PsiH.repl), lower=0.01, upper=0.5)  
  minPsi.repl[i] <- optH.repl$objective
  optDelta.repl[i] <- optH.repl$minimum
  
  PsiH.split <- function(delta) { Psi.approx( xi.bar=heur.des.split(delta), delta, Tmc, f, A, sigma2.UB=1, Tmax, tau2=tau2 ) } 
  optH.split <- optimize(Vectorize(PsiH.split), lower=0.01, upper=0.5)  
  optH.split   # best objective 4.8577 - better to split than to replicate
  minPsi.split[i] <- optH.split$objective
  optDelta.split[i] <- optH.split$minimum
}

pdf("../figures/optDelta-3factor.pdf", width=5, height=5)
plot(pc, optDelta.split, type="l", col="grey", 
     ylab="Optimal delta", 
     xlab="Relative std dev of discrepancy (randomly selected point)",
     sub="(3-factor example)")
lines(pc, optDelta.repl, col=1)
abline(h=0.2780762, lty=2)
legend("bottomright", c("Retain replicates", "delta = 0.278", "Split replicates"), lty=c(1,2,1), col=c(1,1,"grey"))
dev.off()

pdf("../figures/efficiency-plot-3factor.pdf", width=5, height=5)
plot(pc, minPsi.repl/minPsi.split,type="l", 
     xlab="Relative std dev of discrepancy (randomly selected point)", 
     ylab="Efficiency (splitting vs replication)",
     main="Splitting is more efficient than replication",
     sub="(3-factor example)")
dev.off()

results_3factor <- cbind(pc, optDelta.split, optDelta.repl, minPsi.split, minPsi.repl, eff=minPsi.repl/minPsi.split)
# write.csv(results_3factor, file="results_3factor.csv")


# investigating unusual leveling off in graphs of optimal delta

# plot figure of maximal principal eigenvalue of bias-sensitivity matrix for varying delta 
# also identify the maximizing t

numDeltas = 500
deltas <- seq(from=0.01,to=0.40,length=numDeltas) # grid of delta values 

Tmax <- rbind( 2*randomLHS(50,3)-1, as.matrix( expand.grid(c(-1,1), c(-1,1),c(-1,1)) ) )

lambdaMax <- rep(NA,numDeltas) # space to store max eigenvalue for each delta
worstIdx <- rep(NA,numDeltas) # space to store index of maximizing t for each delta 

lambdas <- rep(NA,nrow(Tmax))  # space for eigenvalue for each different t (fixed delta)

# NB this loop takes some time 
for (i in 1:numDeltas) {
  for (j in 1:nrow(Tmax)) {
    lambdas[j] <- (biassq.approx.repl( Tmax[j,], heur.des.repl(deltas[i]), r=r, deltas[i], 1, f, A) -1) * deltas[i]^3
  }
  worstIdx[i] = which.max(lambdas)
  lambdaMax[i] = max(lambdas)
}
 
pdf("../figures/3factor-repl-eigenvalue.pdf", width=5, height=5)
par(mar=c(5,5,4,2)+0.1)
plot(deltas,lambdaMax,
     type="l", 
     xlab=expression(delta), 
     ylab=expression(max[t%in%T]*lambda[max](dot(K)[d(t)])))
tChanges <- abs(diff(worstIdx))>0 # at which delta does the maximizing t change?
abline(v=deltas[tChanges][1], lty=2) # add vertical lines to the plot at such values
legend("topleft", legend=c("delta = 0.278"), lty=2)
dev.off()
                    