#
#  Supporting material for "Replication in random translation designs for model-robust prediction" 
# by Waite, T.W. (2023+)


# Section 4: Numerical examples 
## Section 4.1: Example 1 - 6 points, 1 factors, quadratic model

source("obj-funs-replication.R")

# set up model
f=function(x) { c(1,x,x^2) }
A=diag(c(4,4/3,4/5)); A[1,3] <- 4/3 -> A[3,1]


# V-optimal design 

des <- c(-1,-1,0,0,1,1)

heur.des.split <- function(delta) {
  matrix(c(-1+delta/2, -1+3*delta/2, -delta/2, delta/2, 1-3*delta/2, 1-delta/2), ncol=1)
}

heur.des.repl <- function(delta) {
  matrix(c(-1+delta/2, 0, 1-delta/2), ncol=1)
}

r <- c(2,2,2)


library(lhs)
Tmc <- 2*randomLHS(n=100, k=1)-1 
Tmax <- rbind( 2*randomLHS(50,1)-1, 1, -1 )


pc <- seq(from=0.0,to=0.5,length=50) 
tau2s <- 2^1 * 1 * pc^2


minPsi.split <- rep(NA, length(tau2s))
optDelta.split <- rep(NA, length(tau2s))
MIV.split <- rep(NA,length(tau2s))
MISB.split <- rep(NA,length(tau2s))

minPsi.repl <- rep(NA,length(tau2s))
optDelta.repl <- rep(NA, length(tau2s))
MIV.repl <- rep(NA,length(tau2s))
MISB.repl <- rep(NA,length(tau2s))


for (i in 1:length(tau2s)) {
  #
  # Risk bound for heuristic strategy, as a function of delta
  # - optimize delta and plot the risk bound
  #
  cat("Case ",i,"/",length(tau2s),"\r")
  tau2 <- tau2s[i]
  
  PsiH.repl <- function(delta) { Psi.approx.repl( xi.bar=heur.des.repl(delta), r, delta, Tmc, f, A, sigma2.UB=1, Tmax, tau2=tau2 ) } 
  optH.repl <- optimize(Vectorize(PsiH.repl), lower=0.0001, upper=2/3)  
  minPsi.repl[i] <- optH.repl$objective
  optDelta.repl[i] <- optH.repl$minimum
  MIV.repl[i] <- MIV.approx.repl(xi.bar=heur.des.repl(optDelta.repl[i]), r, optDelta.repl[i], Tmc, f, A, 1)
  MISB.repl[i] <- minPsi.repl[i] - MIV.repl[i]
  
  PsiH.split <- function(delta) { Psi.approx( xi.bar=heur.des.split(delta), delta, Tmc, f, A, sigma2.UB=1, Tmax, tau2=tau2 ) } 
  optH.split <- optimize(Vectorize(PsiH.split), lower=0.0001, upper=2/6)  
  optH.split   # best objective 4.8577 - better to split than to replicate
  minPsi.split[i] <- optH.split$objective
  optDelta.split[i] <- optH.split$minimum
  MIV.split[i] <- MIV.approx(xi.bar=heur.des.split(optDelta.split[i]),   optDelta.split[i], Tmc, f, A, 1)
  MISB.split[i] <-   minPsi.split[i] - MIV.split[i]
}

results_1factor <- cbind(pc, optDelta.split, optDelta.repl, minPsi.split, minPsi.repl, eff=minPsi.repl/minPsi.split)
# write.csv(results_1factor, file="results_1factor.csv")


par(mar=c(5,5,4,2)+0.1, family="serif")
xlab = "Relative std dev of discrepancy (randomly selected point)"

pdf(file = "../figures/optDelta-1factor.pdf", width=5, height=5)
plot(pc, optDelta.repl, type="l", 
     xlab=xlab, 
     ylab="Optimal value of delta", 
     sub="(1-factor example)",
     lwd=2)
lines(pc, optDelta.split, col="grey", lwd=2)
legend("bottomright", legend=c("Retain replicates", "Split replicates"), col=c(1,"grey"), lty=c(1,1))
dev.off()


pdf(file = "../figures/efficiency-plot-1factor.pdf", width=5, height=5)
plot(pc, minPsi.repl/minPsi.split,type="l", 
     xlab=xlab, 
     ylab="Efficiency (splitting vs replication)",
     main="Splitting is more efficient than replication",
     sub="(1-factor example)")
dev.off()

pdf(file = "../figures/1factor-reduced-bias.pdf", width=5, height=5)
plot(pc, MIV.repl,type="l", ylim=c(0,4), 
     lwd=2, col=1, ylab="", 
     xlab=xlab,
     main="Splitting is better due to a reduced bias term",
     sub="(1-factor example)")
lines(pc, MIV.split, lwd=2, col="grey")
lines(pc, MISB.repl, lwd=2,  col=1, lty=2)
lines(pc, MISB.split, lwd=2, col="grey",lty=2)
legend("topleft", 
       legend=c("MIV (retain)", "MIV (split)", "max. MISB (retain)", "max. MISB (split)"),
       col=c(1, "grey", 1, "grey"),
       lty=c(1,1,2,2), 
       lwd=2)
dev.off()

pdf(file = "../figures/1factor-principal-eigenvalue.pdf", width=5, height=5)
par(mar=c(5,5,4,2)+0.1)
plot(pc, (MISB.split/tau2s-1)*optDelta.split,type="l", col="grey", ylim=c(0,3), 
     main="Principal eigenvalue reduced by splitting",
     sub="(1-factor example)",
     xlab=xlab,
     ylab=expression(max[t %in% T ]*lambda[max](dot(K)[d(t)])))
lines(pc, (MISB.repl/tau2s-1)*optDelta.repl, col=1)
legend("topright", c("Retain replicates", "Split replicates"), lty=c(1,1), col=c(1,"grey"))
dev.off()
