#
#  Supporting material for "Replication in random translation designs for model-robust prediction" 
# by Waite, T.W (2023+)

# Section 4 - Numerical examples
# Test Problem 4 - 3 factors, full quadratic model
# Script to find traditional V-optimal design, on which the random translation strategies are based 
#  saves result to design_3factor.csv

source("obj-funs-replication.R")
source("coord-descent.R")

f <- function(x) { c(1,x[1],x[2],x[3], x[2]*x[3], x[1]*x[3], x[1]*x[2], x[1]^2, x[2]^2, x[3]^2)}

#
# Caclulate A matrix
#
A <- matrix(0, ncol=10, nrow=10)

# represent model terms via their exponents 
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

# compute integrals of x1^idx[1], x2^idx[2], x3^idx[3]
INT <- function(idx) { 2/(idx+1) * ( (idx+1) %% 2 == 1) }

# Aij = product int f[i] f[j] dx
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

# NB also tried a more thorough search with point exchange, obtained a design wit the same objective function value to 4dp

