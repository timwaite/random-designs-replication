#
#  Supporting material for "Replication in random translation designs for model-robust prediction"  
# by Waite, T.W. (2023+)

tr <- function(X) { sum(diag(X)) }

Fmat <- function(xi, f) {
  t(apply(xi,1,f))
}

M <- function(xi, f) {
  F1 <- Fmat(xi,f)
  t(F1)%*%F1
}

Vobfun <- function(design,f, A) {
  M1 <- M(design,f)
  tr(solve(M1, A))
}

VobfunVectorWrapper <- function(x,f, n,q, A) {
  xi <- matrix(x, nrow=n, ncol=q) 
  Vobfun(xi, f, A)
}


d <- function(xi.bar, t) { 
  xi.bar + matrix(t, nrow=nrow(xi.bar), ncol=ncol(xi.bar), byrow=T) 
}


K <- function(xi,f,A) {
  F1 <- Fmat(xi,f)
  M1 <- t(F1)%*%F1
  MiFt <- solve( M1, t(F1) )
  t(MiFt)%*%A%*%MiFt
}


MIV.approx <- function( xi.bar, delta, Tmc, f, A, sigma2.UB) {		
  tmp <- 0 
  for (i in 1:nrow(Tmc)) {
    xi <- d(xi.bar, Tmc[i,] *delta/2 )
    tmp <- tmp + tr(solve(M(xi,f),A)) /nrow(Tmc)
  }
  sigma2.UB * tmp
}

biassq.approx <- function( t, xi.bar, delta, tau2, f, A) {
  q <- ncol(xi.bar)
  xi <- d(xi.bar, t*delta/2)
  K1 <- K(xi, f, A)
  ans <- tau2 + tau2 *max(eigen(K1,symmetric=T)$values) / (delta^q) 
  ans
}

Psi.approx <- function( xi.bar, delta,  Tmc, f, A, sigma2.UB, Tmax, tau2) {
  
  MIV1 <- MIV.approx(xi.bar, delta, Tmc, f, A, sigma2.UB)
  
  biassq1 <- rep(0,nrow(Tmax))
  for (i in 1:nrow(Tmax)) {
    biassq1[i] <- biassq.approx( Tmax[i,], xi.bar, delta, tau2, f, A)
  }
  
  MIV1 + max(biassq1)
}


Fdot <- function(xi,r,f) {
  t(apply(xi,1,f)) * r
}

Mdot <- function(xi,r,f) {
  Fm <- Fdot(xi,r,f)
  t(Fm)%*%diag(1/r)%*%Fm 
}


Kdot <- function(xi,r,f,A) {
  F1 <- Fdot(xi,r,f)
  M1 <- t(F1)%*%diag(1/r)%*%F1
  MiFt <- solve( M1, t(F1) )
  t(MiFt)%*%A%*%MiFt
}



MIV.approx.repl <- function( xi.bar, r, delta, Tmc, f, A, sigma2.UB) {		
  tmp <- 0 
  for (i in 1:nrow(Tmc)) {
    xi <- d(xi.bar, Tmc[i,] *delta/2 )
    tmp <- tmp + tr(solve(Mdot(xi,r,f),A)) /nrow(Tmc)
  }
  sigma2.UB * tmp
}


biassq.approx.repl <- function( t, xi.bar, r, delta, tau2, f, A) {
  q <- ncol(xi.bar)
  xi <- d(xi.bar, t*delta/2)
  K1 <- Kdot(xi, r, f, A)
  ans <- tau2 + tau2 *max(eigen(K1,symmetric=T)$values) / (delta^q) 
  ans
}

Psi.approx.repl <- function( xi.bar, r, delta,  Tmc, f, A, sigma2.UB, Tmax, tau2) {
  
  MIV1 <- MIV.approx.repl(xi.bar, r, delta, Tmc, f, A, sigma2.UB)
  
  biassq1 <- rep(0,nrow(Tmax))
  for (i in 1:nrow(Tmax)) {
    biassq1[i] <- biassq.approx.repl( Tmax[i,], xi.bar, r, delta, tau2, f, A)
  }
  
  MIV1 + max(biassq1)
}
