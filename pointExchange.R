#
#  Supporting material for "Replication in random translation designs for model-robust prediction" 
# by Waite, T.W. (2023+)


pointExchange <- function(initialDesign, candidateSet, objfun, tol=1e-6, passMax=100, ...) {
  # initialDesign here should be a design matrix 
  # objfun takes a design matrix as input and returns the objective function value 
  # - NB this is different than coord.descent, where the function must take a vector as input, not a matrix 
  
  bestObjective <- Inf 
  stopFlag <- F
  passCtr <- 1
  
  design <- initialDesign
  designTmp <- initialDesign
  
  while (!stopFlag & passCtr <= passMax) {
    
    prepassObjective <- objfun(design,...)
    
    for (i in 1:nrow(design)) {
      for (j in 1:nrow(candidateSet)) {
        
        designTmp[i, ] <- candidateSet[j,]
        objTmp <- objfun(designTmp, ...)
        
        if (objTmp < bestObjective) {
          bestObjective <- objTmp
          design <- designTmp
        } else {
          designTmp <- design
        }
      }
    }
    
    postpassObjective <- objfun(design,...)
    if (prepassObjective-postpassObjective < tol) stopFlag=T
    passCtr <- passCtr + 1 
  }
  
  return(list(design=design, objective=bestObjective))
}



# to do: refactoring of Vobfun, separate vector-to-matrix conversion from obj fun evaluation so can use same obj fun in pointExchange
# incorporate random starts 


