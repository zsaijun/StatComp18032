## lnorm
lnorm <- function(a, p = 1) {
  (sum(abs(a) ^ p)) ^ (1. / p)
}

# Autoregressive covariance structure
CorrAR <- function(p,rho){
  Sigma <- matrix(nrow=p,ncol=p,NA)
  for(i in 1:p){
    for(j in 1:p){
      Sigma[i,j] <- rho^(abs(i-j))
    }
  }
  Sigma
}

# Compound symmetry covariance structure
CorrCS <- function(p,rho){
  Sigma <- matrix(nrow=p,ncol=p,rho)
  diag(Sigma) <- 1
  Sigma
}

#' @title Simulation model 1
#'
#' @description similar to the the simulation model in SRRR (Lisha Chen & Jianhua Z. Huang 2012), JASS
#' @import MASS
#' @param Q parameter to specify seeds
#' @param n,p,q model dimensions
#' @param p0,q0 parameters for model sparsity
#' @param nrank model rank
#' @param s2n signal to noise ratio
#' @param rho_x correlation parameter in the generation of predictors
#' @param rho_e correlation parameter in the generation of random errors
#'
#' @return similated model and data
#' @examples
#' \dontrun{
#' data<-Simulation1(Q=1e3+10,n=50,p=20,q=10,rho_x=0.0,rho_e=0.0,s2n=1,p0=10,q0=5,nrank=3)
#' }
#' @export
Simulation1<-function(Q,n,p,q,rho_x,rho_e,s2n,p0,q0,nrank){
   set.seed(Q)
   Xsigma<-CorrAR(p,rho_x)
   X <- mvrnorm(n,rep(0,p),Xsigma)
   E <- mvrnorm(n,rep(0,q),CorrCS(q,rho_e))
   A1 <- matrix(ncol=nrank,nrow=q0,rnorm(q0*nrank))
   A0 <- matrix(ncol=nrank,nrow=q-q0,0)
   A <- rbind(A1,A0)
   B1 <- matrix(ncol=nrank,nrow=p0,rnorm(p0*nrank))
   B0 <- matrix(ncol=nrank,nrow=p-p0,0)
   B <- rbind(B1,B0)
   C <- B%*%t(A)
   svdC<-svd(C)
   C3 <- svdC$u[,nrank]%*%t(svdC$v[,nrank])*svdC$d[nrank]
   Y3 <- X%*%C3
   sigma <- sqrt(sum(as.numeric(Y3)^2)/sum(as.numeric(E)^2))/s2n
   E<- E*sigma
   Y=X%*%C+E
   list(Y=Y,X=X,C=C,Xsigma=Xsigma)
   }
