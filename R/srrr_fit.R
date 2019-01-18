
#' @title Fitting Sparse Reduced-Rank Regression
#'
#' @description Given a response matrix and a covariate matrix, fits sparse reduced-Rank
#' regression for a specified rank.
#' @import remMap
#' @param Y a matrix of response (n by q)
#' @param X a matrix of covariate (n by p)
#' @param lam penalty parameter for B
#' @param nrank an integer specifying the desired rank
#' @param epsilon convergence tolerance
#' @param maxit maximum number of iterations
#' @return a list of model parameters
#' @examples
#' \dontrun{
#' data<-Simulation1(Q=1e3+10,n=50,p=20,q=10,rho_x=0.0,rho_e=0.0,s2n=1,p0=10,q0=5,nrank=3)
#' Result<-srrr.fit(Y=data$Y,X=data$X,lam=1,nrank=3,epsilon=1e-3,maxit=1e4)
#' }
#' @export
srrr.fit <- function(Y,X,lam,nrank,epsilon,maxit)
{
  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Y);r<-nrank;lam<-lam*sqrt(log(p)/n)
  B<-matrix(0.1,p,r);A<-matrix(0.1,q,r)
  VALUE<-rep(0,maxit)
  for (i in 1:maxit) { ##    while(diff[iter] > conv & iter < maxit){
    ## 1. B step
    B0<-B
    B<-remMap(X, Y%*%A,lamL1=0, lamL2=lam, C.m=NULL)$phi

    ## 2. A step
    A0<-A
    Z1<-t(Y)%*%X%*%B
    Z1svd<-svd(Z1)
    A<-Z1svd$u%*%t(Z1svd$v)

    a<-0
    for(j in 1:(p)){a<-a+norm(B[j,],type="2")}
    VALUE[i]<-norm(Y-X%*%B%*%t(A),type="F")^2+lam*a
    norm(Y-X%*%B%*%t(A),type="F")^2+lam*a
    # del <- sqrt(sum((B0-B)^2) + sum((B_0-B_)^2)+sum((A0-A)^2)+sum((A_0-A_)^2))
    if(i>1) {if (abs(VALUE[i]-VALUE[i-1]) < epsilon) break
    }
  }
  list(B=B,A=A,step=i)
}
