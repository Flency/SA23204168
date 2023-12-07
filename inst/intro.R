## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----warning=FALSE------------------------------------------------------------
library(MASS)
library(SA23204168)
set.seed(1)
n=100
p=200
m=20
true.beta<-c(-0.6,-0.6,-0.6,rep(0,p-3))
Mu <- rep(0,p) 
Sigma <- diag(1,p)
X<-mvrnorm(n, Mu, Sigma)
y1 <- (1+X %*% true.beta)^3+sin(X%*%true.beta)+rnorm(n)
y2 <- X%*%true.beta+rnorm(n)
y<-cbind(y1,y2)
covX<-cov(X)
eigencovX <- eigen(covX)
sqcovX <- eigencovX$vectors%*%sqrt(diag(pmax(eigencovX$values,0)))%*%t(eigencovX$vectors)
initPi<-true.beta%*%t(true.beta)
initH<-sqcovX%*%initPi%*%sqcovX
beta<-sdr(X=X,y=y,m=m,lambda=10,K=1,nu=0.01,alpha=0.01,epsilon=1e-2,maxiter=10,trace=FALSE,initPi=initPi,initH=initH)

print(beta[1:6])

## -----------------------------------------------------------------------------
x <- matrix(c(1, 2, 4, 9), 2, 2)
print(lmax(x))

