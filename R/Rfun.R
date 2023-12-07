#' @import boot
#' @import bootstrap
#' @import DAAG
#' @import coda
#' @import microbenchmark
#' @import MASS
#' @importFrom Rcpp evalCpp
#' @importFrom stats rbeta rbinom
#' @useDynLib SA23204168
NULL

#' @title A Gibbs sampler using R
#' @description A Gibbs sampler using R
#' @param N the number of samples
#' @param a a param of beta distribution
#' @param b another param of beta distribution
#' @param n sample size
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' rnR <- gibbsR(100,2,3,10)
#' par(mfrow=c(2,1));
#' plot(rnR[,1],type='l')
#' plot(rnR[,2],type='l')
#' }
#' @export
gibbsR <- function(N,a,b,n) {
  mat <- matrix(nrow = N, ncol = 2)
  x <- y <- 0
  for (i in 1:N) {
    x <- rbinom(1, n, y)
    y <- rbeta(1, x+a, n-x+b)
    mat[i, ] <- c(x, y)
  }
  mat
}

######################################################
# Update Pi
diff_f<-function(X_diff,y_matrix,Pi){
  M<-X_diff%*%Pi%*%t(X_diff)
  nonzero_index<-which(diag(M)!=0)
  X_diff<-X_diff[nonzero_index,]
  y_matrix<-y_matrix[nonzero_index,nonzero_index]
  M<-X_diff%*%Pi%*%t(X_diff)
  negative_half_power_matrix<-diag(diag(M)^(-1/2),ncol=ncol(M),nrow=nrow(M))
  Lambda<-negative_half_power_matrix%*%y_matrix/(2*nrow(X_diff))
  f<-t(X_diff)%*%Lambda%*%X_diff
  return(f)
}

######################################################
updatePi <- function(f,covx,sqcovx,H,Gamma,U,Lambda,alpha,nu,lambda,Pi,tau,delta){
  A <- (2*delta+tau)*Pi+f+alpha*(U-Lambda)- nu*covx%*%Pi%*%covx+nu*sqcovx%*%(H-Gamma)%*%sqcovx
  B <- lambda
  temp<-1/(2*delta+alpha+tau)*Soft(A,B)
  return(temp)
}

######################################################
# Update H
######################################################
updateH <- function(sqcovx,Gamma,nu,Pi,K){
  
  temp <- Gamma + sqcovx%*%Pi%*%sqcovx
  temp <- (temp+t(temp))/2
  svdtemp <- eigen(temp)
  d <- svdtemp$values
  p <- length(d)
  
  if(sum(pmin(1,pmax(d,0)))<=K){
    dfinal <- pmin(1,pmax(d,0))
    return(svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors))
  }
  
  fr <- function(x){
    sum(pmin(1,pmax(d-x,0)))
  }
  # Vincent Vu Fantope Projection
  knots <- unique(c((d-1),d))
  knots <- sort(knots,decreasing=TRUE)
  temp <- which(sapply(knots,fr)<=K)
  lentemp <- tail(temp,1)
  a=knots[lentemp]
  b=knots[lentemp+1]
  fa <- sum(pmin(pmax(d-a,0),1))
  fb <- sum(pmin(pmax(d-b,0),1))
  theta <- a+ (b-a)*(K-fa)/(fb-fa)
  dfinal <- pmin(1,pmax(d-theta,0))
  res <- svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors)
  return(res)
}


make_minibatch<-function(X,y,m){
  n<-nrow(X)
  shuffled_vector <- 1:n
  nl <- length(shuffled_vector) %/% m
  split_vector <- split(shuffled_vector, f = rep(1:m, each = nl, length.out = length(shuffled_vector)))
  split_matrix<-matrix(unlist(split_vector),nc=m,byrow = FALSE)
  # minibatch
  x_array<-array(0, dim=c(nl*(nl-1),p,m))
  y_array<-array(0, dim=c(nl*(nl-1),nl*(nl-1),m))
  covx_array<-array(0,dim=c(p,p,m))
  sqcovx_array<-array(0,dim=c(p,p,m))
  max_eigen_value<-rep(0,m)
  
  for(k in 1:m){
    #minibatch
    combinations <- expand.grid(split_matrix[,k],split_matrix[,k])  # 生成两两组合
    selected<-(combinations$Var1!=combinations$Var2)
    index<-combinations[selected,]
    X_pre<-X[index[,1],]-X[index[,2],]
    x_array[,,k]=X_pre
    y_diff<-y[index[,1],]-y[index[,2],]
    y_diff<-apply(y_diff,1, function(row) norm(row, type = "2"))
    y_diff_mean<-rep(mean(y_diff),length(y_diff))
    y_diff_byk<-matrix(y_diff,length(y_diff)/nl,nl,byrow = FALSE)
    y_diff_kmean<-rep(colMeans(y_diff_byk),each=length(y_diff)/nl)
    y_combined<-y_diff+y_diff_mean-2*y_diff_kmean
    y_pre<-diag(y_combined,nrow=length(y_combined),ncol=length(y_combined))
    y_array[,,k]=y_pre
    covx_array[,,k]=cov(X[split_matrix[,k],])
    eigencovx <- eigen(covx_array[,,k])
    sqcovx_array[,,k] <- eigencovx$vectors%*%sqrt(diag(pmax(eigencovx$values,0)))%*%t(eigencovx$vectors)
    max_eigen_value[k]<-eigencovx$values[1]
  }
  return(list(x_array=x_array,y_array=y_array,covx_array=covx_array,sqcovx_array=sqcovx_array,max_eigen_value=max_eigen_value))
}

######################################################
# Soft-thresholding Operator
######################################################
Soft <- function(a,b){
  if(b<0) stop("Can soft-threshold by a nonnegative quantity only.")
  return(sign(a)*pmax(0,abs(a)-b))
}

#' @title sparse sufficient dimension reduction
#' @description a function estimates the sparse sufficient subspace
#' @param X the n by p predictor, p>n
#' @param y the n by q response
#' @param m minibatct size
#' @param lambda tuning parameter
#' @param K Structural dimension
#' @param initPi initial value
#' @param initH initial value
#' @return p by K basis matrix of sufficient dimension reduction subspace 
#' @examples
#' \dontrun{
#' rnR <- gibbsR(100,2,3,10)
#' par(mfrow=c(2,1));
#' plot(rnR[,1],type='l')
#' plot(rnR[,2],type='l')
#' }
#' @export
sdr <- function(X,y,m,lambda,K,nu=0.01,alpha=0.01,epsilon=1e-2,maxiter=10,trace=FALSE,initPi=NULL,initH=NULL){
  minibatch<-make_minibatch(X,y,m)
  x_array<-minibatch$x_array
  y_array<-minibatch$y_array
  covx_array<-minibatch$covx_array
  sqcovx_array<-minibatch$sqcovx_array
  max_eigen_value<-minibatch$max_eigen_value
  
  criteria <- 1e10
  i <- 1
  
  # Initialize parameters
  U<-oldU <-Pi<- oldPi<-initPi
  H<-initH
  Lambda<-Gamma <- matrix(0,p,p)
  
  time_computation<-time_update<-0
  # While loop for the iterations
  while(criteria > epsilon && i <= maxiter){
    start<-Sys.time()
    new_i<-(i-1)%%m+1
    tau <- 4*nu*max_eigen_value[new_i]^2
    #delta<-2*tau
    nabla_f<-diff_f(x_array[,,new_i],y_array[,,new_i],U)
    eigen_f<-svd(-nabla_f)
    delta<-1/2*max(eigen_f$d)
    Pi <- updatePi(nabla_f,covx_array[,,new_i],sqcovx_array[,,new_i],H,Gamma,U,Lambda,alpha,nu,lambda,Pi,tau,delta)
    eigen_Pi<-svd(Pi+Lambda)
    temp_beta<-eigen_Pi$u[,1]
    U<-eigen_Pi$d[1]*temp_beta%*%t(temp_beta)
    end<-Sys.time()
    time_computation<-time_computation+difftime(end, start, units = "min")
    
    start<-Sys.time()
    
    H <- updateH(sqcovx_array[,,new_i],Gamma,nu,Pi,K)
    Lambda<-Lambda+Pi-U
    Gamma <- Gamma + sqcovx_array[,,new_i]%*%Pi%*%sqcovx_array[,,new_i]-H
    criteria <- sqrt(sum((Pi-oldPi)^2))
    oldU <- U
    oldPi<-Pi
    i <- i+1
    if(trace==TRUE)
    {
      print(i)
      print(criteria)
    }
    
    end<-Sys.time()
    time_update<-time_update+difftime(end, start, units = "min")
  }
  final_beta<-svd(Pi)$u[,1:K]
  return(final_beta)
}
