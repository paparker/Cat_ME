library(Matrix)
library(BayesLogit)
library(mvtnorm)
library(sampling)
library(BayesLogit)
library(foreach)
library(doParallel)
library(gtools)


pgBin <- function(x, y, iter=100, burn=50, wgt=NULL){
  n <- nrow(x)
  p <- ncol(x)
  if(is.null(wgt)) wgt <- rep(1,n)
  beta <- rep(0,p)
  betaOut <- matrix(NA, nrow=iter, ncol=p)
  Binv <- Diagonal(p, 0.1)
  W <- rep(1, n)
  K <- wgt*(y-0.5)
  pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    precBeta <- solve(t(x)%*%Diagonal(length(W),W)%*%x + Binv)
    meanBeta <- t(x)%*%Diagonal(length(W),W)%*%(K/W)
    beta <- betaOut[i,]  <- as.numeric(rmvnorm(1, mean=precBeta%*%meanBeta, sigma=as.matrix(precBeta)))
    
    W <- rpg.gamma(n, wgt, x%*%beta)
    setTxtProgressBar(pb, i)
  }
  return(list(Beta=betaOut[-c(1:burn),]))
}


regME2 <- function(x1, x2, yME, y, iter=100, burn=50, wgt2=NULL,  a1=1, b1=1, a2=1, b2=1, ids){
  x1 <- x1[ids, ]
  y <- y[ids]
  wgt1 <- wgt2[ids]
  wgt1 <- wgt1/mean(wgt1)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  p1 <- ncol(x1)+1
  p2 <- ncol(x2)
  if(is.null(wgt1)) wgt1 <- rep(1,n1)
  if(is.null(wgt2)) wgt2 <- rep(1,n2)
  beta1 <- rep(0, p1)
  beta2 <- rep(0, p2)
  beta1Out <- matrix(NA, nrow=iter, ncol=p1)
  beta2Out <- matrix(NA, nrow=iter, ncol=p2)
  B1inv <- Diagonal(p1, 0.1)
  B2inv <- Diagonal(p2, 0.1)
  W1 <-  rep(1, n1)
  W2 <-  rep(1, n2)
  yT <- yME
  K1 <- wgt1*(y-0.5)
  K2 <- wgt2*(yT-0.5)
  s1 <- 0.95
  s1Out <- rep(NA, iter)
  s2 <- 0.95
  s2Out <- rep(NA, iter)
  yTout <- matrix(NA, nrow=n2, ncol=iter)
  pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    
    ### Top level regression
    X1 <- cbind(x1, yT[ids])
    precBeta <- solve(t(X1)%*%Diagonal(length(W1),W1)%*%X1 + B1inv)
    meanBeta <- t(X1)%*%Diagonal(length(W1),W1)%*%(K1/W1)
    beta1 <- beta1Out[i,]  <- as.numeric(rmvnorm(1, mean=precBeta%*%meanBeta, sigma=as.matrix(precBeta)))
    
    W1 <- rpg.gamma(n1, wgt1, X1%*%beta1)
    
    
    ### Bottom level regression
    precBeta <- solve(t(x2)%*%Diagonal(length(W2),W2)%*%x2 + B2inv)
    meanBeta <- t(x2)%*%Diagonal(length(W2),W2)%*%(K2/W2)
    beta2 <- beta2Out[i,]  <- as.numeric(rmvnorm(1, mean=precBeta%*%meanBeta, sigma=as.matrix(precBeta)))
    
    ### Latent rvs
    W2 <- rpg.gamma(n2, wgt2, x2%*%beta2)
    pll <- plogis(x2%*%beta2)
    p11 <- (pll)^(wgt2)*s2^(wgt2*yME)*(1-s2)^(wgt2*(1-yME))
    p11[ids] <- p11[ids]*dbinom(y, 1, plogis(x1%*%beta1[-c(p1)] + beta1[c(p1)]))
    p00 <- (1-pll)^(wgt2)*(1-s1)^(wgt2*yME)*(s1^(wgt2*(1-yME)))
    p00[ids] <- p00[ids]*dbinom(y, 1, plogis(x1%*%beta1[-c(p1)]))
    yT <- rbinom(n2, 1, p11/(p11+p00))
    yTout[,i] <- yT
    K2 <- wgt2*(yT-0.5)
    
    
    ### Error rates
    s1 <- s1Out[i] <- rbeta(1, a1+sum((1-yME)*(1-yT)*wgt2), b1+ sum(yME*(1-yT)*wgt2))
    s2 <- s2Out[i] <-  rbeta(1, a2+sum((yME)*(yT)*wgt2), b2+ sum((1-yME)*(yT)*wgt2))

    setTxtProgressBar(pb, i)
  }
  return(list(Beta1=beta1Out[-c(1:burn),], Beta2=beta2Out[-c(1:burn),], sp=s2Out[-c(1:burn)], se=s1Out[-c(1:burn)], yT=yTout[,-c(1:burn)]))
}
