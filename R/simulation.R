library(Matrix)
library(BayesLogit)
library(mvtnorm)
library(sampling)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(TruncatedDistributions)
library(LaplacesDemon)
library(pROC)
library(yardstick)
library(readr)
source("R/functions.R")


sim <- 100
expSS <- 2000
###############################
### Informative Sample Case ###
###############################


set.seed(1)
mod1B <- mod2B <- matrix(NA, nrow=3, ncol=sim)
mod2sp <- mod2se <- rep(NA, sim)
ESSmod2 <- compTime <-  rep(NA, sim)
for(s in 1:sim){
  ss <- 100000
  x <- cbind(rnorm(ss), rnorm(ss), rnorm(ss), rnorm(ss))
  beta <- c(0.7,-2, 0.5, -0.3)
  pT <- plogis(x%*%beta)
  yT <- rbinom(ss, 1, pT)
  yC <- yT
  cR <- 0.025
  cR2 <- 0.25
  ru <- runif(ss)
  yC[yT==0] <- 0 + (ru[yT==0] < cR)
  yC[yT==1] <- 1 - (ru[yT==1] < cR2)
  
  cv1 <- rnorm(ss)
  x2T <- cbind(1, cv1,  as.numeric(yT==1))
  x2C <- cbind(1, cv1,  as.numeric(yC==1))
  beta2 <- c(0.2, -0.7, 0.6)
  pT2 <- plogis(x2T%*%beta2)
  resp <- rbinom(ss, 1, pT2)
  
  pik <- inclusionprobabilities(exp(rnorm(ss)-resp)^0.1, expSS)
  ind <- UPpoisson(pik)
  
  xSampTT <- x2T[as.logical(ind),]
  xSampTC <- x2C[as.logical(ind),]
  xSampB <- x[as.logical(ind),]
  yCsamp <- yC[as.logical(ind)]
  respSamp <- resp[as.logical(ind)]
  wgtSamp <- 1/pik[as.logical(ind)]
  wgtSamp <- wgtSamp/sum(wgtSamp)*sum(ind)
  
  mod <- sampler(xSampTC, respSamp, iter=2500, burn=500, wgt=wgtSamp)
  mod1B[,s] <- colMeans(mod$Beta)
  
  t1 <- Sys.time()
  mod2 <- regME2(x1=xSampTC[,-c(3)], x2=xSampB, yME=yCsamp, y=respSamp, iter=2500, burn=500, wgt=wgtSamp)
  compTime[s] <- as.numeric(Sys.time() - t1)
  mod2B[,s] <- colMeans(mod2$Beta1)
  mod2sp[s] <- mean(mod2$sp)
  mod2se[s] <- mean(mod2$se)
  ESSmod2[s] <- ESS(mod2$Beta1[,3])
   
  print(s)
}

rowMeans((mod1B - beta2)^2)
rowMeans((mod2B - beta2)^2)

abs(rowMeans(mod1B) - beta2)
abs(rowMeans(mod2B) - beta2)



coefs <- factor(paste0("b",0:2), labels=c(expression(beta[0]), expression(beta[1]), expression(beta[2])))
mod1DF <- data.frame(t(mod1B))
names(mod1DF) <- coefs
mod1DF$Mod <- "Naive"

mod2DF <- data.frame(t(mod2B))
names(mod2DF) <- coefs
mod2DF$Mod <- "Mixture"


modDF <- rbind(mod1DF, mod2DF) %>% pivot_longer(1:3)


dfNew <- data.frame(name=coefs, value=c(0.2, -0.7, 0.6))

ggplot(modDF %>% filter(name == "beta[2]"), aes(x=Mod, y=value))+
  geom_boxplot(fill='green')+
  geom_hline(data=dfNew %>% filter(name == "beta[2]"), color="red", aes(yintercept=value))+
  facet_wrap(~name, scales='free', labeller = label_parsed)+
  xlab("")+
  ylab("Estimate")+
  theme_classic()+
  ggtitle("Informative Sample")


# Compute ROC values
rocDF <- data.frame(truth=factor(yT[as.logical(ind)], levels=c(0,1)), estimate=rowMeans(mod2$yT), estimate2=yC[as.logical(ind)])
auc1 <- 1-roc_auc(rocDF, truth = truth, estimate) %>% pull(.estimate)
auc2 <- 1-roc_auc(rocDF, truth = truth, estimate2) %>% pull(.estimate)

label1 <- paste0("Mixture (AUC = ", round(auc1, 3), ")")
label2 <- paste0("Naive (AUC = ", round(auc2, 3), ")")

rocDF <- data.frame(truth=factor(yT[as.logical(ind)], levels=c(0,1)), estimate=rowMeans(mod2$yT), estimate2=yC[as.logical(ind)])
roc1 <- roc_curve(rocDF, truth = truth, estimate) %>% mutate(Model = label1)
roc2 <- roc_curve(rocDF, truth = truth, estimate2) %>% mutate(Model = label2)
roc_all <- bind_rows(roc1, roc2)

ggplot(roc_all, aes(x = specificity, y = 1-sensitivity, color=Model)) +
  geom_line(size = 1) +
  geom_abline(lty = 2, color = "gray") +
  labs(x = "False Positive Rate",
       y = "True Positive Rate") +
  theme_minimal()




#################################
### Simple Random Sample Case ###
#################################



set.seed(1)
mod1B <- mod2B <- matrix(NA, nrow=3, ncol=sim)
mod2sp <- mod2se <- rep(NA, sim)
ESSmod2 <- compTime <-  rep(NA, sim)
for(s in 1:sim){
  ss <- 100000
  x <- cbind(rnorm(ss), rnorm(ss), rnorm(ss), rnorm(ss))
  beta <- c(0.7,-2, 0.5, -0.3)
  pT <- plogis(x%*%beta)
  yT <- rbinom(ss, 1, pT)
  yC <- yT
  cR <- 0.025
  cR2 <- 0.25
  ru <- runif(ss)
  yC[yT==0] <- 0 + (ru[yT==0] < cR)
  yC[yT==1] <- 1 - (ru[yT==1] < cR2)
  
  cv1 <- rnorm(ss)
  x2T <- cbind(1, cv1,  as.numeric(yT==1))
  x2C <- cbind(1, cv1,  as.numeric(yC==1))
  beta2 <- c(0.2, -0.7, 0.6)
  pT2 <- plogis(x2T%*%beta2)
  resp <- rbinom(ss, 1, pT2)
  

  pik <- inclusionprobabilities(rep(1, length(resp)), expSS)
  ind <- UPpoisson(pik)
  
  xSampTT <- x2T[as.logical(ind),]
  xSampTC <- x2C[as.logical(ind),]
  xSampB <- x[as.logical(ind),]
  yCsamp <- yC[as.logical(ind)]
  respSamp <- resp[as.logical(ind)]
  wgtSamp <- 1/pik[as.logical(ind)]
  wgtSamp <- wgtSamp/sum(wgtSamp)*sum(ind)
  
  mod <- sampler(xSampTC, respSamp, iter=2500, burn=500, wgt=wgtSamp)
  mod1B[,s] <- colMeans(mod$Beta)
  
  t1 <- Sys.time()
  mod2 <- regME2(x1=xSampTC[,-c(3)], x2=xSampB, yME=yCsamp, y=respSamp, iter=2500, burn=500, wgt=wgtSamp)
  compTime[s] <- as.numeric(Sys.time() - t1)
  mod2B[,s] <- colMeans(mod2$Beta1)
  mod2sp[s] <- mean(mod2$sp)
  mod2se[s] <- mean(mod2$se)
  ESSmod2[s] <- ESS(mod2$Beta1[,3])
  
  print(s)
}

rowMeans((mod1B - beta2)^2)
rowMeans((mod2B - beta2)^2)

abs(rowMeans(mod1B) - beta2)
abs(rowMeans(mod2B) - beta2)



coefs <- factor(paste0("b",0:2), labels=c(expression(beta[0]), expression(beta[1]), expression(beta[2])))
mod1DF <- data.frame(t(mod1B))
names(mod1DF) <- coefs
mod1DF$Mod <- "Naive"

mod2DF <- data.frame(t(mod2B))
names(mod2DF) <- coefs
mod2DF$Mod <- "Mixture"


modDF <- rbind(mod1DF, mod2DF) %>% pivot_longer(1:3)


dfNew <- data.frame(name=coefs, value=c(0.2, -0.7, 0.6))

ggplot(modDF %>% filter(name == "beta[2]"), aes(x=Mod, y=value))+
  geom_boxplot(fill='green')+
  geom_hline(data=dfNew %>% filter(name == "beta[2]"), color="red", aes(yintercept=value))+
  facet_wrap(~name, scales='free', labeller = label_parsed)+
  xlab("")+
  ylab("Estimate")+
  theme_classic()+
  ggtitle("Simple Random Sample")


