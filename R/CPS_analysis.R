library(Matrix)
library(BayesLogit)
library(mvtnorm)
library(sampling)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(readr)
library(kableExtra)
source("R/functions.R")


cps <- read_csv('Data/cps_2008_2016_v2.csv')
cps16 <- cps %>% filter(year==2016)
cps16$entryc[is.na(cps16$entryc)] <- 0
cps16 <- cps16 %>% mutate(inscov = as.factor(inscov),
                          natz = as.factor(natz),
                          entryc = as.factor(entryc),
                          birthreg = as.factor(birthreg),
                          male = as.factor(male),
                          marst = as.factor(marst),
                          ownchild = as.factor(ownchild),
                          educrec = as.factor(educrec),
                          ownhome = as.factor(ownhome),
                          race = as.factor(race),
                          veteran = as.factor(veteran),
                          govtemp = as.factor(govtemp),
                          disability = as.factor(disability),
                          work = as.factor(work),
                          statereg = as.factor(statereg)) %>%
  filter(age < 65 & age >= 18)

cps16$rowid <- 1:nrow(cps16)

cps16 <- cps16 %>% filter(perwt != 0 & !(is.na(marst))& !(is.na(povpct))) %>% mutate(scaledWgt = perwt/mean(perwt))

cps16_2 <- cps16 %>% filter(forborn == 1 )  %>% mutate(scaledWgt = perwt/mean(perwt))

xSamp2 <- scale(model.matrix(~ entryc + birthreg +  age + male + marst + ownchild + hhhsize +
                               educrec + povpct + ownhome , data=cps16_2)[,-c(1)])

xSamp1 <- model.matrix(  ~ race + age + male + marst + ownchild + hhhsize + educrec + veteran + govtemp , data=cps16_2) ## V2

yME <- as.numeric(cps16_2$natz) - 1
ySamp <- as.numeric(cps16_2$inscov) - 1
ids <- which(cps16_2$rowid %in% cps16_2$rowid)


set.seed(1)
##### MIXTURE MODEL #####
modPost <- regME2(x1=xSamp1, x2=xSamp2, yME=yME, y=ySamp, iter=10000,burn=1000, wgt2=cps16_2$scaledWgt, ids)


modDFPost <- data.frame(Coefficient=c(colnames(xSamp1),"natz"), 
                        Estimate=colMeans(modPost$Beta1),
                        SE=apply(modPost$Beta1, 2, sd),
                        Percentile_2.5=apply(modPost$Beta1, 2, quantile, probs=0.025),
                        Percentile_97.5=apply(modPost$Beta1, 2, quantile, probs=0.975))

##### NAIVE MODEL #####
xSamp <- model.matrix(  ~ race + age + male + marst + ownchild + hhhsize + educrec + veteran + govtemp + natz, data=cps16_2)
ySamp <- as.numeric(cps16_2$inscov) - 1
modPostN <- pgBin(x=xSamp, y=ySamp, iter=10000, burn=1000, wgt= cps16_2$scaledWgt)
modDFPostN <- data.frame(Coefficient=colnames(xSamp), 
                         Estimate=colMeans(modPostN$Beta),
                         SE=apply(modPostN$Beta, 2, sd),
                         Percentile_2.5=apply(modPostN$Beta, 2, quantile, probs=0.025),
                         Percentile_97.5=apply(modPostN$Beta, 2, quantile, probs=0.975))

##### TABLES #####
modDFPost[,-1] <- round(modDFPost[,-1], 3)
modDFPost %>%
  kbl(caption="Post-ACA model summary of estimated health insurance regression coefficients using 2016 CPS data.",
      format="latex",
      col.names = c("Coefficient", "Posterior Mean", "Standard Error", "2.5%", "97.5%"))

modDFPostN[,-1] <- round(modDFPostN[,-1], 3)
modDFPostN %>%
  kbl(caption="Post-ACA Naive model summary of estimated health insurance regression coefficients using 2016 CPS data.",
      format="latex",
      col.names = c("Coefficient", "Posterior Mean", "Standard Error", "2.5%", "97.5%"))

##### PLOTS #####

df1 <- data.frame(Naive=modPostN$Beta[,19], Corrected=modPost$Beta1[,19]) %>% 
  pivot_longer(1:2, names_to="Model", values_to="Draw")
ggplot(df1)+
  geom_density(alpha= 0.4, aes(x=Draw, fill=Model))+
  theme_classic(base_size = 18)+
  xlab(expression(beta[c]))+
  ylab("Posterior Density")



xPred <- c(1,0,0,0,mean(cps16_2$age),0,0,0,0,0,0,0,mean(cps16_2$hhhsize),0,0,0,0,0)

df2 <- data.frame(Citizen=c(plogis(c(xPred,1)%*%t(modPost$Beta1))), `Noncitizen`=c(plogis(c(xPred,0)%*%t(modPost$Beta1)))) %>%
  pivot_longer(1:2, names_to = "Type")
ggplot(df2)+
  geom_density(alpha= 0.4, aes(x=value, fill=Type))+
  theme_classic(base_size = 18)+
  xlab("Probability of Health Insurance")+
  ylab("Posterior Density") 



