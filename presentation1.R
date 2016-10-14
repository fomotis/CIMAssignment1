## ----prep, message=FALSE,warning=FALSE,echo=FALSE------------------------
library(tibble)
library(dplyr)
library(purrr)
library(magrittr)
library(knitr)
library(fBasics)
library(fitdistrplus)
library(extRemes)
library(epitools)
knitr::opts_chunk$set(echo=F,tidy=T,external = F,message = FALSE,warning=F,fig.height = 3,fig.width = 3)

## ----explore3------------------------------------------------------------
data("cars")
distance <- cars$dist
set.seed(5)
boot_length <- c(5000, 50000, 100000)
n <- length(distance)
#fitting lognormal
#fit_lnormal <- fitdist(distance,"lnorm")
#aic_fit_lnormal <- (-2*fit_lnormal$loglik) + (2*length(fit_lnormal$estimate))
#fitting gamma
#fit_gamma <- fitdist(distance,"gamma")
#aic_fit_gamma <- (-2*fit_gamma$loglik) + (2*length(fit_gamma$estimate))
#fitting weibull
fit_wei <- fitdist(distance,"weibull")
aic_fit_wei <- (-2*fit_wei$loglik) + (2*length(fit_wei$estimate))
#random draws from a weibull with estimated parameter values
rswei <- sort(rweibull(n,fit_wei$estimate[[1]],fit_wei$estimate[[2]]))
denrswei <- dweibull(rswei,fit_wei$estimate[[1]],fit_wei$estimate[[2]])
#fitting gumbel
fit_gum <- fevd(distance,type="Gumbel",method='MLE')
#fitting generalized extremevalue
fit_gev <- fevd(distance,type="GEV",method='MLE')
hist(distance,col=2,probability = T)
lines(rswei,denrswei,col=1,lwd=2)
#legend("topright","weibull density",lwd=2,col=1)

## ----nonpara3------------------------------------------------------------
set.seed(5)
results3 <- list()
for(j in 1:length(boot_length)){
    rest <- unlist(map(1:boot_length[j],function(i){
    draws <- sample(distance,n,replace = T)
    max(draws)
  }))
  results3[[j]] <- rest
}
serror3 <- round(unlist(map(results3,sd)),2)
#or
#out <- numeric(B1)
#for(i in 1:B1){
#  draws <- sample(distance,n,replace = T)
#  out[[i]] <- max(draws)
#}

## ----parametric3---------------------------------------------------------
set.seed(500)
results32 <- list()
for(j in 1:length(boot_length)){
    rest <- unlist(map(1:boot_length[j],function(i){
    draws <- rweibull(n,fit_wei$estimate[[1]],fit_wei$estimate[[2]])
    max(draws)
  }))
  results32[[j]] <- rest
}
serror32 <- round(unlist(map(results32,sd)),2)

## ----probability---------------------------------------------------------
set.seed(5)
results33 <- list()
for(j in 1:length(boot_length)){
    rest <- unlist(map(1:boot_length[j],function(i){
    draws <- sample(distance,n,replace = T)
    sum(draws>65)/n
  }))
  results33[[j]] <- rest
}
prob33 <- round(unlist(map(results33,mean)),2)

## ----odds----------------------------------------------------------------
resp<-as.factor(c(rep(0,189),rep(1,10845),rep(0,104),rep(1 ,10933)))
trt<-as.factor(c(rep(1,189),rep(1,10845),rep(2,104),rep(2,10933)))
Aspirin.L <- tibble(trt,resp)
Aspirin.1<-table(trt,resp)
#Aspirin.1
row.names(Aspirin.1)=c("Placebo","Aspirin")
colnames(Aspirin.1) <- c("Yes","No")
oddsr <- round(oddsratio(Aspirin.1,method = 'wald')$measure[2,],2)
seoddsr <- sqrt(sum(1/Aspirin.1))

## ----oddsplot------------------------------------------------------------
plot(Aspirin.1,col=heat.colors(2),xlab='Treatment',ylab='Heart Disease')

## ----oddsbootstrapt,cache=TRUE-------------------------------------------
set.seed(5)
n4 <- length(Aspirin.L$trt)
results41 <- list()
#bootstrap algorithm
for(j in 1:(length(boot_length)-1)){
    rest <- map(1:boot_length[j],function(i){
    draws <- sample(1:n4,n4,replace = T)
    data_draws <- Aspirin.L %>% slice(draws)
    tab_data_draws <- table(data_draws)
    or_draws <- (tab_data_draws[1,1]*tab_data_draws[2,2])/(tab_data_draws[1,2]*tab_data_draws[2,1])
    log_or_draws <- log(or_draws)
    se_draws <- sqrt(sum(1/tab_data_draws))
    z_draws <- (mean(log_or_draws) - log(oddsr[[1]]))/se_draws
    return(list(or_star=or_draws,logor_star=log_or_draws,z_star=z_draws))
  })
  results41[[j]] <- rest
}
#separating the 5000 from 50000 replicates
list5000 <- results41[[1]]
list50000 <- results41[[2]]
#orstar and zstar for 5000 replicates
orstar5000 <- unlist(map(1:length(list5000),function(i){
  list5000[[i]]$or_star
}))
#
zstar5000 <- unlist(map(1:length(list5000),function(i){
  list5000[[i]]$z_star
}))
#orstar and zstar for 50000 replicates
orstar50000 <- unlist(map(1:length(list50000),function(i){
  list50000[[i]]$or_star
}))
#
zstar50000 <- unlist(map(1:length(list50000),function(i){
  list50000[[i]]$z_star
}))
#confidence interval and quantiles for 5000 replicates
q5000 <- round(quantile(zstar5000,c(0.025,0.975)),2)
ci5000 <- round(exp(log(oddsr[1]) + (q5000*seoddsr)),2)
c5000 <-  1 - ((sum(zstar5000 < quantile(zstar5000,0.025)) + sum(zstar5000 > quantile(zstar5000,0.975)))/length(zstar5000))
#confidence interval and quantiles for 50000 replicates
q50000 <- round(quantile(zstar50000,c(0.025,0.975)),2)
ci50000 <- round(exp(log(oddsr[1]) + (q50000*seoddsr)),2)
c50000 <-  1 - ((sum(zstar50000 < quantile(zstar50000,0.025)) + sum(zstar50000 > quantile(zstar50000,0.975)))/length(zstar50000))

## ----oddsbootstrappi,cache=T---------------------------------------------
set.seed(5)
results42 <- list()
#bootstrap percentile confidence interval algorithm
for(j in 1:(length(boot_length)-1)){
    rest <- unlist(map(1:boot_length[j],function(i){
    draws <- sample(1:n4,n4,replace = T)
    data_draws <- Aspirin.L %>% slice(draws)
    tab_data_draws <- table(data_draws)
    or_draws <- (tab_data_draws[1,1]*tab_data_draws[2,2])/(tab_data_draws[1,2]*tab_data_draws[2,1])
    #log_or_draws <- log(or_draws)
    #se_draws <- sqrt(sum(1/tab_data_draws))
    #z_draws <- (mean(log_or_draws) - log(oddsr[[1]]))/se_draws
    #return(list(or_star=or_draws,logor_star=log_or_draws,z_star=z_draws))
  }))
  results42[[j]] <- rest
}
pi5000 <- round(quantile(results42[[1]],c(0.025,0.975)),2)
pi50000 <- round(quantile(results42[[2]],c(0.025,0.975)),2)
#coverage
coverage5000 <- 1 - ((sum(results42[[1]] < quantile(results42[[1]],0.025)) + sum(results42[[1]] > quantile(results42[[1]],0.975)))/length(results42[[1]]))
coverage50000 <- 1 - ((sum(results42[[2]] < quantile(results42[[2]],0.025)) + sum(results42[[2]] > quantile(results42[[2]],0.975)))/length(results42[[2]]))

## ----plots3, fig.height=5,fig.width=7------------------------------------
par(mfrow=c(2,2))
for(i in seq_along(results3)){
  hist(results3[[i]],main=paste(boot_length[i],'replicates',sep=' '),xlab=expression(theta^'*'))
}

## ----plots32, fig.height=5,fig.width=7-----------------------------------
par(mfrow=c(2,2))
for(i in seq_along(results32)){
  hist(results32[[i]],main=paste(boot_length[i],'replicates',sep=' '),xlab=expression(theta^'*'))
}

## ----plots33, fig.height=5,fig.width=7-----------------------------------
par(mfrow=c(2,2))
for(i in seq_along(results33)){
  hist(results33[[i]],main=paste(boot_length[i],'replicates',sep=' '),xlab=expression(theta^'*'))
}

## ----plots41, fig.height=5,fig.width=7-----------------------------------
par(mfrow=c(2,2))
hist(orstar5000,xlab=expression(theta^'*'),main='5000 replicates')
#abline(v=quantile(orstar5000),c(0.025,0.975))
hist(orstar50000,xlab=expression(theta^'*'),main='50000 replicates')
#abline(v=quantile(orstar50000),c(0.025,0.975))

## ------------------------------------------------------------------------
purl('presentation1.Rnw')

