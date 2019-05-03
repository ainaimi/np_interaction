# load libraries/functions
#args=(commandArgs(TRUE))
#setwd(".")

packages <- c("foreach","doParallel","doRNG","boot","rmutil","mvtnorm","gam","sandwich","ggplot2",
              "devtools","glmnet","data.table","rpart","ranger","nnet","arm","earth","e1071","tidyverse")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

devtools::install_github("ecpolley/SuperLearner")
library(SuperLearner)

install.packages("tmle",repos='http://lib.stat.cmu.edu/R/CRAN')
library(tmle)

time<-proc.time()
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }
nsim <- 20

cols <- c("glm","glm_strat","AIPW_strat","TMLE_strat","AIPW_auto","TMLE_auto")
cols0 <- paste0(cols,"0")
cols1 <- paste0(cols,"1")
cols <- c(cols0,cols1)
res.est <- data.frame(matrix(nrow=nsim,ncol=length(cols)));colnames(res.est) <- cols; res.se <- res.est

ranger_learner <- create.Learner("SL.ranger", tune=list(min.node.size = c(15,30)))
# xgboost_learner <- create.Learner("SL.xgboost", tune=list(minobspernode = c(15,30)))
# svm_learner <- create.Learner("SL.svm",tune=list(nu = c(.25,.5)))

sl.lib_mu <- c(ranger_learner$names)#,svm_learner$names,xgboost_learner$names)
sl.lib_pi <- c(ranger_learner$names)#,svm_learner$names,xgboost_learner$names)

##  true value
true1<-6
true2<--3
interaction_sim<-function(counter){
  set.seed(counter)
  # data management
  i<-counter
  n<-1500
  p<-2
  cat("Now running iteration",i,"with a sample size of",n,'\n');flush.console()

  # confounders
  sigma<-matrix(0,nrow=p,ncol=p);diag(sigma)<-1
  c <- rmvnorm(n, mean=rep(0,p), sigma=sigma)

  # design matrix for outcome model

  muMatT<-model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]",collapse="+"),")")))

  parms3<-c(3,3)
  parms4<-c(log(1.5),log(1.5))

  beta<-parms3;beta<-c(120,beta)
  # design matrix for propensity score model
  piMatT<-model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]",collapse="+"),")")))
  piMatT2<-model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]",collapse="+"),")")))
  theta<-c(-.5,parms4)
  theta2<-c(-.5,parms4)
  mu <- muMatT%*%beta
  # propensity score model
  pi <- expit(piMatT%*%theta);
  pi_m <- expit(piMatT2%*%theta2);
  x<-rbinom(n,1,pi)
  m<-rbinom(n,1,pi_m)
  # outcome model: true expsoure effect = 6
  y <- x*6 + m*6 - 3*x*m + mu + rnorm(n,0,6)

  # correct specification
  dat <- data.frame(c,m,x,y); colnames(dat)[1:ncol(c)] <- paste("c",1:ncol(c),sep="")

  C <- muMatT[,-1];colnames(C) <- paste("C",1:ncol(muMatT[,-1]), sep="");
  X <- x
  M <- m
  Y <- y

  C0 <- subset(C,M==0)
  C1 <- subset(C,M==1)

  X0 <- subset(X,M==0)
  X1 <- subset(X,M==1)

  Y0 <- subset(Y,M==0)
  Y1 <- subset(Y,M==1)

  mod1 <- lm(y~x+m+x*m+c1+c2,data=dat)


  mod2_0 <- lm(y~x+c1+c2,data=subset(dat,m==0))
  mod2_1 <- lm(y~x+c1+c2,data=subset(dat,m==1))

  folds<-c(10)
  cat("Number of cross-validation folds is",folds,'\n');flush.console()
  tmle_strat0 <- tmle(Y0,X0,C0,family="gaussian",
                      Q.SL.library=sl.lib_mu,g.SL.library=sl.lib_pi,
                      V=folds)

  tmle_strat1 <- tmle(Y1,X1,C1,family="gaussian",
                      Q.SL.library=sl.lib_mu,g.SL.library=sl.lib_pi,
                      V=folds)

  # ps.sl.coefNPT <- cbind(n,t(tmleNPT$g$coef))
  # mu.sl.coefNPT <- cbind(n,t(tmleNPT$Qinit$coef))

  pihat_0 <- tmle_strat0$g$g1W
  pihat_1 <- tmle_strat1$g$g1W

  pihat_0 <- ifelse(pihat_0 < 0.025,0.025,pihat_0)
  pihat_0 <- ifelse(pihat_0 > 1-0.025,1-0.025,pihat_0)

  pihat_1 <- ifelse(pihat_1 < 0.025,0.025,pihat_1)
  pihat_1 <- ifelse(pihat_1 > 1-0.025,1-0.025,pihat_1)

  muhat_0  <- tmle_strat0$Qinit$Q[,2]*X0+tmle_strat0$Qinit$Q[,1]*(1-X0)
  muhat_1  <- tmle_strat1$Qinit$Q[,2]*X1+tmle_strat1$Qinit$Q[,1]*(1-X1)

  muhat1_0  <- tmle_strat0$Qinit$Q[,2]
  muhat1_1  <- tmle_strat1$Qinit$Q[,2]

  muhat0_0  <- tmle_strat0$Qinit$Q[,1]
  muhat0_1  <- tmle_strat1$Qinit$Q[,1]


  # compute estimators
  res.est$glm0[i] <- coef(mod1)["x"]
  res.est$glm1[i] <- coef(mod1)["x"] + coef(mod1)["x:m"]

  res.est$glm_strat0[i] <- coef(mod2_0)["x"]
  res.est$glm_strat1[i] <- coef(mod2_1)["x"]

  lower_bound <- min(y) - max(y)
  upper_bound <- max(y) - min(y)

  aipw_strat0 <- mean((((2*X0-1)*(Y0 - muhat_0))/((2*X0-1)*pihat_0 + (1-X0)) + muhat1_0 - muhat0_0)) ## inside the mean is EIE_aipw
  aipw_strat1 <- mean((((2*X1-1)*(Y1 - muhat_1))/((2*X1-1)*pihat_1 + (1-X1)) + muhat1_1 - muhat0_1))
  aipw_strat0 <- ifelse(aipw_strat0>upper_bound,upper_bound,aipw_strat0)
  aipw_strat0 <- ifelse(aipw_strat0<lower_bound,lower_bound,aipw_strat0)

  aipw_strat1 <- ifelse(aipw_strat1>upper_bound,upper_bound,aipw_strat1)
  aipw_strat1 <- ifelse(aipw_strat1<lower_bound,lower_bound,aipw_strat1)

  res.est$AIPW_strat0[i]  <- aipw_strat0
  res.est$AIPW_strat1[i]  <- aipw_strat1

  res.est$TMLE_strat0[i]<-tmle_strat0$estimates$ATE$psi
  res.est$TMLE_strat1[i]<-tmle_strat1$estimates$ATE$psi

  # compute closed-form SEs
  res.se$glm0[i] <- sqrt(vcov(mod1)["x","x"])
  res.se$glm1[i] <- sqrt(vcov(mod1)["x","x"] + vcov(mod1)["m","m"] - 2*vcov(mod1)["x","m"])

  res.se$glm_strat0[i] <- sqrt(vcov(mod2_0)["x","x"])
  res.se$glm_strat1[i] <- sqrt(vcov(mod2_1)["x","x"])

  aipw_strat0_se <- sd((((2*X0-1)*(Y0 - muhat_0))/((2*X0-1)*pihat_0 + (1-X0)) + muhat1_0 - muhat0_0))/sqrt(nrow(C0))
  aipw_strat1_se <- sd((((2*X1-1)*(Y1 - muhat_1))/((2*X1-1)*pihat_1 + (1-X1)) + muhat1_1 - muhat0_1))/sqrt(nrow(C1))

  res.se$AIPW_strat0[i]  <- aipw_strat0_se
  res.se$AIPW_strat1[i]  <- aipw_strat1_se

  res.se$TMLE_strat0[i] <- sqrt(tmle_strat0$estimates$ATE$var.psi)
  res.se$TMLE_strat1[i] <- sqrt(tmle_strat1$estimates$ATE$var.psi)


  # print updating results

  res.est0 <- res.est %>% select(ends_with("0"))
  res.est1 <- res.est %>% select(ends_with("1"))

  res.se0 <- res.se %>% select(ends_with("0"))
  res.se1 <- res.se %>% select(ends_with("1"))

  res.cov0 <- res.est0-1.96*res.se0 < true1 & true1 < res.est0+1.96*res.se0
  res.cov1 <- res.est1-1.96*res.se1 < true2 & true2 < res.est1+1.96*res.se1

  res.width <- (res.est+1.96*res.se) - (res.est-1.96*res.se)
  tmp <- data.frame(rbind(c(n,apply(res.est0-true1,2,mean,na.rm=T),apply(res.est1-true2,2,mean,na.rm=T)),
                          c(n,apply((res.est0-true1)^2,2,mean,na.rm=T),apply((res.est1-true2)^2,2,mean,na.rm=T)),
                          c(n,apply(res.cov0,2,mean,na.rm=T),apply(res.cov1,2,mean,na.rm=T)),
                          c(n,apply(res.width,2,mean,na.rm=T))))
  tmp.se <- data.frame(rbind(c(n,apply(res.se0,2,mean,na.rm=T),apply(res.se1,2,mean,na.rm=T))))
  rownames(tmp)<-c("bias","rmse","cov","width");colnames(tmp)[1]<-"N";print(round(tmp,3));cat('\n');flush.console()
  colnames(tmp.se)[1]<-"N"
  setDT(tmp, keep.rownames = TRUE)[];colnames(tmp)[1] <- "type"

  if(i==1){
    write.table(tmp,"./output/results.txt",sep="\t",row.names=F)
    write.table(tmp.se,"./output/results_se.txt",sep="\t",row.names=F)
  } else{
    write.table(tmp,"./output/results.txt",sep="\t",row.names=F,col.names=F,append=T)
    write.table(tmp.se,"./output/results_se.txt",sep="\t",row.names=F,col.names=F,append=T)
  }
  return(tmp)
}

interaction_sim(3)

start_time <- Sys.time()
cores<-detectCores()-2
print(cores)
results<-mclapply(1:5, function(x) interaction_sim(x),mc.cores=cores,mc.set.seed=F)
## run time:
Sys.time() - start_time
