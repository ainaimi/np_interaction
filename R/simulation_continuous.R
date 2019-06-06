###################################### PRELIMINARY RESULTS ######################################
# Now running iteration 10 with a sample size of 1500                                           
# Number of cross-validation folds is 10                                                        
# N   glm0 glm_strat0 AIPW_strat0 TMLE_strat0 AIPW_auto0 TMLE_auto0   glm1 glm_strat1
# bias  1500 -0.519     -0.546      -0.117      -0.233     -4.954     -0.301  5.472      5.531
# rmse  1500  0.269      0.299       0.014       0.054     24.546      0.091 29.942     30.587
# cov   1500  1.000      1.000       1.000       1.000        NaN        NaN  0.000      0.000
# width 1500  1.693      1.705       1.280       1.277        NaN        NaN  1.874      2.115
# AIPW_strat1 TMLE_strat1 AIPW_auto1 TMLE_auto1
# bias        5.401       5.280     12.325      6.057
# rmse       29.173      27.878    151.893     36.682
# cov         0.000       0.000        NaN        NaN
# width       1.533       1.530        NaN        NaN
################################################################################################

# load libraries/functions
#args=(commandArgs(TRUE))
#setwd(".")

packages <- c("foreach","doParallel","boot","rmutil","mvtnorm","gam","sandwich","ggplot2",
              "devtools","glmnet","data.table","rpart","ranger","nnet","arm","earth","e1071","tidyverse")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

library(SuperLearner)
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
xgboost_learner <- create.Learner("SL.xgboost", tune=list(minobspernode = c(15,30)))
glmnet_learner <- create.Learner("SL.glmnet", tune=list(alpha = c(0,.5,1)))
svm_learner <- create.Learner("SL.svm",tune=list(nu = c(.25,.5)))

sl.lib_mu <- c(ranger_learner$names,svm_learner$names,glmnet_learner$names)
sl.lib_pi <- c(ranger_learner$names,svm_learner$names,glmnet_learner$names)

##  true value
true1<-6
true2<-3
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
  parms5<-c(3,3)
  
  beta<-parms3;beta<-c(120,beta)
  # design matrix for propensity score model
  piMatT<-model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]",collapse="+"),")")))
  piMatT2<-model.matrix(as.formula(paste("~(",paste("c[,",1:ncol(c),"]",collapse="+"),")")))
  theta<-c(-.5,parms4)
  theta2<-c(15,parms5)
  mu <- muMatT%*%beta
  # propensity score model
  pi <- expit(piMatT%*%theta);
  mu_m <- piMatT2%*%theta2
  x <- rbinom(n,1,pi)
  m <- (runif(n,0,mu_m)*13)/max(mu_m)
  # outcome model
  y <- x*6 + m*2.5 + x*(4*sqrt(9*m)*as.numeric(m<2) + as.numeric(m>=2)*(abs(m-6)^(2))) + mu + rnorm(n,0,6) 
  
  # mm <- seq(0,max(m),by=.1)
  # xx <- 1
  # y_true1 <- xx*6 + mm*2.5 + xx*(4*sqrt(9*mm)*as.numeric(mm<2) + as.numeric(mm>=2)*(abs(mm-6)^(2)))
  # xx <- 0
  # y_true0 <- xx*6 + mm*2.5 + xx*(4*sqrt(9*mm)*as.numeric(mm<2) + as.numeric(mm>=2)*(abs(mm-6)^(2)))
  # 
  # plot_dat <- tibble(psi=y_true1-y_true0,m=mm)
  # ggplot(plot_dat) + geom_point(aes(y=psi,x=m))
  
  # correct specification
  dat <- data.frame(c,m,x,y); colnames(dat)[1:ncol(c)] <- paste("c",1:ncol(c),sep="")
  
  C <- muMatT[,-1];colnames(C) <- paste("C",1:ncol(muMatT[,-1]), sep="");
  X <- x
  M <- m
  Z <- cbind(C,M)
  Y <- y
  
  folds<-c(10)
  cat("Number of cross-validation folds is",folds,'\n');flush.console()
  tmle <- tmle(Y,X,Z,family="gaussian",
               Q.SL.library=sl.lib_mu,g.SL.library=sl.lib_pi,
               V=folds)
  

  
  # ps.sl.coefNPT <- cbind(n,t(tmleNPT$g$coef))
  # mu.sl.coefNPT <- cbind(n,t(tmleNPT$Qinit$coef))
  
  pihat <- tmle_strat0$g$g1W
  
  ## ASHLEY STOPPED HERE
  
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
  
  ## Auto versions
  #### TMLE version
  Z <- cbind(M,C)
  tmle_res <- tmle(Y,X,Z,family="gaussian",
                   Q.SL.library=sl.lib_mu,g.SL.library=sl.lib_pi,
                   V=folds)
  
  pihat <- tmle_res$g$g1W
  pihat <- ifelse(pihat < 0.025,0.025,pihat)
  pihat <- ifelse(pihat > 1-0.025,1-0.025,pihat)
  
  muhat        <- tmle_res$Qinit$Q[,2]*X+tmle_res$Qinit$Q[,1]*(1-X)
  muhat_star   <- tmle_res$Qstar[,2]*X+tmle_res$Qstar[,1]*(1-X)
  muhat1  <- tmle_res$Qinit$Q[,2]
  muhat0  <- tmle_res$Qinit$Q[,1]
#<<<<<<< HEAD


#  aipw_EFF <- as.numeric((((2*X-1)*(Y - muhat))/((2*X-1)*pihat + (1-X)) + muhat1 - muhat0)) ## inside the mean is EIE_aipw
#=======
  
  # Take the AIPW EFF funciton and replace muhat
  tmle_EFF <- as.numeric((((2*X-1)*(Y - muhat_star))/((2*X-1)*pihat + (1-X)) + tmle_res$Qstar[,2] - tmle_res$Qstar[,1]))
  #aipw_EFF <- as.numeric((((2*X-1)*(Y - muhat))/((2*X-1)*pihat + (1-X)) + muhat1 - muhat0))
  
  aipw_EFF <- as.numeric((((2*X-1)*(Y - muhat))/((2*X-1)*pihat + (1-X)) + tmle_res$Qinit$Q[,2] - tmle_res$Qinit$Q[,1]))
  ## NOTE: this part -> {+ tmle_res$Qstar[,2] - tmle_res$Qstar[,1]} is the tmle_EFF
  ## What about muhat?
  
#>>>>>>> a635b7d546acd9e4911e3baeae093332c1a0d0b6
  
  ## we have to figure out what the best way to do this is:
  Y_tmle_eff=tmle_EFF
  Y_aipw_eff=aipw_EFF
  X=data.frame(M,C1=C[,1],C2=C[,2])
  
  tmle_eff <- SuperLearner(Y=Y_tmle_eff,
                           X=X,
                           family = gaussian(),
                           SL.library=sl.lib_mu,
                           method = "method.NNLS")
  aipw_eff <- SuperLearner(Y=Y_aipw_eff,
                           X=X,
                           family = gaussian(),
                           SL.library=sl.lib_mu,
                           method = "method.NNLS")
  
  
  
  mu_tmle_eff0 <- predict(tmle_eff,newdata=subset(X,M==0),onlySL=T)$pred
  mu_tmle_eff1 <- predict(tmle_eff,newdata=subset(X,M==1),onlySL=T)$pred
  
  mu_aipw_eff0 <- predict(aipw_eff,newdata=subset(X,M==0),onlySL=T)$pred
  mu_aipw_eff1 <- predict(aipw_eff,newdata=subset(X,M==1),onlySL=T)$pred
  
  
#<<<<<<< HEAD
  mean(mu_aipw_eff0)
  mean(mu_aipw_eff1)
#=======
  # compute estimators
  res.est$glm0[i] <- coef(mod1)["x"]
  res.est$glm1[i] <- coef(mod1)["x"] + coef(mod1)["x:m"]
  
  res.est$glm_strat0[i] <- coef(mod2_0)["x"]
  res.est$glm_strat1[i] <- coef(mod2_1)["x"]
  
  lower_bound <- min(y) - max(y)
  upper_bound <- max(y) - min(y)
  
  ## ASK: Do we take the sd(aipw_strat)? including the mean?
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
  
  res.est$AIPW_auto0[i] <- mean(mu_aipw_eff0)
  res.est$AIPW_auto1[i] <- mean(mu_aipw_eff1)
  
  res.est$TMLE_auto0[i] <- mean(mu_tmle_eff0)
  res.est$TMLE_auto1[i] <- mean(mu_tmle_eff1)
  
  # SEs corrections
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
  
  ### AUTO TMLE / AIPW SEs
  
  
  # print updating results
  
  res.est0 <- res.est %>% select(ends_with("0"))
  res.est1 <- res.est %>% select(ends_with("1"))
  
  res.se0 <- res.se %>% select(ends_with("0"))
  res.se1 <- res.se %>% select(ends_with("1"))
  
  res.cov0 <- res.est0-1.96*res.se0 < true1 & true1 < res.est0+1.96*res.se0
  res.cov1 <- res.est1-1.96*res.se1 < true2 & true2 < res.est1+1.96*res.se1
  
  ### AUTO TMLE / AIPW Cov
  
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

results<-lapply(1:5, function(x) interaction_sim(x))

start_time <- Sys.time()
cores<-detectCores()-2
print(cores)
results<-mclapply(1:5, function(x) interaction_sim(x),mc.cores=cores,mc.set.seed=F)
## run time:
Sys.time() - start_time

res <- do.call(rbind,results)

apply(res[type=="bias",-c("type","N")],2,mean)
