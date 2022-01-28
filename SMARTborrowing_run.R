#code used to yield simulation results from: Within-Trial Data Borrowing for Sequential Multiple Assignment Randomized Trials

#libraries
library(basket)
library(rjags)
library(parallel)
library(foreach)
library(itertools)
library(GenSA)
library(tidyverse)
cores=10

#first load function file
source("path/SMARTborrowing_functions.R")


#function that draws from the mixture for beta
rmixture = function(wts, betas){
  ##random generation of the indices
  id = sample(1:length(wts),prob=wts,size=nrow(betas),replace=TRUE)  
  id = cbind(1:nrow(betas),id)
  betas[id]
}


# #simulation

doSim <- function(seed, n, p1, p2, b0, b1, b2, b3, v0, v1, v2, v3, method, pi_omega) {
  set.seed(seed)
  r1 <- rbinom(250, size=1, p=p1)
  r2 <- rbinom(250, size=1, p=p2)
  
  dat <- data.frame(id=1:500, response=c(r1,r2), t1=c(rep(1,250), rep(2, 250)))
  #first treatment assignment is 1 or 2
  #second trt assignment is same for responders and 3 or 4 for nonresponders
  dat$t2 <- NA
  dat$t2[dat$response==1 & dat$t1==1] <- 1
  dat$t2[dat$response==1 & dat$t1==2] <- 2
  
  ind2 <- sample(1:nrow(dat[dat$response==0 & dat$t1==1,]), size = floor(0.5 * nrow(dat[dat$response==0 & dat$t1==1,])))
  dat$t2[dat$response==0 & dat$t1==1][ind2] <- 3
  dat$t2[dat$response==0 & dat$t1==1][-ind2] <- 4
  
  ind2 <- sample(1:nrow(dat[dat$response==0 & dat$t1==2,]), size = floor(0.5 * nrow(dat[dat$response==0 & dat$t1==2,])))
  dat$t2[dat$response==0 & dat$t1==2][ind2] <- 3
  dat$t2[dat$response==0 & dat$t1==2][-ind2] <- 4
  
  dat$psuccess <- NA
  
  dat$psuccess[dat$response==1] <- b0 + b1*1*(dat$t1[dat$response==1]==1)
  dat$psuccess[dat$response==0] <- v0 + v1*1*(dat$t1[dat$response==0]==1) + v2*1*(dat$t2[dat$response==0]==3) + v3*(dat$t1[dat$response==0]==1) * (dat$t2[dat$response==0]==3)
  
  dat$psuccess[dat$psuccess>1] <- 0.99
  dat$psuccess[dat$psuccess<0] <- 0.01
  
  dat$y <- NA
  for(i in unique(dat$psuccess)) {
    dat$y[dat$psuccess==i] <- rbinom(n=sum(dat$psuccess==i), size=1, p=i)
  }
  
  dat$baskets <- NA
  dat$baskets[dat$t1==1 & dat$t2==1] <- "11"
  dat$baskets[dat$t1==1 & dat$t2==3] <- "13"
  dat$baskets[dat$t1==1 & dat$t2==4] <- "14"
  dat$baskets[dat$t1==2 & dat$t2==2] <- "22"
  dat$baskets[dat$t1==2 & dat$t2==3] <- "23"
  dat$baskets[dat$t1==2 & dat$t2==4] <- "24"
  
  dat_wide <- data.frame(baskets = c("11", "13", "14", "22", "23", "24"))
  dat_wide$size <- as.vector(table(dat$baskets))
  dat_wide$y <- as.vector(table(dat$y, dat$baskets)[2,])
  dat_wide$truepsuccess <- tapply(dat$psuccess, dat$baskets, mean)
  
  #wide1 is responders, wide 2 is nonresponders
  dat_wide1 <- dat_wide[dat_wide$baskets %in% c("11", "12", "21", "22"),]
  dat_wide2 <- dat_wide[dat_wide$baskets %in% c("13", "14", "23", "24"),]
  
  #response rates beta-bin model
  
  #first group
  
  model_string <- "model{
  
  # Likelihood
  Y ~ dbinom(theta,n)
  
  # Prior
  theta ~ dbeta(1, 1)
}"

  model <- jags.model(textConnection(model_string), 
                      data = list(Y=sum(dat$response[dat$t1==1]),n=length(dat$response[dat$t1==1])))
  
  
  update(model, 5000, progress.bar="none"); # Burnin for 5000 samples
  
  samp <- coda.samples(model, 
                       variable.names=c("theta"), 
                       n.iter=100000, progress.bar="none")
  
  r1 <- as.matrix(samp[[1]])
  
  
  
  #second group
  
  model <- jags.model(textConnection(model_string), 
                      data = list(Y=sum(dat$response[dat$t1==2]),n=length(dat$response[dat$t1==2])))
  
  
  update(model, 5000, progress.bar="none"); # Burnin for 5000 samples
  
  samp <- coda.samples(model, 
                       variable.names=c("theta"), 
                       n.iter=100000, progress.bar="none")
  
  r2 <- as.matrix(samp[[1]])
  
  
  if(method=="new"){
    #run basket package
    
    #responders
    #prior matrix: conservative
    mat <- matrix(pi_omega, nrow=nrow(dat_wide1), ncol=nrow(dat_wide1))
    diag(mat) <- 1
    
    res1 <- mem_mcmc2(dat_wide1$y, dat_wide1$size, dat_wide1$baskets, prior = mat, mcmc_burnin = 1000, mcmc_iter = 10000, shape1 = 1, shape2 = 1)

    
    a <- 1
    b <- 1
    y_11 <- dat_wide1$y[1]
    y_22 <- dat_wide1$y[2]
    n_11 <- dat_wide1$size[1]
    n_22 <- dat_wide1$size[2]
    
    #model 1: no borrowing
    betas <- NULL
    #posterior for 11 with no borrowing
    betas <- cbind(betas, rbeta(100000, y_11+a, b+n_11-y_11))
    #posterior for 22 with no borrowing
    betas <- cbind(betas, rbeta(100000, y_22+a, b+n_22-y_22))
    
    
    #model 2: borrowing
    #posterior for pooling
    betas <- cbind(betas, rbeta(100000, y_11+y_22+a, b+n_11+n_22-y_11-y_22))
    
    wts_r <- c(1-res1$basket$pep[1,2], res1$basket$pep[1,2])
    
    #draw with ties
    X11 <- rmixture(wts=wts_r, betas=betas[,c(1,3)])
    X22 <- rmixture(wts=wts_r, betas=betas[,c(2,3)])
    s1 <- data.frame(X11, X22)
    
    
    #nonresponders
    #prior matrix: conservative
    mat <- matrix(pi_omega, nrow=nrow(dat_wide2), ncol=nrow(dat_wide2))
    diag(mat) <- 1
    
    res2 <- mem_mcmc2(dat_wide2$y, dat_wide2$size, dat_wide2$baskets, prior = mat, mcmc_burnin = 1000, mcmc_iter = 10000, shape1 = 1, shape2 = 1)
    
    #their weights
    frs <- lapply(res2$map_list, function (x) lapply(res2$mem_samp, function(y) unlist(prod(x==y))))
    frs <- lapply(frs, function(x) mean(unlist(x)))
    
    models1 <- res2$basket$models
    models2 <- res2$basket$models[,c(2, 1, 3, 4)]
    models3 <- res2$basket$models[,c(2, 3, 1, 4)]
    models4 <- res2$basket$models[,c(2, 3, 4, 1)]
    
    wts_nr1 <- wts_nr2 <- wts_nr3 <- wts_nr4 <- rep(NA, 8)
    for (k in 1:nrow(res2$basket$models)) {
      f1<-lapply(res2$map_list, function(x) prod(x[1,]==models1[k,]))
      wts_nr1[k] <- sum(unlist(f1)*unlist(frs))
      
      f2<-lapply(res2$map_list, function(x) prod(x[2,]==models2[k,]))
      wts_nr2[k] <- sum(unlist(f2)*unlist(frs))
      
      f3<-lapply(res2$map_list, function(x) prod(x[3,]==models3[k,]))
      wts_nr3[k] <- sum(unlist(f3)*unlist(frs))
      
      f4<-lapply(res2$map_list, function(x) prod(x[4,]==models4[k,]))
      wts_nr4[k] <- sum(unlist(f4)*unlist(frs))
      
      
    }
    
    y_13 <- dat_wide2$y[1]
    y_14 <- dat_wide2$y[2]
    y_23 <- dat_wide2$y[3]
    y_24 <- dat_wide2$y[4]
    n_13 <- dat_wide2$size[1]
    n_14 <- dat_wide2$size[2]
    n_23 <- dat_wide2$size[3]
    n_24 <- dat_wide2$size[4]
    
    #each basket only needs 8 posteriors though, of these:
    ps1 <- rbeta(100000, y_13+a, b+n_13-y_13)
    ps2 <- rbeta(100000, y_14+a, b+n_14-y_14)
    ps3 <- rbeta(100000, y_23+a, b+n_23-y_23)
    ps4 <- rbeta(100000, y_24+a, b+n_24-y_24)
    ps12 <- rbeta(100000, y_13+y_14+a, b+n_13+n_14-y_13-y_14)
    ps13 <- rbeta(100000, y_13+y_23+a, b+n_13+n_23-y_13-y_23)
    ps14 <- rbeta(100000, y_13+y_24+a, b+n_13+n_24-y_13-y_24)
    ps23 <- rbeta(100000, y_14+y_23+a, b+n_14+n_23-y_14-y_23)
    ps24 <- rbeta(100000, y_14+y_24+a, b+n_14+n_24-y_14-y_24)
    ps34 <- rbeta(100000, y_23+y_24+a, b+n_23+n_24-y_23-y_24)
    ps123 <- rbeta(100000, y_13+y_14+y_23+a, b+n_13+n_14+n_23-y_13-y_14-y_23)
    ps124 <- rbeta(100000, y_13+y_14+y_24+a, b+n_13+n_14+n_24-y_13-y_14-y_24)
    ps134 <- rbeta(100000, y_13+y_23+y_24+a, b+n_13+n_23+n_24-y_13-y_23-y_24)
    ps234 <- rbeta(100000, y_14+y_23+y_24+a, b+n_14+n_23+n_24-y_14-y_23-y_24)
    ps1234 <- rbeta(100000, y_13+y_14+y_23+y_24+a, b+n_13+n_14+n_23+n_24-y_13-y_14-y_23-y_24)
    
    
    X13 <- rmixture(wts=wts_nr1, 
                    betas=data.frame(ps1, ps12, ps13, ps14, ps123, ps124, ps134, ps1234))
    
    X14 <- rmixture(wts=wts_nr2, 
                    betas=data.frame(ps2, ps12, ps23, ps24, ps123, ps124, ps234, ps1234))
    
    X23 <- rmixture(wts=wts_nr3, 
                    betas=data.frame(ps3, ps13, ps23, ps34, ps123, ps134, ps234, ps1234))
    
    X24 <- rmixture(wts=wts_nr4, 
                    betas=data.frame(ps4, ps14, ps24, ps34, ps124, ps134, ps234, ps1234))
    
    s2 <- data.frame(X13, X14, X23, X24)
    
    
    
    #calculate the posterior for overall expected outcome of all DTR
    
    finalsamples1 <- s1$X11*r1 + s2$X13*(1-r1)
    finalsamples2 <- s1$X11*r1 + s2$X14*(1-r1)
    finalsamples3 <- s1$X22*r2 + s2$X23*(1-r2)
    finalsamples4 <- s1$X22*r2 + s2$X24*(1-r2)
    
    pointest_r1 <- mean(finalsamples1)
    sdev_r1 <- sd(finalsamples1)
    CI_r1 <- c(quantile(finalsamples1, 0.025), quantile(finalsamples1, 0.975))
    
    pointest_r2 <- mean(finalsamples1-finalsamples2)
    sdev_r2 <- sd(finalsamples1-finalsamples2)
    CI_r2 <- c(quantile(finalsamples1-finalsamples2, 0.025), quantile(finalsamples1-finalsamples2, 0.975))
    
    
    maxvec <- pmax(finalsamples1, finalsamples2, finalsamples3, finalsamples4)
    probsbest_r3 <- c(sum(finalsamples1 == maxvec), sum(finalsamples2 == maxvec), sum(finalsamples3 == maxvec), sum(finalsamples4 == maxvec))/length(finalsamples1)
    probsbest2_r4 <- which.max(probsbest_r3)
    
    p12 <- mean(finalsamples1-finalsamples2>-0.01 & finalsamples1-finalsamples2<0.01)
    p13 <- mean(finalsamples1-finalsamples3>-0.01 & finalsamples1-finalsamples3<0.01)
    p14 <- mean(finalsamples1-finalsamples4>-0.01 & finalsamples1-finalsamples4<0.01)
    
    p23 <- mean(finalsamples2-finalsamples3>-0.01 & finalsamples2-finalsamples3<0.01)
    p24 <- mean(finalsamples2-finalsamples4>-0.01 & finalsamples2-finalsamples4<0.01)
    
    p34 <- mean(finalsamples3-finalsamples4>-0.01 & finalsamples3-finalsamples4<0.01)
    probssame_r5 <- c(p12, p13, p14, p23, p24, p34)
    
    # #clustering using MEM:
    # 
    #clustering using MEM weights (analytically):
    #first consider borrowing on response rates:
    mat <- matrix(pi_omega, nrow=2, ncol=2)
    diag(mat) <- 1
    dat_resp <- data.frame(baskets = c(1, 2), size = c(nrow(dat[dat$t1==1,]), nrow(dat[dat$t1==2,])), r = c(sum(dat$response[dat$t1==1]), sum(dat$response[dat$t1==2])))
    resr <- mem_mcmc2(dat_resp$r, dat_resp$size, dat_resp$baskets, prior = mat, mcmc_burnin = 1000, mcmc_iter = 10000, shape1 = 1, shape2 = 1)
    
    probrsame <- resr$basket$pep[1,2]
    
    p12 <- res2$basket$pep["13", "14"]
    p13 <- res1$basket$pep["11", "22"] * res2$basket$pep["13", "23"] * probrsame
    p14 <- res1$basket$pep["11", "22"] * res2$basket$pep["13", "24"] * probrsame
    
    p23 <- res1$basket$pep["11", "22"] * res2$basket$pep["14", "23"] * probrsame
    p24 <- res1$basket$pep["11", "22"] * res2$basket$pep["14", "24"] * probrsame
    
    p34 <- res2$basket$pep["23", "24"]
    
    probssame_r6 <- c(p12, p13, p14, p23, p24, p34)
    
    pointest_all <- c(mean(finalsamples1), mean(finalsamples2), mean(finalsamples3), mean(finalsamples4))
    
    
  } else if(method=="old"){
    #compare to the traditional method
    
    #beta bin model for each group separately
    model <- jags.model(textConnection(model_string), 
                        data = list(Y=sum(dat$y[dat$baskets=="11"]),n=length(dat$y[dat$baskets=="11"])))
    
    update(model, 1000, progress.bar="none"); # Burnin for 1000 samples
    samp <- coda.samples(model, variable.names=c("theta"), n.iter=100000, progress.bar="none")
    yr11 <- as.matrix(samp[[1]])
    
    model <- jags.model(textConnection(model_string),
                        data = list(Y=sum(dat$y[dat$baskets=="13"]),n=length(dat$y[dat$baskets=="13"])))
    
    update(model, 1000, progress.bar="none"); # Burnin for 1000 samples
    samp <- coda.samples(model, variable.names=c("theta"), n.iter=100000, progress.bar="none")
    yr13 <- as.matrix(samp[[1]])
    
    model <- jags.model(textConnection(model_string),
                        data = list(Y=sum(dat$y[dat$baskets=="14"]),n=length(dat$y[dat$baskets=="14"])))
    
    update(model, 1000, progress.bar="none"); # Burnin for 1000 samples
    samp <- coda.samples(model, variable.names=c("theta"), n.iter=100000, progress.bar="none")
    yr14 <- as.matrix(samp[[1]])
    
    model <- jags.model(textConnection(model_string),
                        data = list(Y=sum(dat$y[dat$baskets=="22"]),n=length(dat$y[dat$baskets=="22"])))
    
    update(model, 1000, progress.bar="none"); # Burnin for 1000 samples
    samp <- coda.samples(model, variable.names=c("theta"), n.iter=100000, progress.bar="none")
    yr22 <- as.matrix(samp[[1]])
    
    model <- jags.model(textConnection(model_string),
                        data = list(Y=sum(dat$y[dat$baskets=="23"]),n=length(dat$y[dat$baskets=="23"])))
    
    update(model, 1000, progress.bar="none"); # Burnin for 1000 samples
    samp <- coda.samples(model, variable.names=c("theta"), n.iter=100000, progress.bar="none")
    yr23 <- as.matrix(samp[[1]])
    
    model <- jags.model(textConnection(model_string),
                        data = list(Y=sum(dat$y[dat$baskets=="24"]),n=length(dat$y[dat$baskets=="24"])))
    
    update(model, 1000, progress.bar="none"); # Burnin for 1000 samples
    samp <- coda.samples(model, variable.names=c("theta"), n.iter=100000, progress.bar="none")
    yr24 <- as.matrix(samp[[1]])
    
    #calculate the posterior for overall expected outcome of all 4 DTRs
    finalsamples1_trad <- yr11*r1 + yr13*(1-r1)
    finalsamples2_trad <- yr11*r1 + yr14*(1-r1)
    finalsamples3_trad <- yr22*r1 + yr23*(1-r1)
    finalsamples4_trad <- yr22*r1 + yr24*(1-r1)
    
    pointest_r1 <- mean(finalsamples1_trad)
    sdev_r1 <- sd(finalsamples1_trad)
    CI_r1 <- c(quantile(finalsamples1_trad, 0.025), quantile(finalsamples1_trad, 0.975))
    
    pointest_r2 <- mean(finalsamples1_trad-finalsamples2_trad)
    sdev_r2 <- sd(finalsamples1_trad-finalsamples2_trad)
    CI_r2 <- c(quantile(finalsamples1_trad-finalsamples2_trad, 0.025), quantile(finalsamples1_trad-finalsamples2_trad, 0.975))
    
    
    maxvec <- pmax(finalsamples1_trad, finalsamples2_trad, finalsamples3_trad, finalsamples4_trad)
    probsbest_r3 <- c(sum(finalsamples1_trad == maxvec), sum(finalsamples2_trad == maxvec), sum(finalsamples3_trad == maxvec), sum(finalsamples4_trad == maxvec))/length(finalsamples1_trad)
    probsbest2_r4 <- which.max(probsbest_r3)
    
    #clustering
    p12 <- mean(finalsamples1_trad-finalsamples2_trad>-0.01 & finalsamples1_trad-finalsamples2_trad<0.01)
    p13 <- mean(finalsamples1_trad-finalsamples3_trad>-0.01 & finalsamples1_trad-finalsamples3_trad<0.01)
    p14 <- mean(finalsamples1_trad-finalsamples4_trad>-0.01 & finalsamples1_trad-finalsamples4_trad<0.01)
    
    p23 <- mean(finalsamples2_trad-finalsamples3_trad>-0.01 & finalsamples2_trad-finalsamples3_trad<0.01)
    p24 <- mean(finalsamples2_trad-finalsamples4_trad>-0.01 & finalsamples2_trad-finalsamples4_trad<0.01)
    
    p34 <- mean(finalsamples3_trad-finalsamples4_trad>-0.01 & finalsamples3_trad-finalsamples4_trad<0.01)
    
    probssame_r5 <- c(p12, p13, p14, p23, p24, p34)
    probssame_r6 <- rep(NA, 6)
    
    pointest_all <- c(mean(finalsamples1_trad), mean(finalsamples2_trad), mean(finalsamples3_trad), mean(finalsamples4_trad))
    
  }
  
  return(list(pointest_r1=pointest_r1, sdev_r1=sdev_r1, CI_r1=CI_r1, pointest_r2=pointest_r2, sdev_r2=sdev_r2, CI_r2=CI_r2, probsbest_r3=probsbest_r3, probsbest2_r4=probsbest2_r4, probssame_r5=probssame_r5, probssame_r6=probssame_r6, pointest_all=pointest_all))
  }



#run the simulation
nsim=1000
RNGkind("L'Ecuyer-CMRG")
pars=list(
  
  #sc 1 (2,4)
  list(n=500, p1=0.3, p2=0.45, b0 = 0.6, b1= 0, v0 = 0.4, v1= 0, v2= 0, v3= 0, pi_omega=0.05,  method="old"),
  list(n=500, p1=0.3, p2=0.45, b0 = 0.6, b1= 0, v0 = 0.4, v1= 0, v2= 0, v3= 0, pi_omega=0.05,  method="new"),
  list(n=500, p1=0.3, p2=0.45, b0 = 0.6, b1= 0, v0 = 0.4, v1= 0, v2= 0, v3= 0, pi_omega=0.1,  method="new"),
  
  #sc 2 (2,3)
  list(n=500, p1=0.3, p2=0.45, b0 = 0.6, b1= 0, v0 = 0.1, v1= 0.4, v2= 0.4, v3= -0.4, pi_omega=0.05,  method="old"),
  list(n=500, p1=0.3, p2=0.45, b0 = 0.6, b1= 0, v0 = 0.1, v1= 0.4, v2= 0.4, v3= -0.4, pi_omega=0.05,  method="new"),
  list(n=500, p1=0.3, p2=0.45, b0 = 0.6, b1= 0, v0 = 0.1, v1= 0.4, v2= 0.4, v3= -0.4, pi_omega=0.1,  method="new"),
  
  #sc 3 (0,3) NEW
  list(n=500, p1=0.3, p2=0.4, b0 = 0.6, b1= 0.25, v0 = 0.4, v1= -0.3, v2= 0, v3= 0.3, pi_omega=0.05,  method="old"),
  list(n=500, p1=0.3, p2=0.4, b0 = 0.6, b1= 0.25, v0 = 0.4, v1= -0.3, v2= 0, v3= 0.3, pi_omega=0.05,  method="new"),
  list(n=500, p1=0.3, p2=0.4, b0 = 0.6, b1= 0.25, v0 = 0.4, v1= -0.3, v2= 0, v3= 0.3, pi_omega=0.1,  method="new"),
  
  #sc 4 (0,2)
  list(n=500, p1=0.3, p2=0.45, b0 = 0.6, b1= 0.25, v0 = 0.05, v1= 0.15, v2= 0.45, v3= -0.15, pi_omega=0.05,  method="old"),
  list(n=500, p1=0.3, p2=0.45, b0 = 0.6, b1= 0.25, v0 = 0.05, v1= 0.15, v2= 0.45, v3= -0.15, pi_omega=0.05,  method="new"),
  list(n=500, p1=0.3, p2=0.45, b0 = 0.6, b1= 0.25, v0 = 0.05, v1= 0.15, v2= 0.45, v3= -0.15, pi_omega=0.1,  method="new"),
  
  #sc 5 (0,0)
  list(n=500, p1=0.3, p2=0.45, b0 = 0.6, b1= 0.25, v0 = 0.05, v1= 0.15, v2= 0.35, v3= 0.1, pi_omega=0.05,  method="old"),
  list(n=500, p1=0.3, p2=0.45, b0 = 0.6, b1= 0.25, v0 = 0.05, v1= 0.15, v2= 0.35, v3= 0.1, pi_omega=0.05,  method="new"),
  list(n=500, p1=0.3, p2=0.45, b0 = 0.6, b1= 0.25, v0 = 0.05, v1= 0.15, v2= 0.35, v3= 0.1, pi_omega=0.1,  method="new"),
  
  #sc 6 (0,0)
  list(n=500, p1=0.3, p2=0.45, b0 = 0.6, b1= 0.05, v0 = 0.1, v1= 0.2, v2= 0.25, v3= -0.05, pi_omega=0.05,  method="old"),
  list(n=500, p1=0.3, p2=0.45, b0 = 0.6, b1= 0.05, v0 = 0.1, v1= 0.2, v2= 0.25, v3= -0.05, pi_omega=0.05,  method="new"),
  list(n=500, p1=0.3, p2=0.45, b0 = 0.6, b1= 0.05, v0 = 0.1, v1= 0.2, v2= 0.25, v3= -0.05, pi_omega=0.1,  method="new")
  
  
)



summ1 <- matrix(nrow=length(pars),ncol=4)
summ2 <- matrix(nrow=length(pars),ncol=4)
summ3 <- matrix(nrow=length(pars),ncol=4)
summ4 <- matrix(nrow=length(pars),ncol=4)
summ5 <- matrix(nrow=length(pars),ncol=6)
summ6 <- matrix(nrow=length(pars),ncol=6)
summ7 <- matrix(nrow=length(pars),ncol=4)

is.between <- function(x, a, b) {
  (x - a)  *  (b - x) >= 0
}


pb= txtProgressBar(min = 0, max = length(pars), style = 3, char=":)")
bootdf <- list()
system.time(for(i in 1:length(pars)) {
  ## Define the parameters
  par <- pars[[i]] ####
  
  #### Do the simulation with these settings
  #set.seed(1)
  sim.results <- t(mclapply(1:nsim,doSim,n=par$n, p1=par$p1, p2=par$p2, b0=par$b0, b1=par$b1, v0=par$v0, v1=par$v1, v2=par$v2, v3=par$v3, pi_omega=par$pi_omega, method = par$method, mc.cores=cores, mc.silent=T))
  
  result <- unlist(sim.results)
  result <- matrix(result, ncol=29, byrow=T)
  bootdf[[i]] <- result
  
  result1 <- result[,1:4]
  result2 <- result[,5:8]
  result3 <- result[,9:12]
  result4 <- as.vector(result[,13])
  result5 <- result[,14:19]
  result6 <- result[,20:25]
  result7 <- result[,26:29]
  
  truth1 <- (par$b0+par$b1)*par$p1 + (par$v0+par$v1+par$v2 + par$v3)*(1-par$p1)
  
  cvrg1 <- is.between(truth1, result1[,3], result1[,4])
  ## Summarize the simulation results
  summ1[i,1] <- sum(cvrg1==T)/nsim
  summ1[i,2] <- mean(result1[,1] - truth1)
  summ1[i,3] <- sum((result1[,1] - truth1)**2)/nsim
  summ1[i,4] <- mean(result1[,2])
  
  truth2 <- ((par$b0+par$b1)*par$p1 + (par$v0+par$v1+par$v2 + par$v3)*(1-par$p1))  -  ((par$b0+par$b1)*par$p1 + (par$v0+par$v1)*(1-par$p1))
  cvrg2 <- is.between(truth2, result2[,3], result2[,4])
  ## Summarize the simulation results
  summ2[i,1] <- sum(cvrg2==T)/nsim
  summ2[i,2] <- mean(result2[,1] - truth2)
  summ2[i,3] <- sum((result2[,1] - truth2)**2)/nsim
  summ2[i,4] <- mean(result2[,2])
  
  ## Summarize the simulation results
  summ3[i,] <- colMeans(result3)
  
  
  ## Summarize the simulation results
  summ4[i,] <- c(sum(result4==1), sum(result4==2), sum(result4==3), sum(result4==4))/nsim
  
  ## Summarize the simulation results
  summ5[i,] <- colMeans(result5)
  
  summ6[i,] <- colMeans(result6)
  
  summ7[i,] <- colMeans(result7)
  
  #sink('output2not.txt')
  print("estimating one DTR outcome")
  print(summ1[1:i,])
  print("difference between 2 DTRs")
  print(summ2[1:i,])
  print("how many times selected as optimal")
  print(summ4[1:i,])
  print("used for clustering graph")
  print(summ6[1:i,])
  
  setTxtProgressBar(pb, i)
  #sink()
})
close(pb)




########################################################
######## PLOTS  ################################################
########################################################

library(ggplot2)
library(gridExtra)
library(tidyverse)

########################################################
######## 1 ################################################
########################################################

#relationship graph

#to create the matrix for clustering:
pc <- summ6[2,]
matc2 <- matrix(c(1, pc[1:3], NA, 1, pc[4:5], NA, NA, 1, pc[6], NA, NA, NA, 1), byrow=T, ncol=4)
matc2 <- round(matc2, 3)

color_by = c("mean_est")
layout = c("fr")
pep_cutoff = 0

pep <- matc2
pep[lower.tri(pep)] = t(pep)[lower.tri(pep)]
rownames(pep) <- c("DTR1", "DTR2","DTR3","DTR4")
colnames(pep) <- c("DTR1", "DTR2","DTR3","DTR4")

node_attrs <- tibble::tibble(name = c("DTR1", "DTR2","DTR3","DTR4"),
                             mean_est = summ7[2,])
graph <- tidygraph::as_tbl_graph(pep, directed = FALSE) %>%
  tidygraph::activate("nodes") %>% tidygraph::left_join(node_attrs,
                                                        by = "name")
legend_name <- c(mean_est = "Mean Trt\nSuccess")
graph <- graph %>% tidygraph::activate("edges")

p1r <- ggraph::ggraph(graph, layout = layout, weights = .data$weight) +
  ggraph::geom_edge_link(color = "gray", alpha = 0.6, aes(width = .data$weight)) +
  ggraph::geom_node_point(aes(color = .data[[color_by]]),
                          size = 7) + ggraph::geom_node_label(aes(label = .data$name),
                                                              nudge_y = 0.15, label.size = NA, hjust = "inward", alpha = 0.6) +
  ggplot2::scale_color_viridis_c(option = "plasma", name = legend_name[color_by]) +
  ggraph::scale_edge_width_continuous(name = "PEP", range=c(0.4,8), limits=c(0, 0.8), breaks=c(0.2, 0.4, 0.6)) #can add range=c(0,10) to play with width


########################################################
######## 2 ################################################
########################################################
#relationship graph

#to create the matrix for clustering:
pc <- summ6[5,]
matc2 <- matrix(c(1, pc[1:3], NA, 1, pc[4:5], NA, NA, 1, pc[6], NA, NA, NA, 1), byrow=T, ncol=4)
matc2 <- round(matc2, 4)

color_by = c("mean_est")
layout = c("fr")
pep_cutoff = 0

pep <- matc2
pep[lower.tri(pep)] = t(pep)[lower.tri(pep)]
rownames(pep) <- c("DTR1", "DTR2","DTR3","DTR4")
colnames(pep) <- c("DTR1", "DTR2","DTR3","DTR4")

node_attrs <- tibble::tibble(name = c("DTR1", "DTR2","DTR3","DTR4"),
                             mean_est = summ7[5,])
graph <- tidygraph::as_tbl_graph(pep, directed = FALSE) %>%
  tidygraph::activate("nodes") %>% tidygraph::left_join(node_attrs,
                                                        by = "name")
legend_name <- c(mean_est = "Mean Trt\nSuccess")
graph <- graph %>% tidygraph::activate("edges")

p2r <- ggraph::ggraph(graph, layout = layout, weights = .data$weight) +
  ggraph::geom_edge_link(color = "gray", alpha = 0.6, aes(width = .data$weight)) +
  ggraph::geom_node_point(aes(color = .data[[color_by]]),
                          size = 7) + ggraph::geom_node_label(aes(label = .data$name),
                                                              nudge_y = 0.15, label.size = NA, hjust = "inward", alpha = 0.6) +
  ggplot2::scale_color_viridis_c(option = "plasma", name = legend_name[color_by]) +
  ggraph::scale_edge_width_continuous(name = "PEP", range=c(0.5,8), limits=c(0, 0.7), breaks=c(0.2, 0.4, 0.8)) #can add range=c(0,10) to play with width

p2r <- ggraph::ggraph(graph, layout = layout, weights = .data$weight) +
  ggraph::geom_edge_link(color = "gray", alpha = 0.6, aes(width = .data$weight)) +
  ggraph::geom_node_point(aes(color = .data[[color_by]]),
                          size = 7) + ggraph::geom_node_label(aes(label = .data$name),
                                                              nudge_y = 0.15, label.size = NA, hjust = "inward", alpha = 0.6) +
  ggplot2::scale_color_viridis_c(option = "plasma", name = legend_name[color_by]) +
  ggraph::scale_edge_width_continuous(name = "PEP", range=c(0.5,8), breaks=c(0.2, 0.4, 0.6)) #can add range=c(0,10) to play with width


########################################################
######## 3 ################################################
########################################################
#relationship graph

#to create the matrix for clustering:
pc <- summ6[8,]
matc2 <- matrix(c(1, pc[1:3], NA, 1, pc[4:5], NA, NA, 1, pc[6], NA, NA, NA, 1), byrow=T, ncol=4)
matc2 <- round(matc2, 3)

color_by = c("mean_est")
layout = c("fr")
pep_cutoff = 0

pep <- matc2
pep[lower.tri(pep)] = t(pep)[lower.tri(pep)]
rownames(pep) <- c("DTR1", "DTR2","DTR3","DTR4")
colnames(pep) <- c("DTR1", "DTR2","DTR3","DTR4")

node_attrs <- tibble::tibble(name = c("DTR1", "DTR2","DTR3","DTR4"),
                             mean_est = summ7[8,])
graph <- tidygraph::as_tbl_graph(pep, directed = FALSE) %>%
  tidygraph::activate("nodes") %>% tidygraph::left_join(node_attrs,
                                                        by = "name")
legend_name <- c(mean_est = "Mean Trt\nSuccess")
graph <- graph %>% tidygraph::activate("edges")

p3r <- ggraph::ggraph(graph, layout = layout, weights = .data$weight) +
  ggraph::geom_edge_link(color = "gray", alpha = 0.6, aes(width = .data$weight)) +
  ggraph::geom_node_point(aes(color = .data[[color_by]]),
                          size = 7) + ggraph::geom_node_label(aes(label = .data$name),
                                                              nudge_y = 0.15, label.size = NA, hjust = "inward", alpha = 0.6) +
  ggplot2::scale_color_viridis_c(option = "plasma", name = legend_name[color_by]) +
  ggraph::scale_edge_width_continuous(name = "PEP", range=c(0.5,8), limit=c(0, 0.8), breaks=c(0.2, 0.4, 0.6)) #can add range=c(0,10) to play with width


########################################################
######## 4 ################################################
########################################################
#relationship graph

#to create the matrix for clustering:
pc <- summ6[11,]
matc2 <- matrix(c(1, pc[1:3], NA, 1, pc[4:5], NA, NA, 1, pc[6], NA, NA, NA, 1), byrow=T, ncol=4)
matc2 <- round(matc2, 2)

color_by = c("mean_est")
layout = c("fr")
pep_cutoff = 0

pep <- matc2
pep[lower.tri(pep)] = t(pep)[lower.tri(pep)]
rownames(pep) <- c("DTR1", "DTR2","DTR3","DTR4")
colnames(pep) <- c("DTR1", "DTR2","DTR3","DTR4")

node_attrs <- tibble::tibble(name = c("DTR1", "DTR2","DTR3","DTR4"),
                             mean_est = summ7[11,])
graph <- tidygraph::as_tbl_graph(pep, directed = FALSE) %>%
  tidygraph::activate("nodes") %>% tidygraph::left_join(node_attrs,
                                                        by = "name")
legend_name <- c(mean_est = "Mean Trt\nSuccess")
graph <- graph %>% tidygraph::activate("edges")

p4r <- ggraph::ggraph(graph, layout = layout, weights = .data$weight) +
  ggraph::geom_edge_link(color = "gray", alpha = 0.6, aes(width = .data$weight)) +
  ggraph::geom_node_point(aes(color = .data[[color_by]]),
                          size = 7) + ggraph::geom_node_label(aes(label = .data$name),
                                                              nudge_y = 0.15, label.size = NA, hjust = "inward", alpha = 0.6) +
  ggplot2::scale_color_viridis_c(option = "plasma", name = legend_name[color_by]) +
  ggraph::scale_edge_width_continuous(name = "PEP", range=c(0.5,8), breaks=seq(0, 1, 0.05)) #can add range=c(0,10) to play with width


########################################################
######## 5 ################################################
########################################################
#relationship graph

#to create the matrix for clustering:
pc <- summ6[14,]
matc2 <- matrix(c(1, pc[1:3], NA, 1, pc[4:5], NA, NA, 1, pc[6], NA, NA, NA, 1), byrow=T, ncol=4)
matc2 <- round(matc2, 4)

color_by = c("mean_est")
layout = c("fr")
pep_cutoff = 0

pep <- matc2
pep[lower.tri(pep)] = t(pep)[lower.tri(pep)]
rownames(pep) <- c("DTR1", "DTR2","DTR3","DTR4")
colnames(pep) <- c("DTR1", "DTR2","DTR3","DTR4")

node_attrs <- tibble::tibble(name = c("DTR1", "DTR2","DTR3","DTR4"),
                             mean_est = summ7[14,])
graph <- tidygraph::as_tbl_graph(pep, directed = FALSE) %>%
  tidygraph::activate("nodes") %>% tidygraph::left_join(node_attrs,
                                                        by = "name")
legend_name <- c(mean_est = "Mean Trt\nSuccess")
graph <- graph %>% tidygraph::activate("edges")

p5r <- ggraph::ggraph(graph, layout = layout, weights = .data$weight) +
  ggraph::geom_edge_link(color = "gray", alpha = 0.6, aes(width = .data$weight)) +
  ggraph::geom_node_point(aes(color = .data[[color_by]]),
                          size = 7) + ggraph::geom_node_label(aes(label = .data$name),
                                                              nudge_y = 0.15, label.size = NA, hjust = "inward", alpha = 0.6) +
  ggplot2::scale_color_viridis_c(option = "plasma", name = legend_name[color_by]) +
  ggraph::scale_edge_width_continuous(name = "PEP", range=c(0.5,8)) #can add range=c(0,10) to play with width


########################################################
######## 6 ################################################
########################################################

#relationship graph

#to create the matrix for clustering:
pc <- summ6[17,]
matc2 <- matrix(c(1, pc[1:3], NA, 1, pc[4:5], NA, NA, 1, pc[6], NA, NA, NA, 1), byrow=T, ncol=4)
matc2 <- round(matc2, 3)

color_by = c("mean_est")
layout = c("fr")
pep_cutoff = 0

pep <- matc2
pep[lower.tri(pep)] = t(pep)[lower.tri(pep)]
rownames(pep) <- c("DTR1", "DTR2","DTR3","DTR4")
colnames(pep) <- c("DTR1", "DTR2","DTR3","DTR4")

node_attrs <- tibble::tibble(name = c("DTR1", "DTR2","DTR3","DTR4"),
                             mean_est = summ7[17,])
graph <- tidygraph::as_tbl_graph(pep, directed = FALSE) %>%
  tidygraph::activate("nodes") %>% tidygraph::left_join(node_attrs,
                                                        by = "name")
legend_name <- c(mean_est = "Mean Trt\nSuccess")
graph <- graph %>% tidygraph::activate("edges")

p6r <- ggraph::ggraph(graph, layout = layout, weights = .data$weight) +
  ggraph::geom_edge_link(color = "gray", alpha = 0.6, aes(width = .data$weight)) +
  ggraph::geom_node_point(aes(color = .data[[color_by]]),
                          size = 7) + ggraph::geom_node_label(aes(label = .data$name),
                                                              nudge_y = 0.15, label.size = NA, hjust = "inward", alpha = 0.6) +
  ggplot2::scale_color_viridis_c(option = "plasma", name = legend_name[color_by]) +
  ggraph::scale_edge_width_continuous(name = "PEP", range=c(0.5,1), breaks=seq(0, 1, 0.02)) #can add range=c(0,10) to play with width



#put them in a grid
grid.arrange(p1r, p2r, p3r, p4r, p5r, p6r, nrow = 3)
