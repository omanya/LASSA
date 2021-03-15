#FUNCTION etilde
#   return the tau-expectile estimate with given partition
#   Y: data matrix p times n. p is the dimension, n number of observations
#   lab: label vector n times 1. Specifies the partitions. 1 for positive, -1 for negative.
#   alpha: -0.5 < alpha < 0.5. alpha + 0.5 is the expectile level tau.
#output: 

etilde <- function(Y, lab, alpha){
    p = dim(Y)[1]
    nplus <- sum(lab == 1)
    nminus <- sum(lab == -1)
    if(p > 1){
        sumplus <- apply(as.matrix(Y[,lab == 1]), 1, sum)
        summinus <- apply(as.matrix(Y[,lab == -1]), 1, sum)
    } else{
        sumplus <- sum(Y[,lab == 1])
        summinus <- sum(Y[,lab == -1])
    }
    top <- (0.5 + alpha)*sumplus + (0.5-alpha)*summinus
    bottom <- (0.5+alpha)*nplus + (0.5-alpha)*nminus
    return(top/bottom)
}

#FUNCTION expectile. iteratively find the tau-expectile of a sequence of numbers in R

expectile <- function(Y, alpha,nurexp=TRUE){
    Y=matrix(Y,1,length(Y))
    n = dim(Y)[2]
    lab <- ifelse(c(1:n) > n/2, 1, -1)        
    change = TRUE
    while(change){
        change = FALSE
        ehat <- etilde(Y, lab, alpha)
        newlab <- ifelse(Y > ehat, 1, -1)
        change <- ifelse(sum(lab!=newlab)>0, TRUE, FALSE)
        lab <- newlab
    }
  	if (nurexp==FALSE){
  		out=lab
  	} else {
  		out=ehat
  	}
    return(out)
}

#FUNCTION: return the first eigenvector psi given the partition
psihat <- function(Y, lab, alpha){
    n = dim(Y)[2]
    m <- etilde(Y, lab, alpha)
    #compute the covariance matrix normalized by m instead of the mean
    Yplus <- Y[,lab == 1]
    Yminus <- Y[,lab == -1]    
    Cplus = (0.5+alpha)/n*(Yplus-m)%*%t(Yplus-m)
    Cminus = (0.5-alpha)/n*(Yminus-m)%*%t(Yminus-m)
    C = Cplus + Cminus
    return(eigen(C)$vectors[,1])
}

#FUNCTION compute tau-variance of the projection
tauvar<-function(Y,psi, alpha){
    obs <- t(psi)%*%Y
    ehat <- expectile(obs, alpha)
    sumplus<-sum((obs[obs>ehat]-ehat)^2)*(0.5 +alpha)
    summinus<-sum((obs[obs<=ehat]-ehat)^2)*(0.5-alpha)
    vt<-1/length(obs)*(sumplus+summinus)
    return(vt)
}

#FUNCTION: given a direction psi, compute the expectile of the data in this direction. Return the lab.
expdir <- function(Y, psi, alpha){
    #compute the inner product
    obs <- t(psi)%*%Y
    ehat <- expectile(obs, alpha)
    lab <- ifelse(obs > ehat, 1, -1)
    return(lab)
}

#FUNCTION: principalDirection main function. 
#find the principal direction by iteratively updating psi and the weights
prdir <- function(Y, alpha, reset.opt = "random", ini.opt = "mean", lab.ini = NA, iter.tol = 10, reset.tol = 50,nc=NA,silent=T){
    p = dim(Y)[1]
    n = dim(Y)[2]
    #initiate some label vectors
    if(ini.opt == "mean"){
        lab = reinit(n, opt = "mean", Y = Y,nc=nc)
    }else{
        if(is.na(lab.ini)){
            lab = reinit(n, opt = "random",nc=nc)
        }else{
            lab = lab.ini
        }
    }
    change = TRUE
    iter = 1
    reset = 0
	  conv=1
    while(change){
        change = FALSE
        iter = iter + 1
        psi <- psihat(Y, lab, alpha)        
        newlab <- expdir(Y, psi, alpha)
        differ <- sum(lab!=newlab)
        change <- ifelse(differ > 0, TRUE, FALSE)
        lab <- newlab        
        #can get stuck in local min. In this case, initialize 
        #to a random label
        if(iter > iter.tol){
            lab <- reinit(n, opt = reset.opt, lab = lab, alpha = alpha, Y = Y, psi = psi)
            iter = 1
            reset = reset+1
		        change=TRUE
            if(reset > reset.tol){
                if(!silent){print(paste("WARNING: exceed reset level for alpha = ", alpha))}
                conv=0 #algorithm did not converge
		            change=FALSE
		        }            
        }
    }
	  if(!silent){print(paste("iter = ", iter, "reset = ", reset))}
    return(list(psi,conv,differ,lab))
}

predir <- function(Y, alpha, iter.tol = 10, reset.tol = 10,silent=T){
  p = dim(Y)[1]
  n = dim(Y)[2]
  #initiate some label vectors
  lab = reinit(n, opt = 'random')
  #initialize vars
  change = TRUE
  iter = 1
  reset = 0
  conv=1
  while(change){
    change = FALSE
    iter = iter + 1
    psi <- psihat(Y, lab, alpha)        
    newlab <- expdir(Y, psi, alpha)
    differ <- sum(lab!=newlab)
    change <- ifelse(differ > 0, TRUE, FALSE)
    lab <- newlab        
    #can get stuck in local min. In this case, initialize 
    #to a random label
    if(iter > iter.tol){
      lab <- reinit(n, opt = 'random')
      iter = 1
      reset = reset+1
      change=TRUE
      if(reset > reset.tol){
        if(!silent){print(paste("WARNING: exceed reset level for alpha = ", alpha))}
        conv=0 #algorithm did not converge
        change=FALSE
      }            
    }
  }
  if(!silent){print(paste("iter = ", iter, "reset = ", reset))}
  return(list(psi,conv,differ,lab))
}

#FUNCTION: options for reset and initialize the labels
#random
#flip far away
#jiggle tau
reinit <- function(n, opt = "random", lab = NA, alpha = NA, Y = NA, psi = NA,nc=NA){
    if(opt == "random"){
        lab=ifelse(runif(n,0,1)>0.5, 1, -1)
    }
    if(opt == "far"){
        n = length(lab)
        lab=lab*ifelse(runif(n,0,1)>0.8, 1, -1)
    }
    if(opt == "jiggle"){
        delta = (0.5 - abs(alpha))/10
        alpha.new = alpha + runif(1,-delta, delta)
        lab=expdir(Y, psi, alpha.new)
    }
    if(opt == "mean"){
        C = cov(t(Y))
        pc1 = eigen(C)$vectors[,1]
        lab=expdir(Y, pc1, 0)
    }
  	#check if at least one element of lab has the opposite sign, if not change sign of one elem
  	if (all(lab==1) || all(lab==-1)){
  		ind <-sample.int(n, size = nc, replace = FALSE)
  		lab[ind]<-lab[ind]*(-1)
  	}
	  return(lab)
}


#FUNCTION: compute the first k principal components
#added pec.ini parameter for initial lab
pec.k.ini <- function(Y, alpha, nk=2, reset.opt = "random", ini.opt = "mean", lab.ini = NA, pec.ini=NA,iter.tol = 10, reset.tol = 50,silent=T){
    p = dim(Y)[1]
    n = dim(Y)[2]
    psi.mat = matrix(NA, nrow = p, ncol = nk)
    mu.vec = rep(NA, nk)
    resid.obs <- Y
    conv=1
    cconv=NULL
    for(k in 1:nk){
        print(paste("computing for k = ", k))
        if(!is.na(pec.ini)[1]){
          lab.ini<-expdir(Y,pec.ini[,k],alpha)
        }
        outp <- prdir(resid.obs, alpha, reset.opt = reset.opt, ini.opt = ini.opt, lab.ini = lab.ini, iter.tol = iter.tol, reset.tol = reset.tol,nc=nk,silent=silent)
        psi<-outp[[1]]; conv=outp[[2]]; lab=outp[[4]]; cconv=c(cconv,conv)
        psi.mat[,k] <- psi
        obs = t(psi)%*%resid.obs
        mu.vec[k] <- expectile(obs, alpha)
        resid.obs = resid.obs - as.matrix(psi)%*%obs	
    }
   # psi.mat = psi.mat[,rev(1:dim(psi.mat)[2])] #WHY?
    return(list(mu.vec, psi.mat, min(cconv),lab))
}


#FUNCTION: compute the first k principal components with restarts
pec.k <- function(Y, alpha, nk=2, iter.tol = 10, reset.tol = 10,silent=T, restarts=100){
  p = dim(Y)[1]
  n = dim(Y)[2]
  psi.mat = matrix(NA, nrow = p, ncol = nk)
  mu.vec = rep(NA, nk)
  resid.obs <- Y
  conv=1
  cconv=NULL
  for(k in 1:nk){
    print(paste("computing for k = ", k))
    vt<-0
    rs<-0
    while(rs<=restarts){
      outp <- predir(resid.obs, alpha, iter.tol = iter.tol, reset.tol = reset.tol,silent=silent)
      vt_new<-tauvar(resid.obs,outp[[1]],alpha)
      if(vt_new>vt){
          psi<-outp[[1]]; conv=outp[[2]]; lab=outp[[4]]; cconv=c(cconv,conv)
          vt<-vt_new
      }
     rs<-rs+1 
    }
    psi.mat[,k] <- psi
    obs = t(psi)%*%resid.obs
    mu.vec[k] <- expectile(obs, alpha)
    resid.obs = resid.obs - as.matrix(psi)%*%obs	
  }
  psi.mat = psi.mat[,rev(1:dim(psi.mat)[2])]
  return(list(mu.vec, psi.mat, min(cconv),lab))
}


#---- FUNCTION: compare basis
#This is to compare the TopDown output with BottomUp and PrincipalExpectile
compareBasis <- function
### Project basis.old onto the span of a fixed basis basis.fix.
(basis.fix, 
###Fixed basis
basis.old
###Basis to be re-expressed in the span of basis.fix
){
    qrobj <- qr.solve(basis.fix, basis.old)
    basis.new <- basis.old%*%solve(qrobj)
    return(basis.new)
}

