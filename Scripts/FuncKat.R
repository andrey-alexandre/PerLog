##############
# Predito MC #
##############
predi.kat<- function(beta.vet,rho,covar,resp, li=0.000000001, ls=0.999999999){
  n<-length(resp)
  p<-theta<-NULL
  # cat(length(covar[1,]), length(beta.vet), '\n')
  for(t in 1:n){
    theta[t]<-exp(sum(covar[t,]*beta.vet))/(1+exp(sum(covar[t,]*beta.vet)))
    # theta[t] <- 1/(1+exp(-sum(covar[t,]*beta.vet)))
  }
  # cat(theta)
  theta[theta<=li]<-li;theta[theta>=ls]<-ls
  p[1] <- theta[1]
  c1 <- theta[2:n]*(1-theta[2:n])
  c2 <- theta[1:(n-1)]/(1-theta[1:(n-1)])
  c3 <- (1-theta[1:(n-1)])/theta[1:(n-1)]
  
  # c1[c1==-Inf]<--999999999999999999999999;c1[c1==Inf]<-999999999999999999999999
  # c2[c2==-Inf]<--999999999999999999999999;c2[c2==Inf]<-999999999999999999999999
  # c3[c3==-Inf]<--999999999999999999999999;c3[c3==Inf]<-999999999999999999999999
  
  p[2:n] <- theta[2:n] - rho*sqrt(c1)*(sqrt(c2)-resp[1:(n-1)]*(sqrt(c2)+sqrt(c3)))
  p[p<=li]<-li;p[p>=ls]<-ls
  # return(list(p=p, theta=theta))
  return(p)
}


################################
# Funcao da logverossimilhanca #
################################
loglik.kat <- function(par,out){
  resp<-out$resp; covar<-out$covar
  li<-out$li; ls<-out$ls
  n<-length(resp); k<-dim(covar)[2]
  beta.vet <- par[1:k]
  # rho <- exp(par[k+1])
  rho <- par[k+1]
  
  p<-theta<-NULL
  
  for(t in 1:n){
    theta[t]<-exp(sum(covar[t,]*beta.vet))/(1+exp(sum(covar[t,]*beta.vet)))
    # theta[t] <- 1/(1+exp(-sum(covar[t,]*beta.vet)))
  }
  theta[theta<=li]<-li;theta[theta>=ls]<-ls
  
  p[1] <- theta[1]
  
  c1 <- theta[2:n]*(1-theta[2:n])
  c1aux <- theta[1:(n-1)]*(1-theta[1:(n-1)]); c1aux<-c(0, c1aux)
  c2 <- theta[1:(n-1)]/(1-theta[1:(n-1)])
  c3 <- (1-theta[1:(n-1)])/theta[1:(n-1)]
  c4 <- (theta[1:(n-1)]^(-.5))*((1-theta[1:(n-1)])^(-1.5))
  c5 <- (theta[1:(n-1)]^(-1.5))*((1-theta[1:(n-1)])^(-.5))
  
  # c1[c1==-Inf]<--999999999999999999999999;c1[c1==Inf]<-999999999999999999999999
  # c2[c2==-Inf]<--999999999999999999999999;c2[c2==Inf]<-999999999999999999999999
  # c3[c3==-Inf]<--999999999999999999999999;c3[c3==Inf]<-999999999999999999999999
  # c4[c4==-Inf]<--999999999999999999999999;c4[c4==Inf]<-999999999999999999999999
  # c5[c5==-Inf]<--999999999999999999999999;c5[c5==Inf]<-999999999999999999999999
  
  p[2:n] <- theta[2:n] - rho*sqrt(c1)*(sqrt(c2)-resp[1:(n-1)]*(sqrt(c2)+sqrt(c3)))
  
  p[p<=li]<-li;p[p>=ls]<-ls
  
  dpdvt0 <- 1-.5*rho*(c1^(-.5))*(1-2*theta[2:n])*(sqrt(c2)-resp[1:(n-1)]*(sqrt(c2)+sqrt(c3)))
  dpdvt1 <- .5*rho*sqrt(c1)*(c4-resp[1:(n-1)]*(c4-c5))
  dpdvt0 <- c(1,dpdvt0)
  dpdvt1 <- c(0,dpdvt1)
  
  dldp <- (resp-p)/(p*(1-p))
  dvt0db <- theta*(1-theta)*covar
  dvt1db <- c1aux*covar
  
  dpdb <- dpdvt0*dvt0db + dpdvt1*dvt1db
  dldb <- as.vector(dldp)*dpdb
  
  dpdrho <- -sqrt(c1)*(sqrt(c2)-resp[1:(n-1)]*(sqrt(c2)+sqrt(c3))); dpdrho <- c(0, dpdrho)
  dldrho <- dldp*dpdrho
  
  # gradient
  grad <- c(colSums(dldb), sum(dldrho))
  
  #hessian
  d2ldp2 <- -resp/(p^2)-(1-resp)/((1-p)^2)
  d2ldrho2 <- d2ldp2*(dpdrho^2)
  
  d2ldb2 <- matrix(0, k, k)
  for (t in 1:n) {
    d2ldb2 <- d2ldb2 + d2ldp2[t]*dpdb[t,]%*%t(dpdb[t,])
  }
  
  d2pdrhodvt0 <- -.5*(c1^(-.5))*(1-2*theta[2:n])*(sqrt(c2)-resp[1:(n-1)]*(sqrt(c2)+sqrt(c3)))
  d2pdrhodvt0 <- c(0,d2pdrhodvt0)
  d2pdrhodvt1 <- .5*sqrt(c1)*(c4-resp[1:(n-1)]*(c4-c5))
  d2pdrhodvt1 <- c(0,d2pdrhodvt1)
  d2pdbdrho <- d2pdrhodvt0*dvt0db + d2pdrhodvt1*dvt1db
  
  d2ldbdrho <- as.vector(d2ldp2)*dpdb*dpdrho + as.vector(dldp)*d2pdbdrho
  
  hess <- cbind(rbind(d2ldb2, colSums(d2ldbdrho)), c(colSums(d2ldbdrho), sum(d2ldrho2)))
  
  logver<- -sum(resp*log(p/(1-p))+log(1-p))
  attr(logver, 'gradient') <- -grad/n
  attr(logver, 'hessian') <- -hess
  return(logver)
}


gr.loglik.kat<-function(par,out) return(attr(loglik.kat(par,out), 'gradient'))
hess.loglik.kat<-function(par,out) return(attr(loglik.kat(par,out), 'hessian'))


logistic_MC.kat<-function(resp, covar, par0, trace=2, li=0.000000001, ls=0.999999999){
  out<-list(covar=covar, resp=resp, li=li, ls=ls)
  otim <- optim(par = par0, fn = loglik.kat, gr = gr.loglik.kat, out=out, method = "L-BFGS-B",
                lower = c(rep(-Inf, dim(covar)[2]), -0.999999),
                upper = c(rep(Inf, dim(covar)[2]), 0.999999), hessian = FALSE,
                control = list(trace=trace))
  # seotim<-sqrt(diag(solve(otim$hessian)))
  seotim <- sqrt(diag(solve(hess.loglik.kat(par = otim$par, out = out))))
  mat.res<-cbind(otim$par, seotim, 2*(1-pnorm(abs(otim$par/seotim))))
  colnames(mat.res)<-c('Estimate','Std.Error','p-value')
  # row.names(mat.res)[dim(covar)[2]+1]<-'lambda'
  row.names(mat.res)[1:dim(covar)[2]]<-colnames(covar)
  row.names(mat.res)[dim(covar)[2]+1]<-'rho'
  return(mat.res)
}

#####################
# Gerar o y de y~X  #
#####################
gen.kat <- function(X, beta, rho, li=0.000000001, ls=0.999999999){
  n <- dim(X)[1]
  s <- length(rho)
  theta <- p <- NULL
  
  theta <- exp(X%*%beta)/(1 + exp(X%*%beta))
  theta[theta<=li] <- li
  theta[theta>=ls] <- ls
  p[1] <- theta[1]
  y <- rbinom(1, 1, p[1])
  for (t in 2:n) {
    t_m <- ifelse(t%%s == 0, s, t%%s)
    c1  <-  theta[t]*(1-theta[t])
    c2  <-  theta[t-1]/(1-theta[t-1])
    if(y[t-1]==0) p[t]  <-  theta[t] - rho[t_m]*sqrt(c1*c2) else p[t]  <-  theta[t] + rho[t_m]*sqrt(c1/c2)
    # p[t]  <-  theta[t] - rho*sqrt(c1*c2) + y[t-1]*rho*sqrt(c1)*(sqrt(1/c2)+sqrt(c2))
    if(p[t]>1) p[t] <- 1
    if(p[t]<0) p[t] <- 0
    y[t] <- rbinom(n = 1, size = 1, prob = p[t])
  }
  return(y)
}

#################
# PFP minimo MC #
#################

min.pfp.MC.kat <- function(corte0=seq(0.01,0.99,by=0.01), corte1=seq(0.01,0.99,by=0.01),
                       pfn.target=.05, y, prob, li=0.000000001, ls=0.999999999){
  # INICIO FUNCAO
  sens <- espec <- pfp <- pfn <- acc <- matrix(data = NA, nrow = length(corte0), ncol = length(corte1))
  n <- length(y)
  nenhum.menor.que.pfn.target=TRUE;pfp.min <- 1
  for(l in 1:length(corte0)){
    for(j in 1:length(corte1)){
      classif <-  rep(0,n)
      for(i in 2:n){
        if(y[i-1]==0){
          if(prob[i]>corte0[l]){classif[i] <- 1}
        }else{
          if(prob[i]>corte1[j]){classif[i] <- 1}
        }
      }
      A <- table(y[2:n], classif[2:n]) 
      if(ncol(A)==1){                                                 
        sens[l,j] <- classif[2]                                                      
        espec[l,j] <-  1-sens[l,j]                                                     
        pfp[l,j] <- A[1,]/sum(A[,1])
        pfn[l,j] <- A[2,]/sum(A[,1])
        if(classif[2] == 0) acc[l,j] <- A[1,1]/sum(A)
        if(classif[2] == 1) acc[l,j] <- A[2,1]/ sum(A)
      }else{                                                                   
        espec[l,j] <-  A[1,1]/sum(A[1,])
        sens[l,j] <-  A[2,2]/sum(A[2,])
        pfp[l,j] <- A[1,2]/sum(A[,2])
        pfn[l,j] <- A[2,1]/sum(A[,1])
        acc[l,j] <- sum(diag(A))/sum(A)
      }
      # if((l==1)&(j==1)){pfp.min <- pfp[l,j]}
      if( (pfn[l,j]<=pfn.target) & (pfp[l,j]<pfp.min) ){
        corte0min <- corte0[l]
        corte1min <- corte1[j]
        pfp.min <- pfp[l,j]
        pfn.min <- pfn[l,j]
        acc.min <- acc[l,j]
        nenhum.menor.que.pfn.target=FALSE
      }
    }
  }
  if(nenhum.menor.que.pfn.target){
    ind <- which(pfn == min(pfn), arr.ind = TRUE)[1,]; l <- ind[1]; j <- ind[2]
    pfn.min <- pfn[l,j];pfp.min <- pfp[l,j]; acc.min <- acc[l,j]
    corte0min <- corte0[l];corte1min <- corte1[j]
    warning(paste('Nao ha PFN < ',pfn.target,' (target)!\n Escolhendo PFP com PFN minimo!', sep=''))
  }
  return(list(pfp.min=pfp.min, pfn.target=pfn.target, pfn.min=pfn.min, acc.min = acc.min,
              corte0min=corte0min, corte1min=corte1min))
}

##################
# PFP minimo glm #
##################

min.pfp.glm.kat <- function(corte=seq(0.01,0.99,by=0.01), pfn.target=.05, y, prob,
                        li=0.000000001, ls=0.999999999){
  sens <- espec <- pfp <- pfn <- acc <- NULL
  pfn.min <- NULL
  nenhum.menor.que.pfn.target=TRUE;pfp.min <- 1
  for(j in 1:length(corte)){
    classif <-  rep(0,(length(y)))
    for(i in 1:length(y)){ if(prob[i]>corte[j]) classif[i] <- 1 }
    A <- table(y,classif)
    if(ncol(A)==1){                                                 
      sens[j] <- classif[2]                                                      
      espec[j] <-  1-sens[j]                                                     
      pfp[j] <- A[1,]/sum(A[,1])
      pfn[j] <- A[2,]/sum(A[,1])
      if(classif[2] == 0) acc[j] <- A[1,1]/sum(A)
      if(classif[2] == 1) acc[j] <- A[2,1]/ sum(A)
    }else{                                                                   
      espec[j] <-  A[1,1]/sum(A[1,])
      sens[j] <-  A[2,2]/sum(A[2,])
      pfp[j] <- A[1,2]/sum(A[,2])
      pfn[j] <- A[2,1]/sum(A[,1])
      acc[j] <- sum(diag(A))/sum(A)
    }
    if( (pfn[j]<=pfn.target) & (pfp[j]<pfp.min) ){
      cortemin <- corte[j]
      pfp.min <- pfp[j]
      pfn.min <- pfn[j]
      acc.min <- acc[j]
      nenhum.menor.que.pfn.target=FALSE
    }
  }
  if(nenhum.menor.que.pfn.target){
    pfn.min <- min(pfn)
    cortemin <- corte[which.min(pfn)]
    pfp.min <- pfp[which.min(pfn)]
    acc.min <- acc[which.min(pfn)]
    warning(paste('Nao ha PFN < ',pfn.target,' (target)!\n Escolhendo PFP com PFN minimo!', sep=''))
  }
  return(list(pfp.min=pfp.min, pfn.target=pfn.target, pfn.min=pfn.min, acc.min = acc.min,
              cortemin=cortemin))
}

############
# PrevisÃ£o #
############

#prever<- m-n
#Xprev<- X1[n:m,]
prev.glm.kat <- function(fit, newY=NULL, newX, corte, li=0.000000001, ls=0.999999999){#Only new covariates!!
  prob.glm.p  <- predi.katct(object = fit, newx = newX, type=c("response"))
  n.desc <- dim(newX)[1]
  classif.glm <-  rep(0,n.desc)
  for(i in 1:n.desc){if(prob.glm.p[i]>corte) classif.glm[i] <- 1 }
  if(!is.null(newY)){
    A <- table(newY,classif.glm)
    if(ncol(A)==1){
      if(classif.glm[1]==1) sens.glm <- 1; acc.glm <- A[2,1]/sum(A)
      if(classif.glm[1]==0) sens.glm <- 0; acc.glm <- A[1,1]/ sum(A)
      espec.glm <- 1-sens.glm
      
    }else{
      espec.glm <-  A[1,1]/sum(A[1,])
      sens.glm <-  A[2,2]/sum(A[2,])
      acc.glm <- sum(diag(A))/sum(A)
    }
    if(sens.glm==0){sens.glm <- li}
    if(sens.glm==1){sens.glm <- ls}
    if(espec.glm==0){espec.glm <- li}
    if(espec.glm==1){espec.glm <- ls}
    p <- mean(newY)
    phat.glm <- mean(classif.glm)
    pfp.glm <- 1-(p*sens.glm)/(p*sens.glm+(1-p)*(1-espec.glm))
    pfn.glm <- 1-(espec.glm*(1-p))/((1-sens.glm)*p+espec.glm*(1-p))
    return(list(prev=classif.glm, pfp=pfp.glm, pfn=pfn.glm, acc = acc.glm))
  }
  else{
    return(list(prev=classif.glm))
  }
}

prev.MC.kat<-function(beta, rho, newY, newX, LastY, LastX, corte0, corte1){
  y.desc <- c(LastY, newY)
  X.desc <- rbind(LastX, newX)
  n.desc <- length(y.desc)
  prob.p <- predi.kat(beta.vet=beta, rho = rho, covar=X.desc, resp=y.desc)
  classif.MC<- rep(0,(n.desc-1))
  for(i in 2:(n.desc-1)){
    if(y.desc[i-1]==0){
      if(prob.p[i]>corte0){classif.MC[i-1]<-1}
    }else{
      if(prob.p[i]>corte1){classif.MC[i-1]<-1}
    }
  }
  A<-table(newY,classif.MC) 
  if(ncol(A)==1){                                                 
    sens.MC<-classif.MC[1]                                                      
    espec.MC<- 1-sens.MC
    if(classif[2] == 0) acc.MC <- A[1,1]/sum(A)
    if(classif[2] == 1) acc.MC <- A[2,1]/ sum(A)                                                  
  }else{
    # cortemin<-corte
    espec.MC<- A[1,1]/sum(A[1,])                                             
    sens.MC<- A[2,2]/sum(A[2,])
    acc.MC <- sum(diag(A))/sum(A)                                              
  }
  p<-mean(newY)
  phat.MC<-mean(classif.MC)
  pfp.MC<-1-(p*sens.MC)/(p*sens.MC+(1-p)*(1-espec.MC))
  pfn.MC<-1-(espec.MC*(1-p))/((1-sens.MC)*p+espec.MC*(1-p))
  return(list(prev=classif.MC, pfp=pfp.MC, pfn=pfn.MC, acc = acc.MC))
}

