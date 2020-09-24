#################
### LIBRARIES ###
#################
library(tidyverse)
library(lbfgs)
##############
# Predito MC #
##############
predi <-  function(beta.vet,rho,covar,resp, li=0.000000001, ls=0.999999999){
  n <- length(resp)
  s <- length(rho)
  per <- ifelse(2:n%%s == 0, s, 2:n%%s)
  p <- theta <- NULL
  # cat(length(covar[1,]), length(beta.vet), '\n')
  for(t in 1:n){
    theta[t] <- exp(sum(covar[t,]*beta.vet))/(1+exp(sum(covar[t,]*beta.vet)))
    # theta[t]  <-  1/(1+exp(-sum(covar[t,]*beta.vet)))
  }
  # cat(theta)
  theta[theta<=li] <- li;theta[theta>=ls] <- ls
  p[1]  <-  theta[1]
  c1  <-  theta[2:n]*(1-theta[2:n])
  c2  <-  theta[2:n-1]/(1-theta[2:n-1])
  c3  <-  (1-theta[2:n-1])/theta[2:n-1]
  
  # c1[c1==-Inf] <- -999999999999999999999999;c1[c1==Inf] <- 999999999999999999999999
  # c2[c2==-Inf] <- -999999999999999999999999;c2[c2==Inf] <- 999999999999999999999999
  # c3[c3==-Inf] <- -999999999999999999999999;c3[c3==Inf] <- 999999999999999999999999
  
  p[2:n]  <-  theta[2:n] - rho[per]*sqrt(c1)*(sqrt(c2)-resp[2:n-1]*(sqrt(c2)+sqrt(c3)))
  p[p<=li] <- li;p[p>=ls] <- ls
  # return(list(p=p, theta=theta))
  return(p)
}

#######################################
# Cria a matriz de derivada periódica #
#######################################
permatrix <- function(vector, period_vec, period){
  matrix <- data.frame(A = vector, per = c(1, period_vec)) %>% 
    mutate(ID = row_number())%>% 
    tidyr::spread(key = per, value = 'A') %>% select(-ID) %>% 
    select_at(.vars = 1:period, .funs = function(x){paste0('Per', x)}) %>% 
    mutate_all(.funs = function(x) tidyr::replace_na(x, replace = 0)) %>% as.matrix
  return(matrix)
}

################################
# Funcao da logverossimilhanca #
################################
loglik  <-  function(par_,out,w=NULL){
  resp <- out$resp; covar <- out$covar
  li <- out$li; ls <- out$ls;
  lamb <- out$lamb; alpha <- out$alpha;
  par <- par_/alpha
  n <- length(resp); k <- dim(covar)[2]
  beta  <-  par[1:k]
  # rho  <-  exp(par[k+1])
  rho  <-  par[(k+1):length(par)]
  s <- length(rho)
  per <- ifelse(2:n%%s == 0, s, 2:n%%s)
  if(is.null(w)) w <- rep(1, n-1) else w <- w[2:n]
  p <- NULL
  
  theta <- exp(covar%*%beta)/(1+exp(covar%*%beta)); theta <- as.vector(theta)
  theta[theta<=li] <- li;theta[theta>=ls] <- ls
  
  p[1]  <-  theta[1]
  
  c1  <-  theta[2:n]*(1-theta[2:n])
  c1aux1  <-  theta[2:n-1]*(1-theta[2:n-1]); c1aux <- c(0, c1aux1)
  c2  <-  theta[2:n-1]/(1-theta[2:n-1])
  c3  <-  (1-theta[2:n-1])/theta[2:n-1]
  c4  <-  (theta[2:n-1]^(-.5))*((1-theta[2:n-1])^(-1.5))
  c5  <-  (theta[2:n-1]^(-1.5))*((1-theta[2:n-1])^(-.5))
  
  # c1[c1==-Inf] <- -999999999999999999999999;c1[c1==Inf] <- 999999999999999999999999
  # c2[c2==-Inf] <- -999999999999999999999999;c2[c2==Inf] <- 999999999999999999999999
  # c3[c3==-Inf] <- -999999999999999999999999;c3[c3==Inf] <- 999999999999999999999999
  # c4[c4==-Inf] <- -999999999999999999999999;c4[c4==Inf] <- 999999999999999999999999
  # c5[c5==-Inf] <- -999999999999999999999999;c5[c5==Inf] <- 999999999999999999999999
  A <- sqrt(c1*c1aux[-1])
  B <- sqrt(c1/c1aux[-1])
  # p[2:n] <- theta[2:n] + rho[per]*B*(resp[2:n-1] - theta[2:n-1])
  
  p[2:n] <- theta[2:n] + rho[per]*sqrt(c1)*(resp[2:n-1]-theta[2:n-1])/sqrt(c1aux1)#p_{y_{t-1}}, t>=2
  
  # p[2:n]  <-  theta[2:n] - rho[per]*sqrt(c1)*(sqrt(c2)-resp[2:n-1]*(sqrt(c2)+sqrt(c3)))
  
  p[p<=li] <- li;p[p>=ls] <- ls
  
  dldp <- ((resp - p)/(p*(1-p)))[2:n]
  dldp_aux <- ((resp[1] - p[1])/(p[1]*(1-p[1])))
  
  dBdv0 <- .5*(1-2*theta[2:n])/A
  dBdv1 <- -.5*(1-2*theta[2:n-1])*B/c1aux[-1]
  # 
  # dpdv0 <- 1 + rho[per]*(resp[2:n-1]-theta[2:n-1])*dBdv0
  # dpdv1 <- rho[per]*(dBdv1*(resp[2:n-1]-theta[2:n-1])-B)
  
  dpdv0 <- 1+.5*rho[per]*(1-2*theta[2:n])*(resp[2:n-1]-theta[2:n-1])/sqrt(c1*c1aux1)
  dpdv1 <- -.5*rho[per]*sqrt(c1)*(resp[2:n-1]-theta[2:n-1])*(1-2*theta[2:n-1])/(c1aux1^1.5)
  
  dv0db <- c1*covar[2:n,]
  dv0db_aux <- theta[1]*(1-theta[1])*covar[1,,drop=FALSE]
  
  dv1db <- c1aux1*covar[2:n-1,]
  
  dpdb <- dpdv0*dv0db + dpdv1*dv1db #b = beta, B = B mesmo
  
  dldb <- dldp*dpdb
  
  mat.aux <- NULL
  for(v in 1:s) mat.aux <- cbind(mat.aux, ifelse(per==v, 1, 0))
  dpdr <- (resp[2:n-1]-theta[2:n-1])*sqrt(c1)*mat.aux/sqrt(c1aux1)
  
  dldr <- dldp*dpdr
  
  # #dpdvt0  <-  1-.5*rho[per]*(c1^(-.5))*(1-2*theta[2:n])*(sqrt(c2)-resp[2:n-1]*(sqrt(c2)+sqrt(c3))) # Katiuski
  # dpdvt0  <-  1-.5 * rho[per] * (c2^(-.5)) * (c1^(-.5)) * (resp[2:n-1] + c2 * (resp[2:n-1]-1)) # Andrey
  # #dpdvt1  <-  .5*rho[per]*sqrt(c1)*(c4-resp[2:n-1]*(c4-c5)) # Katiuski
  # dpdvt1  <-  .5 * rho[per] * resp[2:n-1] * (c1^(.5)) * ( ( (2 * theta[2:n-1]-1) * sqrt(c1aux[-1])  ) / (c1aux[-1]^2) ) - .5 * rho[per] * (c1^(-.5)) * ( theta[2:n] / ((1-theta[2:n-1])) )
  # dpdvt0  <-  c(1,dpdvt0)
  # dpdvt1  <-  c(0,dpdvt1)
  # 
  # dldp  <-  (resp-p)/(p*(1-p))
  # dvt0db  <-  theta*(1-theta)*covar
  # dvt1db  <-  c1aux*apply(X = covar, MARGIN = 2, FUN = lag);  dvt1db[1,] <- rep(0, k)  
  # 
  # dpdb  <-  dpdvt0*dvt0db + dpdvt1*dvt1db
  # dldb  <-  as.vector(dldp)*dpdb
  # 
  # #dpdrho  <-  -sqrt(c1)*(sqrt(c2)-resp[2:n-1]*(sqrt(c2)+sqrt(c3))); dpdrho  <-  c(0, dpdrho) # Katiuski
  # dpdrho  <- sqrt(c1) * ( (resp[2:n-1]-1) * c2 + resp[2:n-1] * c3 ); dpdrho  <-  c(0, dpdrho) # Andrey
  # dpdrho <- permatrix(vector = dpdrho, period_vec = per, period = s)
  # dldrho  <-  dldp*dpdrho
  
  # del <- 0.5
  # dpendb <- -lamb*alpha*as.numeric(beta < -del) +
  #   lamb*alpha*as.numeric(beta > del) +
  #   lamb*alpha*beta*as.numeric(beta <= del & beta >= -del)/del
  # # dpendb <- 2*lamb*alpha*beta
  # dpendrho <- rep(0, s)
  # gradient
  grad  <-  c(colMeans(w*dldb), colMeans(w*dldr))/alpha
  
  #hessian
  # d2ldp2 <-  - ( 1 - p - p ^2) / ( (1-p)^2 * p^2 ); d2ldp2 <- d2ldp2[2:n-1]
  d2ldp2 <- (resp[2:n]*(2*p[2:n] - 1) - p[2:n]^2)/(p[2:n]^2*(1 - p[2:n])^2)
  d2ldp2_aux <- (resp[1]*(2*p[1] - 1) - p[1]^2)/(p[1]^2*(1 - p[1])^2)
  
  d2ldr2 <- array(data = 0, dim = c(s, s, n-1))
  for(i in 2:n-1){
    d2ldr2[, , i] <- w[i] * d2ldp2[i] * t(dpdr[i, , drop = F]) %*% dpdr[i, , drop = F]
  }
  
  d2Bdv02 <- - A^(-1) * ( 1 + .25 * (c1^(-1)) * (1 - 2 * theta[2:n]) ^ 2 )
  d2pdv02 <- rho[per]*(resp[2:n-1]-theta[2:n-1])*d2Bdv02
  
  d2Bdv12 <- B * (c1aux[-1] ^ -1) *  (1 + .75 * (c1aux[-1] ^ -1) * (1 - 2 * theta[2:n-1]) ^ 2 )  
  d2pdv12 <- rho[per]*( d2Bdv12 * (resp[2:n-1]-theta[2:n-1]) - 2 * dBdv1 )
  
  d2Bdv0dv1 <- -.25 * (A^(-1)) * (c1aux[-1]^(-1))*(1 - 2 * theta[2:n]) * (1 - 2 * theta[2:n-1])
  d2pdv0dv1 <- rho[per]*( d2Bdv0dv1 * (resp[2:n-1]-theta[2:n-1]) - dBdv0 )
  
  
  d2ldb2 <- C <-  array(data = 0, dim = c(k, k, n-1))
  for(i in 2:n-1){
    d2ldb2[, , i] <-
      w[i] * d2ldp2[i] * t(dpdb[i, , drop = F]) %*% dpdb[i, , drop = F] + #(I)
      w[i] * dldp[i] * ( #(II)
        (t(dv0db[i, , drop = F]) %*% dv0db[i, , drop = F]) * d2pdv02[i] +
          (t(dv0db[i, , drop = F]) %*% dv1db[i, , drop = F]) * d2pdv0dv1[i] +
          (t(dv1db[i, , drop = F]) %*% dv0db[i, , drop = F]) * d2pdv0dv1[i] +
          (t(dv1db[i, , drop = F]) %*% dv1db[i, , drop = F]) * d2pdv12[i] +
          dpdv0[i] * theta[i+1]*(1-theta[i+1])*(1-2*theta[i+1])*(t(covar[i+1,,drop=FALSE])%*%covar[i+1,,drop=FALSE]) +
          dpdv1[i] * theta[i]*(1-theta[i])*(1-2*theta[i])*(t(covar[i,,drop=FALSE])%*%covar[i,,drop=FALSE])
      )
  }
  
  d2ldbdr <- array(data = 0, dim = c(k, s, n-1))
  for(i in 2:n-1){
    # d2ldbdr[, , i] <- t( (resp[2:n-1]-theta[2:n-1])[i] * dBdv0[i] * mat.aux[i, , drop = F] ) %*% dv0db[i, , drop = F] + 
    # t( ( (resp[2:n-1]-theta[2:n-1]) * dBdv1 - B )[i] * mat.aux[i, , drop = F]  )%*% dv1db[i, , drop = F]
    d2ldbdr[, , i] <- w[i] * d2ldp2[i] *(t(dpdb[i, , drop = F]) %*% dpdr[i, , drop = F]) +
      w[i] * dldp[i] * (
        (resp[i] - theta[i])*dBdv0[i]*(t(dv0db[i,,drop=FALSE]) %*% mat.aux[i,,drop=FALSE]) +
          ((resp[i] - theta[i])*dBdv1[i]-B[i])*(t(dv1db[i,,drop=FALSE]) %*% mat.aux[i,,drop=FALSE])
      )
  }
  
  
  # #d2ldp2  <-  -resp/(p^2)-(1-resp)/((1-p)^2) # Katiuski
  # d2ldp2  <-  - ( p^2 - 2 * p * resp ) / ( (1-p)^2 * p^2 ) # Andrey
  # 
  # d2ldrho2  <-  array(0, dim = c(s, s, n))
  # for (t in 1:n) {
  #   d2ldrho2[ , ,t]  <-  d2ldp2[t]*dpdrho[t,]%*%t(dpdrho[t,])
  # }
  # 
  # d2ldb2  <- array(0, dim = c(k, k, n))
  # for (t in 1:n) {
  #   d2ldb2[ , , t]  <- d2ldp2[t]*dpdb[t,]%*%t(dpdb[t,])
  # }
  # 
  # # d2pdrhodvt0  <-  -.5*(c1^(-.5))*(1-2*theta[2:n])*(sqrt(c2)-resp[2:n-1]*(sqrt(c2)+sqrt(c3))); d2pdrhodvt0  <-  c(0,d2pdrhodvt0); # Katiuski
  # d2pdrhodvt0 <- .5 * ( (1 - 2 *  theta[2:n]) / sqrt(c1) ) * ( resp[2:n-1] * ( sqrt(c2) + sqrt(c3) ) -  theta[2:n-1] / sqrt(c1aux[-1]) ); d2pdrhodvt0 <- c(0, d2pdrhodvt0); # Andrey
  # d2pdrhodvt0 <- permatrix(vector = d2pdrhodvt0, period_vec = per, period = s)
  # # d2pdrhodvt1  <-  .5*sqrt(c1)*(c4-resp[2:n-1]*(c4-c5)); d2pdrhodvt1  <-  c(0,d2pdrhodvt1); # Katiuski
  # d2pdrhodvt1  <-  .5 * ( sqrt(c1) / c1aux[-1] ) * ( resp[2:n-1] * sqrt(c3) + (1 - resp[2:n-1]) * sqrt(c2) ); d2pdrhodvt1  <-  c(0,d2pdrhodvt1); # Katiuski
  # d2pdrhodvt1 <- permatrix(vector = d2pdrhodvt1, period_vec = per, period = s)
  # 
  # d2pdbdrho  <- array(0, dim = c(k, s, n))
  # for (t in 1:n) {
  #   d2pdbdrho[ , , t]  <- t(dvt0db)[,t, drop = F]%*%d2pdrhodvt0[t, , drop = F] + t(dvt1db)[,t, drop = F]%*% d2pdrhodvt1[t, , drop = F]
  # }
  # 
  # d2ldbdrho  <-  array(0, dim = c(k, s, n))
  # for (t in 1:n) {
  #   d2ldbdrho[ , , t]  <- as.vector(d2ldp2)[t] * t(dpdb)[, t, drop = F] %*%  dpdrho[t, , drop = F] + as.vector(dldp)[t]*d2pdbdrho[, , t]
  # }
  # 
  hess  <-  cbind(rbind((rowSums(x = d2ldb2, dims = 2)), t(rowSums(x = d2ldbdr, dims = 2))),
                  rbind( rowSums(x = d2ldbdr, dims = 2), rowSums(x = d2ldr2, dims = 2)))
  # 
  logver <-  -mean(w*resp[2:n]*log(p[2:n])+w*(1-resp[2:n])*log(1-p[2:n]))
  # logver <-  -mean(w*resp[2:n]*log(p[2:n])+w*(1-resp[2:n])*log(1-p[2:n])) + lamb*(sum(alpha * beta^2))
  # logver <-  sum(w*resp[2:n]*log(p[2:n])+w*(1-resp[2:n])*log(1-p[2:n])) - lamb*(sum(alpha * abs(beta)))
  # # cat(logver,'\n')
  attr(logver, 'gradient')  <-  -grad
  attr(logver, 'hessian')  <-  -hess
  return(logver)
}


gr.loglik <- function(par,out,w=NULL) return(attr(loglik(par,out,w), 'gradient'))
hess.loglik <- function(par,out,w=NULL) return(attr(loglik(par,out,w), 'hessian'))

logistic_MC <- function(resp,covar,par0, lamb=0, alpha=1, w0=NULL,trace=0,li=0.001,ls=0.999, maxiter=500){
  out <- list(covar=covar, resp=resp, li=li, ls=ls, lamb=lamb, alpha=alpha)
  if(lamb==0){
    otim  <-  optim(par = par0, fn = loglik, gr = gr.loglik, out=out, w=w0, method = "L-BFGS-B",
                    lower = c(rep(-Inf, dim(covar)[2]), rep(-0.999999, length(par0) - dim(covar)[2])),
                    upper = c(rep(Inf, dim(covar)[2]), rep(0.999999, length(par0) - dim(covar)[2])),
                    hessian = FALSE, control = list(trace=trace))
    seotim <- diag(solve(hess.loglik(par = otim$par, out = out, w=w0)))
    mat.res <- cbind(otim$par, seotim, 2*(1-pnorm(abs(otim$par/seotim))))
    colnames(mat.res) <- c('Estimate','Std.Error','p-value')
  }else{
    par00 <- par0*alpha
    otim <- lbfgs(call_eval = loglik, call_grad = gr.loglik, vars = par00, out=out, w=w0,
                  invisible = as.numeric(trace==0),
                  orthantwise_c = lamb#,
                  # orthantwise_start = 0,
                  # orthantwise_end = dim(covar)[2]
                  )
    mat.res <- cbind(otim$par/alpha)
    colnames(mat.res) <- c('Estimate')
  }
  # otim$par <- ifelse(otim$par < 10^-5, 0, otim$par)
  # seotim <- sqrt(diag(solve(otim$hessian)))
  # seotim  <- sqrt(ifelse(diag(solve(hess.loglik(par = otim$par, out = out, w=w0))) <= 0, li, diag(solve(hess.loglik(par = otim$par, out = out, w=w0)))))
  row.names(mat.res) <- c(colnames(covar), paste0('rho', seq(length(par0)-dim(covar)[2])))
  return(mat.res)
}

#####################
# Gerar o y de y~X  #
#####################
gen <- function(X, beta, rho, li=0.000000001, ls=0.999999999){
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
pfp.pfn.MC <- function(corte, out){
  y <- out$y; prob <- out$prob; s <- out$s
  N <- length(y); R <- N%/%s
  c0 <- corte[1:s]; c1 <- corte[s+1:s]
  yp <- NULL
  for(r in 1:R){
    for(v in 1:s) {
      t <- (r-1)*s+v
      if(t>1){
        if(y[t-1] == 0 & prob[t] <= c0[v]) yp[t] <- 0
        if(y[t-1] == 0 & prob[t] > c0[v]) yp[t] <- 1
        if(y[t-1] == 1 & prob[t] <= c1[v]) yp[t] <- 0
        if(y[t-1] == 1 & prob[t] > c1[v]) yp[t] <- 1
      }
    }
  }
  a <- length(which(y[2:N]==0 & yp[2:N]==0))
  b <- length(which(y[2:N]==0 & yp[2:N]==1))
  c <- length(which(y[2:N]==1 & yp[2:N]==0))
  d <- length(which(y[2:N]==1 & yp[2:N]==1))
  acc <- (a+d)/(a+b+c+d)
  esp <- ifelse(a+b==0, (a+c)/(a+b+c+d), a/(a+b))
  sen <- ifelse(c+d==0, (b+d)/(a+b+c+d), d/(c+d))
  pfp <- ifelse(b+d==0, (a+b)/(a+b+c+d), b/(b+d))
  pfn <- ifelse(a+c==0, (c+d)/(a+b+c+d), c/(a+c))
  return(list(pfp=pfp, pfn=pfn, sen=sen, esp=esp, acc=acc))
}

min.pfp.MC.func <- function(corte, out){
  out1 <- list(y=out$y, prob=out$prob, s=out$s)
  lamb <- out$lamb; esp.target <- out$esp.target
  res <- pfp.pfn.MC(corte = corte, out = out1)
  func <- (1-res$sen) + lamb*max(0, (1-res$esp-(1-esp.target)))
  return(func)
}

min.pfp.MC <- function(corte.inic=NULL, lamb=30, esp.target=.95,
                       y, prob, s, maxiter=500){
  # INICIO FUNCAO
  if(is.null(corte.inic) |
     (!is.null(corte.inic) & length(corte.inic)!= 2*s))
    corte.inic <- rep(.5, 2*s)
  out <- list(y=y, prob=prob, s=s, lamb=lamb, esp.target=esp.target)
  fit <- GA::ga(type = 'real-valued',
                fitness = min.pfp.MC.func,
                out=out,
                maxiter = maxiter,
                monitor = FALSE,
                lower = rep(0.01, 2*s), upper = rep(0.99, 2*s)
  )
  
  c0 <- attr(fit, 'solution')[1,1:s]; c1 <- attr(fit, 'solution')[1, s+1:s]
  names(c0) <- names(c1) <- paste('v=',1:s,sep='')
  val <- pfp.pfn.MC(corte = c(c0, c1),
                    out = list(y=y, prob=prob, s=s))
  
  if(1 - val$esp > 1 - esp.target)
    warning(paste('PFN <= ',esp.target,' (target) não encontrado!\n Escolhendo PFP com PFN minimo!', sep=''))
  
  return(list(pfp=val$pfp,
              pfn=val$pfn,
              acc=val$acc,
              esp=val$esp,
              sen=val$sen,
              c0=c0,
              c1=c1,
              attain.esp.target=(1 - val$esp <= 1 - esp.target),
              esp.target=esp.target))
}



##################
# PFP minimo glm #
##################
pfp.pfn.glm <- function(corte, out){
  y <- out$y; prob <- out$prob
  
  yp <- as.numeric(prob > corte)
  
  a <- length(which(y==0 & yp==0))
  b <- length(which(y==0 & yp==1))
  c <- length(which(y==1 & yp==0))
  d <- length(which(y==1 & yp==1))
  acc <- (a+d)/(a+b+c+d)
  esp <- ifelse(a+b==0, (a+c)/(a+b+c+d), a/(a+b))
  sen <- ifelse(c+d==0, (b+d)/(a+b+c+d), d/(c+d))
  pfp <- ifelse(b+d==0, (a+b)/(a+b+c+d), b/(b+d))
  pfn <- ifelse(a+c==0, (c+d)/(a+b+c+d), c/(a+c))
  
  
  return(list(pfp=pfp, pfn=pfn, sen=sen, esp=esp, acc=acc))
}

min.pfp.glm.func <- function(corte, out){
  out1 <- list(y=out$y, prob=out$prob)
  lamb <- out$lamb; esp.target <- out$esp.target
  res <- pfp.pfn.glm(corte = corte, out = out1)
  func <- (1-res$sen) + lamb*max(0, (1-res$esp-(1-esp.target)))
  return(func)
}

min.pfp.glm <- function(corte.inic=NULL, lamb=30, esp.target=.95,
                        y, prob, maxiter=500){
  # INICIO FUNCAO
  if(missing(corte.inic)) corte.inic <- .5
  out <- list(y=y, prob=prob, lamb=lamb, esp.target=esp.target)
  fit <- GA::ga(type = 'real-valued',
                fitness = min.pfp.glm.func,
                out=out,
                maxiter = maxiter,
                monitor = FALSE,
                lower = 0.01, upper = 0.99
  )
  
  c0 <- attr(fit, 'solution')[1,1]; names(c0) <- "v=1"
  val <- pfp.pfn.glm(corte = c0,
                    out = list(y=y, prob=prob))
  
  if(1 - val$esp > 1 - esp.target)
    warning(paste('PFN <= ',esp.target,' (target) não encontrado!\n Escolhendo PFP com PFN minimo!', sep=''))
  
  return(list(pfp=val$pfp,
              pfn=val$pfn,
              acc=val$acc,
              esp=val$esp,
              sen=val$sen,
              c0=c0,
              attain.esp.target=(1 - val$esp <= 1 - esp.target),
              esp.target=esp.target))
}

############
# Previsão #
############

#prever <-  m-n
#Xprev <-  X1[n:m,]
prev.glm <- function(fit, newY=NULL, newX, corte, li=0.000000001, ls=0.999999999){#Only new covariates!!
  prob.glm.p  <- predict(object = fit, newx = newX, type=c("response"))
  n.desc <- dim(newX)[1]
  
  classif.glm <- as.numeric(prob.glm.p > corte)
  
  if(!is.null(newY)){
    a <- length(which(newY==0 & classif.glm==0))
    b <- length(which(newY==0 & classif.glm==1))
    c <- length(which(newY==1 & classif.glm==0))
    d <- length(which(newY==1 & classif.glm==1))
    acc <- (a+d)/(a+b+c+d)
    esp <- ifelse(a+b==0, (a+c)/(a+b+c+d), a/(a+b))
    sen <- ifelse(c+d==0, (b+d)/(a+b+c+d), d/(c+d))
    pfp <- ifelse(b+d==0, (a+b)/(a+b+c+d), b/(b+d))
    pfn <- ifelse(a+c==0, (c+d)/(a+b+c+d), c/(a+c))
    return(list(prev=classif.glm, pfp=pfp, pfn=pfn, esp=esp, sen=sen, acc = acc))
  }
  else{
    return(list(prev=classif.glm))
  }
}

prev.MC <- function(beta, rho, newY, newX, y, X, corte0, corte1){
  y.f  <-  c(y, newY)
  X.f  <-  rbind(X, newX)
  prob.f  <-  predi(beta.vet=beta, rho = rho, covar=X.f, resp=y.f)

  s <- length(rho)
  N <- length(y.f); R <- N%/%s
  
  N0 <- length(y)+1
  
  classif.MC <- NULL
  for(r in 1:R){
    for(v in 1:s) {
      t <- (r-1)*s+v
      if(t>1){
        if(y.f[t-1] == 0 & prob.f[t] <= corte0[v]) classif.MC[t] <- 0
        if(y.f[t-1] == 0 & prob.f[t] > corte0[v]) classif.MC[t] <- 1
        if(y.f[t-1] == 1 & prob.f[t] <= corte1[v]) classif.MC[t] <- 0
        if(y.f[t-1] == 1 & prob.f[t] > corte1[v]) classif.MC[t] <- 1
      }
    }
  }
  a <- length(which(y.f[N0:N]==0 & classif.MC[N0:N]==0))
  b <- length(which(y.f[N0:N]==0 & classif.MC[N0:N]==1))
  c <- length(which(y.f[N0:N]==1 & classif.MC[N0:N]==0))
  d <- length(which(y.f[N0:N]==1 & classif.MC[N0:N]==1))
  
  acc.MC <- (a+d)/(a+b+c+d)
  esp.MC <- ifelse(a+b==0, (a+c)/(a+b+c+d), a/(a+b))
  sen.MC <- ifelse(c+d==0, (b+d)/(a+b+c+d), d/(c+d))
  pfp.MC <- ifelse(b+d==0, (a+b)/(a+b+c+d), b/(b+d))
  pfn.MC <- ifelse(a+c==0, (c+d)/(a+b+c+d), c/(a+c))
  
  return(list(prev=classif.MC[N0:N], pfp=pfp.MC, pfn=pfn.MC, esp = esp.MC, sen = sen.MC, acc = acc.MC))
}

