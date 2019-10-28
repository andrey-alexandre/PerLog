resp <- y; covar <- X
li <- 0.000000001; ls <- 0.999999999
n <- length(resp); k <- dim(covar)[2]
par <- c(beta, rho)
beta.vet  <-  par[1:k]
# rho  <-  exp(par[k+1])
rho  <-  par[(k+1):length(par)]
s <- length(rho)
per <- ifelse(2:n%%s == 0, s, 2:n%%s)

p <- NULL

theta <- exp(covar%*%beta.vet)/(1+exp(covar%*%beta.vet)); theta <- as.vector(theta)
theta[theta<=li] <- li;theta[theta>=ls] <- ls

p[1]  <-  theta[1]

c1  <-  theta[2:n]*(1-theta[2:n])
c1aux  <-  theta[2:n-1]*(1-theta[2:n-1]); c1aux <- c(0, c1aux)
c2  <-  theta[2:n-1]/(1-theta[2:n-1])
c3  <-  (1-theta[2:n-1])/theta[2:n-1]
c4  <-  (theta[2:n-1]^(-.5))*((1-theta[2:n-1])^(-1.5))
c5  <-  (theta[2:n-1]^(-1.5))*((1-theta[2:n-1])^(-.5))

# c1[c1==-Inf] <- -999999999999999999999999;c1[c1==Inf] <- 999999999999999999999999
# c2[c2==-Inf] <- -999999999999999999999999;c2[c2==Inf] <- 999999999999999999999999
# c3[c3==-Inf] <- -999999999999999999999999;c3[c3==Inf] <- 999999999999999999999999
# c4[c4==-Inf] <- -999999999999999999999999;c4[c4==Inf] <- 999999999999999999999999
# c5[c5==-Inf] <- -999999999999999999999999;c5[c5==Inf] <- 999999999999999999999999

p[2:n]  <-  theta[2:n] - rho[per]*sqrt(c1)*(sqrt(c2)-resp[2:n-1]*(sqrt(c2)+sqrt(c3)))

p[p<=li] <- li;p[p>=ls] <- ls

dpdvt0  <-  1-.5*rho[per]*(c1^(-.5))*(1-2*theta[2:n])*(sqrt(c2)-resp[2:n-1]*(sqrt(c2)+sqrt(c3)))
dpdvt1  <-  .5*rho[per]*sqrt(c1)*(c4-resp[2:n-1]*(c4-c5))
dpdvt0  <-  c(1,dpdvt0)
dpdvt1  <-  c(0,dpdvt1)

dldp  <-  (resp-p)/(p*(1-p))
dvt0db  <-  theta*(1-theta)*covar
dvt1db  <-  c1aux*covar

dpdb  <-  dpdvt0*dvt0db + dpdvt1*dvt1db
dldb  <-  as.vector(dldp)*dpdb

dpdrho  <-  -sqrt(c1)*(sqrt(c2)-resp[2:n-1]*(sqrt(c2)+sqrt(c3))); dpdrho  <-  c(0, dpdrho)
dpdrho <- permatrix(vector = dpdrho, period_vec = per, period = s)
dldrho  <-  dldp*dpdrho

# gradient
grad  <-  c(colSums(dldb), colSums(dldrho))

#hessian
d2ldp2  <-  -resp/(p^2)-(1-resp)/((1-p)^2)
d2ldrho2  <-  array(0, dim = c(s, s, n))
for (t in 1:n) {
  d2ldrho2[ , ,t]  <-  d2ldp2[t]*dpdrho[t,]%*%t(dpdrho[t,])
}

d2ldb2  <- array(0, dim = c(k, k, n))
for (t in 1:n) {
  d2ldb2[ , , t]  <- d2ldp2[t]*dpdb[t,]%*%t(dpdb[t,])
}

d2pdrhodvt0  <-  -.5*(c1^(-.5))*(1-2*theta[2:n])*(sqrt(c2)-resp[2:n-1]*(sqrt(c2)+sqrt(c3))); d2pdrhodvt0  <-  c(0,d2pdrhodvt0);
d2pdrhodvt0 <- permatrix(vector = d2pdrhodvt0, period_vec = per, period = s)
d2pdrhodvt1  <-  .5*sqrt(c1)*(c4-resp[2:n-1]*(c4-c5)); d2pdrhodvt1  <-  c(0,d2pdrhodvt1);
d2pdrhodvt1 <- permatrix(vector = d2pdrhodvt1, period_vec = per, period = s)

d2pdbdrho  <- array(0, dim = c(k, s, n))
for (t in 1:n) {
  d2pdbdrho[ , , t]  <- t(dvt0db)[,t, drop = F]%*%d2pdrhodvt0[t, , drop = F] + t(dvt1db)[,t, drop = F]%*% d2pdrhodvt1[t, , drop = F]
}

d2ldbdrho  <-  matrix(0, k, s)
for (t in 1:n) {
  d2ldbdrho  <-  d2ldbdrho + as.vector(d2ldp2)[t] * t(dpdb)[, t, drop = F] %*%  dpdrho[t, , drop = F] + as.vector(dldp)[t]*d2pdbdrho[, , t]
}

hess  <-  cbind(rbind(rowSums(x = d2ldb2, dims = 2), t(d2ldbdrho)), 
                rbind( d2ldbdrho, rowSums(x = d2ldrho2, dims = 2)))

diag( solve( -hess ) )




out <- list(covar=covar, resp=resp, li=li, ls=ls)
otim  <-  optim(par = par, fn = loglik, gr = gr.loglik, out=out, method = "L-BFGS-B",
                lower = c(rep(-Inf, dim(covar)[2]), -0.999999),
                upper = c(rep(Inf, dim(covar)[2]), 0.999999), hessian = T,
                control = list(trace=0))

diag(solve(hess.loglik(par = otim$par, out = out)))
diag( solve( -hess ) )
diag( solve( otim$hessian ) )
# seotim <- sqrt(diag(solve(otim$hessian)))
seotim  <-  sqrt(diag(solve(hess.loglik(par = otim$par, out = out))))
mat.res <- cbind(otim$par, seotim, 2*(1-pnorm(abs(otim$par/seotim))))
colnames(mat.res) <- c('Estimate','Std.Error','p-value')
# row.names(mat.res)[dim(covar)[2]+1] <- 'lambda'
row.names(mat.res) <- c(colnames(covar), paste0('rho', seq(length(par0)-dim(covar)[2])))
return(mat.res)