resp <- y; covar <- X; par <- c(beta, rho)
li <- 0.000000001; ls <- 0.999999999
n <- length(resp); k <- dim(covar)[2]
beta.vet  <-  par[1:k]
# rho  <-  exp(par[k+1])
rho  <-  par[(k+1):length(par)]
s <- length(rho)
per <- ifelse(2:n%%s == 0, s, 2:n%%s)

p <- theta <- NULL

for(t in 1:n){
  theta[t] <- exp(sum(covar[t,]*beta.vet))/(1+exp(sum(covar[t,]*beta.vet)))
  # theta[t]  <-  1/(1+exp(-sum(covar[t,]*beta.vet)))
}
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
dpdrho_df <- permatrix(vector = dpdrho, period = per)
dldrho  <-  dldp*dpdrho
dldrho_m <- dldp*dpdrho_df

# gradient
grad  <-  c(colSums(dldb), colSums(dldrho_m))

#hessian
d2ldp2  <-  -resp/(p^2)-(1-resp)/((1-p)^2)
d2ldrho2  <-  matrix(0, s, s)
for (t in 1:ceiling(n/s)) {
  d2ldrho2  <-  d2ldrho2 + d2ldp2[t]*dpdrho_df[t,]%*%t(dpdrho_df[t,])
}

d2ldb2  <-  matrix(0, k, k)
for (t in 1:n) {
  d2ldb2  <-  d2ldb2 + d2ldp2[t]*dpdb[t,]%*%t(dpdb[t,])
}

d2pdrhodvt0  <-  -.5*(c1^(-.5))*(1-2*theta[2:n])*(sqrt(c2)-resp[2:n-1]*(sqrt(c2)+sqrt(c3))); d2pdrhodvt0  <-  c(0,d2pdrhodvt0);
d2pdrhodvt0 <- permatrix(vector = d2pdrhodvt0, period = per)
d2pdrhodvt1  <-  .5*sqrt(c1)*(c4-resp[2:n-1]*(c4-c5)); d2pdrhodvt1  <-  c(0,d2pdrhodvt1);
d2pdrhodvt1 <- permatrix(vector = d2pdrhodvt1, period = per)
d2pdbdrho  <-  t(d2pdrhodvt0)%*%dvt0db + t(d2pdrhodvt1)%*%dvt1db

d2ldb2  <-  matrix(0, k, k)
for (t in 1:n) {
  d2ldb2  <-  d2ldb2 + d2ldp2[t]*dpdb[t,]%*%t(dpdb[t,])
  d2pdbdrho  <-  t(d2pdrhodvt0)[,t]%*%dvt0db[] + t(d2pdrhodvt1)%*%dvt1db
  d2ldbdrho  <-  d2ldbdrho + as.vector(d2ldp2)[t] * t(dpdb)[,t] %*% dpdrho_df[t,] + as.vector(dldp)[t]*d2pdbdrho
}



d2ldbdrho <- 
  data.frame(d2ldbdrho, per = c(1, per)) %>% 
  group_by(per) %>% 
  mutate(ID = row_number()) %>% 
  summarise_all(sum) %>% select(-ID, - per)

hess  <-  cbind(rbind(d2ldb2, d2ldbdrho), 
                rbind( t(d2ldbdrho), d2ldrho2))

logver <-  -sum(resp*log(p/(1-p))+log(1-p))
attr(logver, 'gradient')  <-  -grad/n
attr(logver, 'hessian')  <-  -hess
return(logver)