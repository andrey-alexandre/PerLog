# ##########
# # Andrey #
# ##########
# otim  <-  optim(par = c(beta, rho), fn = loglik, gr = gr.loglik, out=list(covar=X1, resp=ytotal, li=0.000000001, ls=0.999999999), method = "L-BFGS-B",
#                 lower = c(rep(-Inf, dim(X1)[2]), rep(-0.999999, length( c(beta, rho)) - dim(X1)[2])),
#                 upper = c(rep(Inf, dim(X1)[2]), rep(0.999999, length( c(beta, rho)) - dim(X1)[2])), hessian = T,
#                 control = list(trace=2))
# 
# otim$hessian 
# hess.loglik(par = c(beta, rho), out = list(covar=X1, resp=ytotal, li=0.000000001, ls=0.999999999))
# 
# ############
# # Katiuski #
# ############
# otim.kat  <-  optim(par = c(beta, mean(rho)), fn = loglik.kat, gr = gr.loglik.kat, out=list(covar=X1, resp=ytotal, li=0.000000001, ls=0.999999999), method = "L-BFGS-B",
#                 lower = c(rep(-Inf, dim(X1)[2]), rep(-0.999999, length( c(beta, mean(rho))) - dim(X1)[2])),
#                 upper = c(rep(Inf, dim(X1)[2]), rep(0.999999, length( c(beta, mean(rho))) - dim(X1)[2])), 
#                 hessian = T, control = list(trace=2))

otim.kat$hessian 
hess.loglik(par = c(beta, mean(rho)), out = list(covar=X1, resp=ytotal, li=0.000000001, ls=0.999999999)) #%>% solve %>% diag %>% sqrt
hess.loglik.kat(par = c(beta, mean(rho)), out = list(covar=X1, resp=ytotal, li=0.000000001, ls=0.999999999))  #%>% solve %>% diag %>% sqrt
numDeriv::hessian(func = loglik.kat, x = c(beta, rho), out = list(covar=X1, resp=ytotal, li=0.000000001, ls=0.999999999)) #%>% solve %>% diag %>% sqrt

gr.loglik(par = c(beta, mean(rho)), out = list(covar=X1, resp=ytotal, li=0.000000001, ls=0.999999999)) 
gr.loglik.kat(par = c(beta, mean(rho)), out = list(covar=X1, resp=ytotal, li=0.000000001, ls=0.999999999)) 
numDeriv::grad(func = loglik.kat, x = c(beta, rho), out = list(covar=X1, resp=ytotal, li=0.000000001, ls=0.999999999))

