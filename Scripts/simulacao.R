#Note ANDREY
folder  <-  "/home/andrey/Projetos/PerLog/"
# PC TRABALHO
#folder <- 'C:/Users/Andrey/Documents/Andrey/TCC/'
# PC Casa Alessandro
# folder <- '~/Dropbox/alessandro/UFES/aulas/2019_1/Orientacao_IC_Andrey/'
source(paste(folder, 'Scripts/funcoes.R', sep = ''))
list_mean <- function(list_, name){
  out <- list_ %>% map(~{.x[name]}) %>% unlist %>% mean
  
  return(out)
}

arg <- commandArgs(trailingOnly = T)
# arg <- c(200, round(sin(2*pi*(1:7)/15)/2, 2), 50)
arg <- c(280, 1, .4, .5, .2, .4, .3, .5, .45, 10, 'Z')

# simulacao 
n <- as.numeric(arg[1]); k <- 2; rho <- as.numeric(arg[-c(1, 2,length(arg)-0:1)]);
beta <- c(rep(-1, k-1),  as.numeric(arg[2]));n.desc <- 7*40; m <- n+n.desc; esp.target <- .9
# X1 <- cbind(1, matrix( runif(n = m*(k-1), min = -1.4, max = -0.6), m, (k-1)))
# X1 <- cbind(1, (.5*sin(2*pi*(1:m)/28)+runif(m,-1.25,-.75)))
# X1 <- cbind(1, matrix( seq(from = -2, to = 1, length.out = m*(k-1)), m, (k-1)))
# colnames(X1) <- c('Intercept', kronecker("X_", 1:(k-1), paste, sep=''))
# X <- X1[1:n,]; Xprev <-  X1[(n+1):m,]
REP <- as.numeric(arg[length(arg)-1])

# plot.ts((X%*%beta)[,1])
# plot.ts(exp(X%*%beta)/(1 + exp(X%*%beta))[,1])

est_glm <- yhat_glm <- std_glm <- est_MCP <- std_MCP <- est_MC <- std_MC <- fit_glmvet <- corte <- NULL
esp.insample.MCP <- esp.insample.MC <- esp.insample.glm <- sen.insample.MCP <- sen.insample.MC <- sen.insample.glm <- acc.insample.MCP  <- acc.insample.MC <- acc.insample.glm <- NULL
esp.glm <- sen.glm <- acc.glm <- esp.MCP <- sen.MCP <- acc.MCP <- esp.MC <- sen.MC <- acc.MC <- NULL
corte0.MC <- corte1.MC <- corte0.MCP <- corte1.MCP <- list()
fit_glmvet <- fit_MCPvet  <- fit_MCvet <- list()
est_glm <- std_glm <- matrix(nrow = REP, ncol = k)
est_MCP <- std_MCP <- matrix(nrow = REP, ncol = (k+length(rho)))
est_MC <- std_MC <- matrix(nrow = REP, ncol = (k+1))
attain.esp.target.glm <- attain.esp.target.MC <- attain.esp.target.MCP <- NULL

# Leitura necessária caso haja problemas com a sessão
{
  if(!dir.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/'))) dir.create(path = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/'), recursive = T)
  if(!dir.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Metricas/', arg[length(arg)], '/'))) dir.create(path = paste0('/home/andrey/Projetos/PerLog/Dados/Metricas/', arg[length(arg)], '/'), recursive = T)
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/resultados.rds'))) load(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/resultados.rds'))
}
r <- 1 + length(esp.insample.glm)

ct.fim.vet <- prob <- fit_MCvet <- fit_MCPvet <- y.desc <- results <- list()
set.seed(325 + length(esp.insample.glm))
while(r <= REP){
  cat('Repl.', r)
  # generate data
  X1 <- cbind(1, matrix( runif(n = m*(k-1), min = 0, max = 2), m, (k-1)))
  colnames(X1) <- c('Intercept', kronecker("X_", 1:(k-1), paste, sep=''))
  X <- X1[1:n,]; Xprev <-  X1[(n+1):m,]
  ytotal <- gen(X = X1, beta = beta, rho = rho)
  y.desc[[r]] <- ytotal[(n+1):m]
  y <- ytotal[1:n]
  # w <- mean(y); w <- ifelse(y == 1, 1/w, 1/(1-w))
  #lassim  <-  cv.glmnet(x=matr, y=O3, family = 'binomial', alpha=1, standardize = TRUE)
  
  # FIT GLM
  
  fit_glm1 <-  glm(y~X-1, family = 'binomial')
  fit_glmvet[[r]] <-  (fit_glm1)
  fit_glm <-  summary(fit_glm1)$coefficients
  #lambda <-  rbind(lambda, lassim$lambda)
  yhat_glm <- rbind(yhat_glm,fitted(fit_glm1))
  
  # FIT MC
  
  fit_MC <- tryCatch(logistic_MC(resp = y, covar = X, par0 = c(beta, mean(rho)), trace = 0),
                     error = function(e){NULL})
  fit_MCvet[[r]] <- fit_MC
  
  #FIT MCP
  
  fit_MCP <- tryCatch(logistic_MC(resp = y, covar = X, par0 = c(beta, rho), trace = 0),
                      error = function(e){NULL})
  fit_MCPvet[[r]] <- fit_MCP
  if(!is.null(fit_MCP) & !is.null(fit_MC)){
    est_glm[r,] <- fit_glm[,1]
    std_glm[r,] <- fit_glm[,2]
    est_MC[r,] <- fit_MC[,1]
    std_MC[r,] <- fit_MC[,2]
    est_MCP[r,] <- fit_MCP[,1]
    std_MCP[r,] <- fit_MCP[,2]
    
    # prob[[r]] <- predi(beta.vet=fit_MCPvet[[r]][1:k,1],
    # rho = fit_MCPvet[[r]][(1:length(rho)+k),1], covar=X, resp=y)
    
    prob.glm  <-   predict(object = fit_glmvet[[r]], newx = X, type=c("response"))
    prob.MC <- predi(beta.vet=fit_MCvet[[r]][1:k,1],
                     rho = fit_MCvet[[r]][(1+k),1], covar=X, resp=y)
    prob.MCP  <-  predi(beta.vet=fit_MCPvet[[r]][1:k,1],
                        rho = fit_MCPvet[[r]][(1:length(rho)+k),1], covar=X, resp=y)
    
    # GLM
    in.sample.glm <- min.pfp.glm(y = y, prob = prob.glm, sen.target = esp.target)
    esp.insample.glm[r] <- in.sample.glm$esp
    sen.insample.glm[r] <- in.sample.glm$sen
    acc.insample.glm[r] <- in.sample.glm$acc
    corte[r] <- in.sample.glm$c
    attain.esp.target.glm[r] <- in.sample.glm$attain.esp.target
    pr.glm <- prev.glm(fit = fit_glmvet[[r]], newY = y.desc[[r]], newX = Xprev,
                       corte = corte[r])
    esp.glm[r] <- pr.glm$esp; sen.glm[r] <- pr.glm$sen; acc.glm[r] <- pr.glm$acc
    
    # MC
    in.sample.MC <- min.pfp.MC(y = y, prob = prob.MC, esp.target = esp.target, s = 1)
    esp.insample.MC[r] <- in.sample.MC$esp
    sen.insample.MC[r] <- in.sample.MC$sen
    acc.insample.MC[r] <- in.sample.MC$acc
    corte0.MC[[r]] <- in.sample.MC$c0
    corte1.MC[[r]] <- in.sample.MC$c1
    attain.esp.target.MC[r] <- in.sample.MC$attain.esp.target
    pr.MC <- prev.MC(beta = fit_MCvet[[r]][1:k,1],
                     rho  = fit_MCvet[[r]][(k+1),1],
                     newY = y.desc[[r]], newX = Xprev,
                     y= y, X = X, corte0 = corte0.MC[[r]], corte1 = corte1.MC[[r]])
    esp.MC[r] <- pr.MC$esp; sen.MC[r] <- pr.MC$sen; acc.MC[r] <- pr.MC$acc
    
    # MCP
    in.sample.MCP <- min.pfp.MC(y = y, prob = prob.MCP, esp.target = esp.target, s = length(rho))
    esp.insample.MCP[r] <- in.sample.MCP$esp
    sen.insample.MCP[r] <- in.sample.MCP$sen
    acc.insample.MCP[r] <- in.sample.MCP$acc
    corte0.MCP[[r]] <- in.sample.MCP$c0
    corte1.MCP[[r]] <- in.sample.MCP$c1
    attain.esp.target.MCP[r] <- in.sample.MCP$attain.esp.target
    pr.MCP <- prev.MC(beta = fit_MCPvet[[r]][1:k,1],
                      rho  = fit_MCPvet[[r]][(k+length(rho)),1],
                      newY = y.desc[[r]], newX = Xprev,
                      y= y, X = X, corte0 = corte0.MCP[[r]], corte1 = corte1.MCP[[r]])
    esp.MCP[r] <- pr.MCP$esp; sen.MCP[r] <- pr.MCP$sen; acc.MCP[r] <- pr.MCP$acc
    #Escrita das variáveis
    {
      results[[r]] <- list(beta=beta, rho=rho, 
                           X1=X1, X=X, Xprev=Xprev, 
                           ytotal=ytotal, y.desc=y.desc[[r]], y=y, 
                           fit_glm=fit_glm1, coef_glm=fit_glm, prob.glm=prob.glm, 
                           esp.insample.glm=in.sample.glm$esp, sen.insample.glm=in.sample.glm$sen, acc.insample.glm=in.sample.glm$acc, 
                           corte0.glm=in.sample.glm$c0, corte1.glm=in.sample.glm$c1, attain.esp.target.glm=in.sample.glm$attain.esp.target, 
                           pr.glm=pr.glm, esp.glm=pr.glm$esp, sen.glm=pr.glm$sen, acc.glm=pr.glm$acc,
                           fit_MC=fit_MC,  prob.MC=prob.MC, prob.MC=prob.MC,
                           esp.insample.MC=in.sample.MC$esp, sen.insample.MC=in.sample.MC$sen, acc.insample.MC=in.sample.MC$acc, 
                           corte0.MC=in.sample.MC$c0, corte1.MC=in.sample.MC$c1, attain.esp.target.MC=in.sample.MC$attain.esp.target, 
                           pr.MC=pr.MC, esp.MC=pr.MC$esp, sen.MC=pr.MC$sen, acc.MC=pr.MC$acc,
                           fit_MCP=fit_MCP, prob.MCP=prob.MCP, prob.MCP=prob.MCP,
                           esp.insample.MCP=in.sample.MCP$esp, sen.insample.MCP=in.sample.MCP$sen, acc.insample.MCP=in.sample.MCP$acc, 
                           corte0.MCP=in.sample.MCP$c0, corte1.MCP=in.sample.MCP$c1, attain.esp.target.MCP=in.sample.MCP$attain.esp.target,
                           pr.MCP=pr.MCP, esp.MCP=pr.MCP$esp, sen.MCP=pr.MCP$sen, acc.MCP=pr.MCP$acc)
      
      save(list = "results", file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/resultados.rds', sep = ''))
      
    }
    r  <-  r+1
    cat('\n')
  }else{
    if(is.null(fit_MC)) cat('Deu ruim - Katiuski!\n')
    if(is.null(fit_MCP)) cat('Deu ruim - Periodic!\n')
  }
}

pf_insample_mat <- rbind(c(list_mean(results, "esp.insample.glm"), list_mean(results, "esp.insample.MC"), list_mean(results, "esp.insample.MCP")),
                         c(list_mean(results, "sen.insample.glm"), list_mean(results, "sen.insample.MC"), list_mean(results, "sen.insample.MCP")),
                         c(list_mean(results, "acc.insample.glm"), list_mean(results, "acc.insample.MC"), list_mean(results, "acc.insample.MCP")))
pf_out_mat <- rbind(c(list_mean(results, "esp.glm"), list_mean(results, "esp.MC"), list_mean(results, "esp.MCP")),
                    c(list_mean(results, "sen.glm"), list_mean(results, "sen.MC"), list_mean(results, "sen.MCP")),
                    c(list_mean(results, "acc.glm"), list_mean(results, "acc.MC"), list_mean(results, "acc.MCP")))
pf_mat <- rbind(pf_insample_mat, pf_out_mat)
rownames(pf_mat) <- c('esp - insample','sen - insample','ACC - insample', 'esp - out', 'sen - out', 'ACC - out')
colnames(pf_mat) <- c('GLM', 'MC','MCP')


est_glm <- results %>% map(~{.x$coef_glm[,1]}) %>% bind_rows()
est_MC <- results %>% map(~{.x$fit_MC[,1]}) %>% bind_rows()
est_MCP <- results %>% map(~{.x$fit_MCP[,1]}) %>% bind_rows()
std_glm <- results %>% map(~{.x$coef_glm[,2]}) %>% bind_rows()
std_MC <- results %>% map(~{.x$fit_MC[,2]}) %>% bind_rows()
std_MCP <- results %>% map(~{.x$fit_MCP[,2]}) %>% bind_rows()

Bias <- cbind(c(colMeans(est_glm)-beta, rep(NA, length(rho))),
              c(colMeans(est_MC)-c(beta, mean(rho)), rep(NA, length(rho)- 1)),
              (colMeans(est_MCP)-c(beta, rho)))
EP <- cbind(c(sqrt(colMeans(est_glm^2)-colMeans(est_glm)^2), rep(NA, length(rho))),
            c(sqrt(colMeans(est_MC^2)-colMeans(est_MC)^2), rep(NA, length(rho)-1)),
            sqrt(colMeans(est_MCP^2)-colMeans(est_MCP)^2))
RMSE <- sqrt(Bias^2+EP^2)
mean.std <- cbind(c(colMeans(std_glm), rep(NA, length(rho))),
                  c(colMeans(std_MC,na.rm = T), rep(NA, length(rho)- 1)),
                  colMeans(std_MCP))
rownames(Bias) <- rownames(EP) <- rownames(RMSE) <- rownames(mean.std) <- names(colMeans(est_MCP))
mat_Bias_RMSE <- cbind(Bias, RMSE); colnames(mat_Bias_RMSE) <- c('Bias - GLM', 'Bias - MC', 'Bias - MCP', 'RMSE - GLM', 'RMSE - MC', 'RMSE - MCP')

mat.EP <- cbind(EP[,1], mean.std[,1], EP[,2], mean.std[,2], EP[,3], mean.std[,3])
colnames(mat.EP) <- c('Emp. EP - GLM', 'Hes. EP - GLM', 'Emp. EP - MC', 'Hes. EP - MC', 'Emp. EP - MCP', 'Hes. EP - MCP')

### RESULTADOS PARA GUARDAR
paste0('n = ', n, ', rho = (', paste(rho, collapse = ', '), ')')
PFN_IN_GLM <- mean(!attain.esp.target.glm)
PFN_IN_MC <- mean(!attain.esp.target.MC)
PFN_IN_MCP <- mean(!attain.esp.target.MCP)
PFN_IN_GLM; PFN_IN_MC; PFN_IN_MCP
pf_mat
mat_Bias_RMSE
mat.EP

save(list = c('PFN_IN_GLM', 'PFN_IN_MC', 'PFN_IN_MCP', 'pf_mat', 'mat_Bias_RMSE', 'mat.EP'), file = paste0('/home/andrey/Projetos/PerLog/Dados/Metricas/', arg[length(arg)], '/resultado.rds', sep = ''))