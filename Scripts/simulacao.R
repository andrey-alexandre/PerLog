#Note ANDREY
folder  <-  "/home/andrey/Projetos/PerLog/"
# PC TRABALHO
#folder <- 'C:/Users/Andrey/Documents/Andrey/TCC/'
# PC Casa Alessandro
#folder <- '/home/alessandro/Dropbox/alessandro/2018_2/Orientacao_Monografia_Katiuski/TCCKatiuski/Monografia/codigo/'
source(paste(folder, 'Scripts/funcoes.R', sep = ''))
# source(paste(folder, 'Scripts/FuncKat.R', sep = ''))
arg <- commandArgs(trailingOnly = T)
# arg <- c(200, round(sin(2*pi*(1:7)/15)/2, 2), 50)
# arg <- c(280, 0, .4, .5, .2, .4, .3, .5, .45, 1000, 'Z')

# simulacao 
n <- as.numeric(arg[1]); k <- 2; rho <- as.numeric(arg[-c(1, 2,length(arg)-0:1)]);
beta <- c(rep(0.5, k-1),  as.numeric(arg[2]));n.desc <- 100;m <- n+n.desc; pfn.target <- .1
X1 <- cbind(1, matrix( runif(n = m*(k-1), min = -1, max = 1), m, (k-1)))
# X1 <- cbind(1, sin(2*pi*(1:m)/21))
colnames(X1) <- c('Intercept', kronecker("X_", 1:(k-1), paste, sep=''))
X <- X1[1:n,]; Xprev <-  X1[(n+1):m,]
REP <- as.numeric(arg[length(arg)-1])

plot.ts(exp(X%*%beta)/(1 + exp(X%*%beta))[,1])

est_glm <- yhat_glm <- std_glm <- est_MCP <- std_MCP <- est_MC <- std_MC <- fit_glmvet <- corte0 <- corte1 <- corte <- NULL
pfp.insample.MCP <- pfp.insample.MC <- pfp.insample.glm <- pfn.insample.MCP <- pfn.insample.MC <- pfn.insample.glm <- acc.insample.MCP  <- acc.insample.MC <- acc.insample.glm <- NULL
pfp.glm <- pfn.glm <- acc.glm <- pfp.MCP <- pfn.MCP <- acc.MCP <- pfp.MC <- pfn.MC <- acc.MC <- NULL
fit_glmvet <- fit_MCPvet  <- fit_MCvet <- list()
est_glm <- std_glm <- matrix(nrow = REP, ncol = k)
est_MCP <- std_MCP <- matrix(nrow = REP, ncol = (k+length(rho)))
est_MC <- std_MC <- matrix(nrow = REP, ncol = (k+1))

# Leitura necessária caso haja problemas com a sessão
{
  if(!dir.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/'))) dir.create(path = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/'), recursive = T)
  if(!dir.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Metricas/', arg[length(arg)], '/'))) dir.create(path = paste0('/home/andrey/Projetos/PerLog/Dados/Metricas/', arg[length(arg)], '/'), recursive = T)
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfp.insample.glm.csv'))) pfp.insample.glm <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfp.insample.glm.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfn.insample.glm.csv'))) pfn.insample.glm <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfn.insample.glm.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_acc.insample.glm.csv'))) acc.insample.glm <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_acc.insample.glm.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfp.glm.csv'))) pfp.glm <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfp.glm.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfn.glm.csv'))) pfn.glm <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfn.glm.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_acc.glm.csv'))) acc.glm <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_acc.glm.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_est_glm.csv'))) est_glm <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_est_glm.csv'), stringsAsFactors = F, header = F)
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_std_glm.csv'))) std_glm <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_std_glm.csv'), stringsAsFactors = F, header = F)
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfp.insample.MC.csv'))) pfp.insample.MC <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfp.insample.MC.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfn.insample.MC.csv'))) pfn.insample.MC <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfn.insample.MC.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_acc.insample.MC.csv'))) acc.insample.MC <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_acc.insample.MC.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfp.MC.csv'))) pfp.MC <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfp.MC.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfn.MC.csv'))) pfn.MC <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfn.MC.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_acc.MC.csv'))) acc.MC <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_acc.MC.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_est_MC.csv'))) est_MC <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_est_MC.csv'), stringsAsFactors = F, header = F)
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_std_MC.csv'))) std_MC <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_std_MC.csv'), stringsAsFactors = F, header = F)
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfp.insample.MCP.csv'))) pfp.insample.MCP <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfp.insample.MCP.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfn.insample.MCP.csv'))) pfn.insample.MCP <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfn.insample.MCP.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_acc.insample.MCP.csv'))) acc.insample.MCP <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_acc.insample.MCP.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfp.MCP.csv'))) pfp.MCP <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfp.MCP.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfn.MCP.csv'))) pfn.MCP <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_pfn.MCP.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_acc.MCP.csv'))) acc.MCP <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_acc.MCP.csv'), stringsAsFactors = F, header = F) %>% unlist()
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_est_MCP.csv'))) est_MCP <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_est_MCP.csv'), stringsAsFactors = F, header = F)
  if(file.exists(paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_std_MCP.csv'))) std_MCP <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/', arg[length(arg)], '/parcial_std_MCP.csv'), stringsAsFactors = F, header = F)
}
r <- 1 + length(pfp.insample.glm)

ct.fim.vet <- prob <- fit_MCvet <- fit_MCPvet <- y.desc <- list()
set.seed(444 + length(pfp.insample.glm))
while(r <= REP){
  cat('Repl.', r)
  ytotal <- gen(X = X1, beta = beta, rho = rho)
  y.desc[[r]] <- ytotal[(n+1):m]
  y <- ytotal[1:n]
  #lassim  <-  cv.glmnet(x=matr, y=O3, family = 'binomial', alpha=1, standardize = TRUE)
  fit_glm1 <-  glm(y~X-1, family = 'binomial')
  fit_glmvet[[r]] <-  (fit_glm1)
  fit_glm <-  summary(fit_glm1)$coefficients
  #lambda <-  rbind(lambda, lassim$lambda)
  yhat_glm <- rbind(yhat_glm,fitted(fit_glm1))
  fit_MC <- tryCatch(logistic_MC(resp = y, covar = X, par0 = c(beta, mean(rho)), trace = 0),
                     error = function(e){NULL})
  fit_MCvet[[r]] <- fit_MC
  # fit_MCP <- logistic_MCP(resp = y, covar = X, par0 = c(beta, rho), trace = 0)
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
    
    prob[[r]] <- predi(beta.vet=fit_MCPvet[[r]][1:k,1],
                       rho = fit_MCPvet[[r]][(1:length(rho)+k),1], covar=X, resp=y)
    
    prob.glm  <-   predict(object = fit_glmvet[[r]], newx = X, type=c("response"))
    prob.MC <- predi(beta.vet=fit_MCvet[[r]][1:k,1],
                         rho = fit_MCPvet[[r]][(1+k),1], covar=X, resp=y)
    prob.MCP  <-  prob[[r]]
    in.sample.glm <- min.pfp.glm(y = y, prob = prob.glm, pfn.target = pfn.target)
    pfp.insample.glm[r] <- in.sample.glm$pfp.min
    pfn.insample.glm[r] <- in.sample.glm$pfn.min
    acc.insample.glm[r] <- in.sample.glm$acc.min
    corte[r] <- in.sample.glm$cortemin
    
    pr.glm <- prev.glm(fit = fit_glmvet[[r]], newY = y.desc[[r]], newX = Xprev,
                       corte = corte[r])
    pfp.glm[r] <- pr.glm$pfp; pfn.glm[r] <- pr.glm$pfn; acc.glm[r] <- pr.glm$acc
    
    in.sample.MC <- min.pfp.MC(y = y, prob = prob.MC, pfn.target = pfn.target)
    pfp.insample.MC[r] <- in.sample.MC$pfp.min
    pfn.insample.MC[r] <- in.sample.MC$pfn.min
    acc.insample.MC[r] <- in.sample.MC$acc.min
    corte0[r] <- in.sample.MC$corte0min
    corte1[r] <- in.sample.MC$corte1min
    
    pr.MC <- prev.MC(beta = fit_MCvet[[r]][1:k,1],
                     rho  = fit_MCvet[[r]][(k+1),1],
                     newY = y.desc[[r]], newX = Xprev,
                     LastY= y[n], LastX = X[n,], corte0 = corte0[r], corte1 = corte1[r])
    pfp.MC[r] <- pr.MC$pfp; pfn.MC[r] <- pr.MC$pfn; acc.MC[r] <- pr.MC$acc
    
    in.sample.MCP <- min.pfp.MC(y = y, prob = prob.MCP, pfn.target = pfn.target)
    pfp.insample.MCP[r] <- in.sample.MCP$pfp.min
    pfn.insample.MCP[r] <- in.sample.MCP$pfn.min
    acc.insample.MCP[r] <- in.sample.MCP$acc.min
    corte0[r] <- in.sample.MCP$corte0min
    corte1[r] <- in.sample.MCP$corte1min
    
    pr.MCP <- prev.MC(beta = fit_MCPvet[[r]][1:k,1],
                      rho  = fit_MCPvet[[r]][(k+length(rho)),1],
                      newY = y.desc[[r]], newX = Xprev,
                      LastY= y[n], LastX = X[n,], corte0 = corte0[r], corte1 = corte1[r])
    pfp.MCP[r] <- pr.MCP$pfp; pfn.MCP[r] <- pr.MCP$pfn; acc.MCP[r] <- pr.MCP$acc
    #Escrita das variáveis 
    {
      write.table(x = pfp.insample.glm, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_pfp.insample.glm.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = pfn.insample.glm, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_pfn.insample.glm.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = acc.insample.glm, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_acc.insample.glm.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = pfp.glm, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_pfp.glm.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = pfn.glm, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_pfn.glm.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = acc.glm, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_acc.glm.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = est_glm, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_est_glm.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = std_glm, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_std_glm.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = pfp.insample.MC, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_pfp.insample.MC.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = pfn.insample.MC, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_pfn.insample.MC.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = acc.insample.MC, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_acc.insample.MC.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = pfp.MC, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_pfp.MC.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = pfn.MC, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_pfn.MC.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = acc.MC, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_acc.MC.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = est_MC, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_est_MC.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = std_MC, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_std_MC.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = pfp.insample.MCP, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_pfp.insample.MCP.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = pfn.insample.MCP, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_pfn.insample.MCP.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = acc.insample.MCP, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_acc.insample.MCP.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = pfp.MCP, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_pfp.MCP.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = pfn.MCP, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_pfn.MCP.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = acc.MCP, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_acc.MCP.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = est_MCP, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_est_MCP.csv'), row.names = F, col.names = F, sep = ',')
      write.table(x = std_MCP, file = paste0(folder,'Dados/Parciais/', arg[length(arg)], '/parcial_std_MCP.csv'), row.names = F, col.names = F, sep = ',')
    }
    r  <-  r+1
    cat('\n')
  }else{
    if(is.null(fit_MC)) cat('Deu ruim - Katiuski!\n')
    if(is.null(fit_MC)) cat('Deu ruim - Periodic!\n')
  }
}

pf_insample_mat <- rbind(c(mean(pfp.insample.glm), mean(pfp.insample.MC), mean(pfp.insample.MCP)),
                         c(mean(pfn.insample.glm), mean(pfn.insample.MC), mean(pfn.insample.MCP)),
                         c(mean(acc.insample.glm), mean(acc.insample.MC), mean(acc.insample.MCP)))
pf_out_mat <- rbind(c(mean(pfp.glm), mean(pfp.MC), mean(pfp.MCP)),
                    c(mean(pfn.glm), mean(pfn.MC), mean(pfn.MCP)),
                    c(mean(acc.glm), mean(acc.MC), mean(acc.MCP)))
pf_mat <- rbind(pf_insample_mat, pf_out_mat)
rownames(pf_mat) <- c('PFP - insample','PFN - insample','ACC - insample', 'PFP - out', 'PFN - out', 'ACC - out')
colnames(pf_mat) <- c('GLM', 'MC','MCP')

Bias <- cbind(c(colMeans(est_glm)-beta, rep(NA, length(rho))), c(colMeans(est_MC)-c(beta, rho[1]), rep(NA, length(rho)- 1)), (colMeans(est_MCP)-c(beta, rho)))
EP <- cbind(c(sqrt(colMeans(est_glm^2)-colMeans(est_glm)^2), rep(NA, length(rho))),
            c(sqrt(colMeans(est_MC^2)-colMeans(est_MC)^2), rep(NA, length(rho)-1)), sqrt(colMeans(est_MCP^2)-colMeans(est_MCP)^2))
RMSE <- sqrt(Bias^2+EP^2)
mean.std <- cbind(c(colMeans(std_glm), rep(NA, length(rho))), c(colMeans(std_MC,na.rm = T), rep(NA, length(rho)- 1)), colMeans(std_MCP))
rownames(Bias) <- rownames(EP) <- rownames(RMSE) <- rownames(mean.std) <- names(colMeans(est_MCP))
mat_Bias_RMSE <- cbind(Bias, RMSE); colnames(mat_Bias_RMSE) <- c('Bias - GLM', 'Bias - MC', 'Bias - MCP', 'RMSE - GLM', 'RMSE - MC', 'RMSE - MCP')

mat.EP <- cbind(EP[,1], mean.std[,1], EP[,2], mean.std[,2], EP[,3], mean.std[,3])
colnames(mat.EP) <- c('Emp. EP - GLM', 'Hes. EP - GLM', 'Emp. EP - MC', 'Hes. EP - MC', 'Emp. EP - MCP', 'Hes. EP - MCP')

### RESULTADOS PARA GUARDAR
paste0('n = ', n, ', rho = (', paste(rho, collapse = ', '), ')')
PFN_IN_GLM <- sum(pfn.insample.glm > pfn.target)/REP
PFN_IN_MC <- sum(pfn.insample.MC > pfn.target)/REP
PFN_IN_MCP <- sum(pfn.insample.MCP > pfn.target)/REP
PFN_IN_GLM; PFN_IN_MC; PFN_IN_MCP
pf_mat
mat_Bias_RMSE
mat.EP

save(list = c('PFN_IN_GLM', 'PFN_IN_MC', 'PFN_IN_MCP', 'pf_mat', 'mat_Bias_RMSE', 'mat.EP'), file = paste0('/home/andrey/Projetos/PerLog/Dados/Metricas/', arg[length(arg)], '/resultado.rds', sep = ''))
# load(paste0('/home/andrey/Projetos/PerLog/Dados/Metricas/C/resultado.rds', sep = ''))
