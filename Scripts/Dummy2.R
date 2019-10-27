#Note ANDREY
folder  <-  "C:/Users/Andrey/Projetos/PerLog/Scripts/"
# PC TRABALHO
#folder <- 'C:/Users/Andrey/Documents/Andrey/TCC/'
# PC Casa Alessandro
#folder <- '/home/alessandro/Dropbox/alessandro/2018_2/Orientacao_Monografia_Katiuski/TCCKatiuski/Monografia/codigo/'


LOOP <- expand.grid(n = c(100, 400), rho1 = c(.8, -.5), rho1 = c(-.4, .6), rho1 = c(.8, -.5))
A <- function(PAR){
  source('C:/Users/Andrey/Projetos/PerLog/Scripts/funcoes.txt')
  # simulacao
  n <- PAR[1]; k <- 5; rho <- unlist(PAR[-1]); beta <- rep(0.5, k);n.desc <- 50;m <- n+n.desc; pfn.target <- .05
  X1 <- cbind(1, matrix(runif(n = (m*(k-1)), min = -1, max = 1), m, (k-1)))
  colnames(X1) <- c('Intercept', kronecker("X_", 1:(k-1), paste, sep=''))
  X <- X1[1:n,]; Xprev <-  X1[(n+1):m,]
  REP <- 1000
  
  est_glm <- yhat_glm <- std_glm <- est_MC <- std_MC <- fit_glmvet <- corte0 <- corte1 <- corte <- NULL
  pfp.insample.MC <- pfp.insample.glm <- pfn.insample.MC <- pfn.insample.glm <- NULL
  pfp.glm <- pfn.glm <- pfp.MC <- pfn.MC <- NULL
  fit_glmvet <- fit_MCvet <- list()
  est_glm <- std_glm <- matrix(nrow = REP, ncol = k)
  est_MC <- std_MC <- matrix(nrow = REP, ncol = (k+length(rho)))
  Resultados <- list()
  
  r <- 1
  
  ct.fim.vet <- prob <-  fit_MCvet <- y.desc <- list()
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
    # fit_MC <- logistic_MC(resp = y, covar = X, par0 = c(beta, rho), trace = 0)
    fit_MC <- tryCatch(logistic_MC(resp = y, covar = X, par0 = c(beta, rho), trace = 0),
                       error = function(e){NULL})
    fit_MCvet[[r]] <- fit_MC
    if(!is.null(fit_MC)){
      est_glm[r,] <- fit_glm[,1]
      std_glm[r,] <- fit_glm[,2]
      est_MC[r,] <- fit_MC[,1]
      std_MC[r,] <- fit_MC[,2]
      
      prob[[r]] <- predi(beta.vet=fit_MCvet[[r]][1:k,1],
                         rho = fit_MCvet[[r]][(1:length(rho)+k),1], covar=X, resp=y)
      
      prob.glm  <-   predict(object = fit_glmvet[[r]], newx = X, type=c("response"))
      prob.MC  <-  prob[[r]]
      in.sample.glm <- min.pfp.glm(y = y, prob = prob.glm, pfn.target = pfn.target)
      pfp.insample.glm[r] <- in.sample.glm$pfp.min
      pfn.insample.glm[r] <- in.sample.glm$pfn.min
      corte[r] <- in.sample.glm$cortemin
      in.sample.MC <- min.pfp.MC(y = y, prob = prob.MC, pfn.target = pfn.target)
      pfp.insample.MC[r] <- in.sample.MC$pfp.min
      pfn.insample.MC[r] <- in.sample.MC$pfn.min
      corte0[r] <- in.sample.MC$corte0min
      corte1[r] <- in.sample.MC$corte1min
      
      pr.glm <- prev.glm(fit = fit_glmvet[[r]], newY = y.desc[[r]], newX = Xprev,
                         corte = corte[r])
      pfp.glm[r] <- pr.glm$pfp; pfn.glm[r] <- pr.glm$pfn
      pr.MC <- prev.MC(beta = fit_MCvet[[r]][1:k,1],
                       rho  = fit_MCvet[[r]][(k+length(rho)),1],
                       newY = y.desc[[r]], newX = Xprev,
                       LastY= y[n], LastX = X[n,], corte0 = corte0[r], corte1 = corte1[r])
      pfp.MC[r] <- pr.MC$pfp; pfn.MC[r] <- pr.MC$pfn
      
      r  <-  r+1
      cat('\n')
    }
    else{
      cat('Deu ruim!\n')
    }
  }
  
  
  pf_insample_mat <- rbind(c(mean(pfp.insample.glm), mean(pfp.insample.MC)),
                           c(mean(pfn.insample.glm), mean(pfn.insample.MC)))
  pf_out_mat <- rbind(c(mean(pfp.glm), mean(pfp.MC)),
                      c(mean(pfn.glm), mean(pfn.MC)))
  pf_mat <- rbind(pf_insample_mat, pf_out_mat)
  rownames(pf_mat) <- c('PFP - insample','PFN - insample', 'PFP - out','PFN - out')
  colnames(pf_mat) <- c('GLM','MC')
  
  Bias <- cbind(c(colMeans(est_glm)-beta, rep(NA, length(rho))), (colMeans(est_MC)-c(beta, rho)))
  EP <- cbind(c(sqrt(colMeans(est_glm^2)-colMeans(est_glm)^2), rep(NA, length(rho))), sqrt(colMeans(est_MC^2)-colMeans(est_MC)^2))
  RMSE <- sqrt(Bias^2+EP^2)
  mean.std <- cbind(c(colMeans(std_glm), rep(NA, length(rho))), colMeans(std_MC))
  rownames(Bias) <- rownames(EP) <- rownames(RMSE) <- rownames(mean.std) <- names(colMeans(est_MC))
  mat_Bias_RMSE <- cbind(Bias, RMSE); colnames(mat_Bias_RMSE) <- c('Bias - GLM','Bias - MC', 'RMSE - GLM','RMSE - MC')
  
  mat.EP <- cbind(EP[,1], mean.std[,1], EP[,2], mean.std[,2])
  colnames(mat.EP) <- c('Emp. EP - GLM', 'Hes. EP - GLM', 'Emp. EP - MC', 'Hes. EP - MC')
  
  ### RESULTADOS PARA GUARDAR
  resultado <- list('A' = list('PFN - In Sample - GLM' = sum(pfn.insample.glm>.05)/REP,
                               'PFN - In Sample - GLM' = sum(pfn.insample.MC>.05)/REP, 
                               'PF - Matriz' = pf_mat,
                               'Bias e RMSE' = mat_Bias_RMSE,
                               'Matriz de EP' = mat.EP ))
  names(resultado) <- paste0('n = ', n, " e rho's = (", paste(rho, collapse = ', '), ")")
  return(resultado)
}

set.seed(555)
cl <- makeCluster(getOption("cl.cores", 16))
system.time(RESULTS <- parApply(cl = cl, LOOP, 1, A))
stopCluster(cl)
