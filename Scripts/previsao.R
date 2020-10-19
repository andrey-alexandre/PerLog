#Note ANDREY
folder  <-  "/home/rstudio/PerLog/"
# PC CASA
#folder<-"/home/alessandro/Dropbox/alessandro/2018_2/Orientacao_Monografia_Katiuski/TCCKatiuski/Monografia/codigo/"

#################
### LIBRARIES ###
#################
library(dplyr)
require(glmnet)
source(paste0(folder, 'Scripts/funcoes.R'))

###########################################################
################## PERFORMING PREDICTION ##################
###########################################################

O3 <- read.table(paste0(folder,"Dados/in_O3.txt"), header = T);O3 <- as.matrix(O3)
O3p <- read.table(paste0(folder,"Dados/out_O3.txt"), header = T);O3p <- as.matrix(O3p)
w <- mean(O3); w <- ifelse(O3 == 1, 1-w, w)
i <- 'PWL'

for(i in c('UL', 'WL', 'PWL')){
  
  sel_matr_RL <-  as.matrix(read.table(paste0(folder,"Dados/in_explain_RL_", i, ".txt"), header = TRUE))
  sel_matrp_RL <-  as.matrix(read.table(paste0(folder,"Dados/out_explain_RL_", i, ".txt"), header = TRUE))
  sel_matr_RLA <-  as.matrix(read.table(paste0(folder,"Dados/in_explain_RLA_", i, ".txt"), header = TRUE))
  sel_matrp_RLA <-  as.matrix(read.table(paste0(folder,"Dados/out_explain_RLA_", i, ".txt"), header = TRUE))
  
  nomes <- c('Intercept',colnames(sel_matr_RLA)[-1])
  colnames(sel_matr_RLA) <- colnames(sel_matrp_RLA) <- nomes
  
  sel_matr_RLAP <-  as.matrix(read.table(paste0(folder,"Dados/in_explain_RLAP_", i, ".txt"), header = TRUE, dec = ',')); sel_matr_RLAP <- apply(sel_matr_RLAP, 2, as.numeric)
  sel_matrp_RLAP <-  as.matrix(read.table(paste0(folder,"Dados/out_explain_RLAP_", i, ".txt"), header = TRUE))
  
  nomes <- c('Intercept',colnames(sel_matr_RLAP)[-1])
  colnames(sel_matr_RLAP) <- colnames(sel_matr_RLAP) <- nomes

  fit_glm <- cv.glmnet(x=sel_matr_RL, y=O3, family = 'binomial', alpha=1,
                       parallel = FALSE, standardize = TRUE,
                       type.measure = 'auc')
  fit_MC <- logistic_MC(resp = O3, covar = sel_matr_RLA, par0 = c(rep(.5, dim(sel_matr_RLA)[2]), 0), w0=w, lamb=0, alpha=1)
  fit_MCP <-logistic_MC(resp = O3, covar = sel_matr_RLAP, par0 = c(rep(.5, dim(sel_matr_RLAP)[2]), rep(0, 7)), w0=w, lamb=0, alpha=1)
  
  # out <- list(covar=sel_matr1, resp=O3, li=0.000000001, ls=0.999999999)
  # par0 = rep(0.5, ncol(sel_matr1) + 7)
  # loglik(par0, out)
  
  #################################################
  ################## FORECASTING ##################
  #################################################
  esp_target <- .8
  
  pglm <- predict(object = fit_glm, newx = sel_matrp_RL, type = 'response')
  minpfpGLM <- min.pfp.glm(y=O3, prob=pglm, esp.target = esp_target)
  prevGLM <- prev.glm(fit=fit_glm, newY=as.matrix(O3p), newX=sel_matrp_RL, corte=minpfpGLM$c)
  
  pMC <- predi(beta.vet=fit_MC[1:length(fit_glm$coefficients)],
               rho=fit_MC[length(fit_glm$coefficients)+1], covar=sel_matr_RLA, resp=O3)
  minpfpMC <- min.pfp.MC(y=O3, prob=pMC, esp.target = esp_target, s=1)
  prevMC <- prev.MC(beta=fit_MC[1:length(fit_glm$coefficients)],
                    rho  = fit_MC[length(fit_glm$coefficients)+1],
                    newY = O3p, newX = sel_matrp_RLA,
                    y= O3[2191], X = sel_matr_RLA[2191,], corte0 = minpfpMC$c0, corte1 = minpfpMC$c1)
  
  pMCP <- predi(beta.vet=fit_MCP[1:length(fit_glm$coefficients)],
                rho=fit_MCP[length(fit_glm$coefficients)+1:7], covar=sel_matr_RLAP, resp=O3)
  minpfpMCP <- min.pfp.MC(y=O3, prob=pMCP, esp.target = esp_target, s=7)
  prevMCP <- prev.MC(beta=fit_MCP[1:length(fit_glm$coefficients)],
                     rho  = fit_MCP[length(fit_glm$coefficients)+1:7],
                     newY = O3p, newX = sel_matrp_RLAP,
                     y= O3[2191], X = sel_matr_RLAP[2191, ], corte0 = minpfpMCP$c0, corte1 = minpfpMCP$c1)
  
  data.frame(Peso=i,
             Modelo = c('RL', 'RLA', 'RLAP'), 
             PFP = c(prevGLM$pfp, prevMC$pfp, prevMCP$pfp), 
             PFN = c(prevGLM$pfn, prevMC$pfn, prevMCP$pfn),
             ESP = c(prevGLM$esp, prevMC$esp, prevMCP$esp),
             SEN = c(prevGLM$sen, prevMC$sen, prevMCP$sen),
             ACC = c(prevGLM$acc, prevMC$acc, prevMCP$acc)) %>% 
    mutate_if(is.numeric, function(x) round(x, 3))%>%
    write.table(paste0(folder,"Dados/out_result.txt"), col.names = TRUE, row.names = FALSE, append = T)
  
  data.frame(Peso=i,
             Modelo = c('RL', 'RLA', 'RLAP'), 
             PFP = c(minpfpGLM$pfp, minpfpMC$pfp, minpfpMCP$pfp), 
             PFN = c(minpfpGLM$pfn, minpfpMC$pfn, minpfpMCP$pfn),
             ESP = c(minpfpGLM$esp, minpfpMC$esp, minpfpMCP$esp),
             SEN = c(minpfpGLM$sen, minpfpMC$sen, minpfpMCP$sen),
             ACC = c(minpfpGLM$acc, minpfpMC$acc, minpfpMCP$acc)) %>% 
    mutate_if(is.numeric, function(x) round(x, 3))%>%
    write.table(paste0(folder,"Dados/in_result.txt"), col.names = TRUE, row.names = FALSE, append = T)
}

