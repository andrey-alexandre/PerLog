#Note ANDREY
folder  <-  "/home/andrey/Projetos/PerLog/"
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

sel_matr <-  as.matrix(read.table(paste(folder,"Dados/in_explan.txt", sep=''), header = TRUE))
sel_matrp <-  as.matrix(read.table(paste(folder,"Dados/out_explan.txt", sep=''), header = TRUE))

nomes <- c('Intercept',colnames(sel_matr))
sel_matr <- cbind(1, sel_matr); sel_matrp <- cbind(1, sel_matrp)
colnames(sel_matr) <- colnames(sel_matrp) <- nomes

v <- 1:ncol(sel_matr)
sel_matr1 <- cbind(sel_matr[,], O3)
colnames(sel_matr1) <- c(colnames(sel_matr[,v]),'O3ar')
sel_matr1

grupos_train <- as.matrix(read.table(paste(folder,"Dados/in_group.txt", sep=''), header = TRUE))
grupos_test <- as.matrix(read.table(paste(folder,"Dados/out_group.txt", sep=''), header = TRUE))

for(i in c('UL', 'WL', 'PWL')){
  set.seed(555)
  if(i == 'UL'){
    ###### tradicional weights 
    ridge_unweighted <- cv.glmnet(x=sel_matr, y=O3, family = 'binomial', alpha=0, foldid = grupos_train,
                                  parallel=FALSE, standardize=TRUE)
    w3<-1/abs(matrix(coef(ridge_unweighted, s=ridge_unweighted$lambda.min)[,1][2:(ncol(sel_matr)+1)]))^1 ## Using gamma = 1
    w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
  }else if(i == 'WL'){
    ###### temporal weights  
    temporal.weights<-(1:length(O3)/length(O3))
    ridge_weighted<- cv.glmnet(x=sel_matr, y=O3, family = 'binomial', alpha=0, weights = temporal.weights,
                               type.measure = 'auc', foldid = grupos_train,
                               parallel=FALSE, standardize=TRUE)
    w3<-1/abs(matrix(coef(ridge_weighted, s=ridge_weighted$lambda.min)[,1][2:(ncol(sel_matr)+1)]))^1 ## Using gamma = 1
    w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
  }else if(i == 'WL'){
    ###### proportional temporal weights 
    weigth_aux <- 1/prop.table(table(O3))
    temporal_weights<-(1:length(O3)/length(O3))
    temporal_weights <- ifelse(O3, temporal_weights*weigth_aux[2], temporal_weights*weigth_aux[1])
    ridge_temporal_weighted<- cv.glmnet(x=sel_matr, y=O3, family = 'binomial', alpha=0, weights = temporal_weights,
                                        type.measure = 'auc', foldid = grupos_train,
                                        parallel=FALSE, standardize=TRUE)
    w3<-1/abs(matrix(coef(ridge_temporal_weighted, s=ridge_temporal_weighted$lambda.min)[,1][2:(ncol(sel_matr)+1)]))^1 ## Using gamma = 1
    w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
  }
  
  fit_glm <- cv.glmnet(x=sel_matr, y=O3, family = 'binomial', alpha=1,
                       parallel = FALSE, standardize = TRUE, foldid = grupos_train,
                       type.measure = 'auc', penalty.factor = w3)
  fit_MC <- logistic_MC(resp = O3, covar = sel_matr[,-1], par0 = c(as.vector(t(coef(fit_glm))), 0), w0=w, alpha=w3)
  fit_MCP <-logistic_MC(resp = O3, covar = sel_matr[,-1], par0 = c(as.vector(t(coef(fit_glm))), rep(0, 7)), w0=w, alpha=w3)
  
  # out <- list(covar=sel_matr1, resp=O3, li=0.000000001, ls=0.999999999)
  # par0 = rep(0.5, ncol(sel_matr1) + 7)
  # loglik(par0, out)
  
  #################################################
  ################## FORECASTING ##################
  #################################################
  pfn_target <- .1
  
  pglm <- predict(object = fit_glm, newx = sel_matrp, type = 'response')
  minpfpGLM <- min.pfp.glm(y=O3, prob=pglm, sen.target = pfn_target)
  prevGLM <- prev.glm(fit=fit_glm, newY=as.matrix(O3p), newX=sel_matrp, corte=minpfpGLM$c)
  
  pMC <- predi(beta.vet=fit_MC[1:length(fit_glm$coefficients)],
               rho=fit_MC[length(fit_glm$coefficients)+1], covar=sel_matr, resp=O3)
  minpfpMC <- min.pfp.MC(y=O3, prob=pMC, pfn.target = pfn_target, s=1)
  prevMC <- prev.MC(beta=fit_MC[1:length(fit_glm$coefficients)],
                    rho  = fit_MC[length(fit_glm$coefficients)+1],
                    newY = O3p, newX = sel_matrp,
                    y= O3[2191], X = sel_matr[2191,], corte0 = minpfpMC$c0, corte1 = minpfpMC$c1)
  
  pMCP <- predi(beta.vet=fit_MCP[1:length(fit_glm$coefficients)],
                rho=fit_MCP[length(fit_glm$coefficients)+1:7], covar=sel_matr, resp=O3)
  minpfpMCP <- min.pfp.MC(y=O3, prob=pMCP, pfn.target = pfn_target, s=7)
  prevMCP <- prev.MC(beta=fit_MCP[1:length(fit_glm$coefficients)],
                     rho  = fit_MCP[length(fit_glm$coefficients)+1:7],
                     newY = O3p, newX = sel_matrp,
                     y= O3[2191], X = sel_matr[2191, ], corte0 = minpfpMCP$c0, corte1 = minpfpMCP$c1)
  
  data.frame(Peso=i,
             Modelo = c('RL', 'RLA', 'RLAP'), 
             PFP = c(prevGLM$pfp, prevMC$pfp, prevMCP$pfp), 
             PFN = c(prevGLM$pfn, prevMC$pfn, prevMCP$pfn),
             ESP = c(prevGLM$esp, prevMC$esp, prevMCP$esp),
             SEN = c(prevGLM$sen, prevMC$sen, prevMCP$sen),
             ACC = c(prevGLM$acc, prevMC$acc, prevMCP$acc)) %>% 
    mutate_if(is.numeric, function(x) round(x, 3))%>%
    write.table(paste0(folder,"Dados/out_result.txt"), col.names = TRUE, row.names = FALSE, append = T)
}

