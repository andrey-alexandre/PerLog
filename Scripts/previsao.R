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

tipo <- 'PWL2'
sel_matr <-  as.matrix(read.table(paste(folder,"Dados/in_cov_", tipo, ".txt", sep=''), header = TRUE))
sel_matrp <-  as.matrix(read.table(paste(folder,"Dados/out_cov_", tipo, ".txt", sep=''), header = TRUE))

nomes <- c('Intercept',colnames(sel_matr))
sel_matr <- cbind(1, sel_matr); sel_matrp <- cbind(1, sel_matrp)
colnames(sel_matr) <- colnames(sel_matrp) <- nomes

v <- 1:ncol(sel_matr)
sel_matr1 <- cbind(sel_matr[,v], O3)
colnames(sel_matr1) <- c(colnames(sel_matr[,v]),'O3ar')
sel_matr1

fit_glm <- glm(O3ar~.-1, data = as.data.frame(sel_matr1), family=binomial)
fit_MC <- logistic_MC(resp = O3, covar = sel_matr[,v], par0 = c(fit_glm$coefficients, 0), w0=w)
fit_MCP <-logistic_MC(resp = O3, covar = sel_matr[,v], par0 = c(fit_glm$coefficients, rep(0, 7)), w0=w)

# out <- list(covar=sel_matr1, resp=O3, li=0.000000001, ls=0.999999999)
# par0 = rep(0.5, ncol(sel_matr1) + 7)
# loglik(par0, out)

#################################################
################## FORECASTING ##################
#################################################
pfn_target <- .1

pglm <- predict(object = fit_glm, newdata = as.data.frame(sel_matr1), type = 'response')
minpfpGLM <- min.pfp.glm(y=O3, prob=pglm, pfn.target = pfn_target)
prevGLM <- prev.glm(fit=fit_glm, newY=as.matrix(O3p), newX=sel_matrp, corte=minpfpGLM$cortemin)

pMC <- predi(beta.vet=fit_MC[1:length(fit_glm$coefficients)],
             rho=fit_MC[length(fit_glm$coefficients)+1], covar=sel_matr[,v], resp=O3)
minpfpMC <- min.pfp.MC(y=O3, prob=pMC, pfn.target = pfn_target)
prevMC <- prev.MC(beta=fit_MC[1:length(fit_glm$coefficients)],
                   rho  = fit_MC[length(fit_glm$coefficients)+1],
                   newY = O3p, newX = sel_matrp[, c(1, 6:10)],
                   LastY= O3p[366], LastX = sel_matrp[366, c(1, 6:10)], corte0 = minpfpMC$corte0min, corte1 = minpfpMC$corte1min)

pMCP <- predi(beta.vet=fit_MCP[1:length(fit_glm$coefficients)],
             rho=fit_MCP[length(fit_glm$coefficients)+1:7], covar=sel_matr[,v], resp=O3)
minpfpMCP <- min.pfp.MC(y=O3, prob=pMCP, pfn.target = pfn_target)

prevMCP <- prev.MC(beta=fit_MCP[1:length(fit_glm$coefficients)],
                  rho  = fit_MCP[length(fit_glm$coefficients)+1:7],
                  newY = O3p, newX = sel_matrp[, c(1, 6:10)],
                  LastY= O3p[366], LastX = sel_matrp[366, c(1, 6:10)], corte0 = minpfpMCP$corte0min, corte1 = minpfpMCP$corte1min)


data.frame(Input = tipo, Modelo = c('GLM', 'MC', 'MCP'), 
           minPFP = c(minpfpGLM$pfp.min, minpfpMC$pfp.min, minpfpMCP$pfp.min), 
           PFP = c(prevGLM$pfp, prevMC$pfp, prevMCP$pfp), 
           PFN = c(prevGLM$pfn, prevMC$pfn, prevMCP$pfn),
           ACC = c(prevGLM$acc, prevMC$acc, prevMCP$acc)) %>% 
  write.csv2(paste0('Dados/', tipo, '.csv'))
