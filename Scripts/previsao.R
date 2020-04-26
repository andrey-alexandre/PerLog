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

sel_matr <-  as.matrix(read.table(paste(folder,"Dados/in_cov_WL7.txt", sep=''), header = TRUE))
sel_matrp <-  as.matrix(read.table(paste(folder,"Dados/out_cov_WL7.txt", sep=''), header = TRUE))

nomes <- c('Intercept',colnames(sel_matr))
sel_matr <- cbind(1, sel_matr); sel_matrp <- cbind(1, sel_matrp)
colnames(sel_matr) <- colnames(sel_matrp) <- nomes
sel_matr1 <- cbind(sel_matr, O3)
colnames(sel_matr1) <- c(colnames(sel_matr),'O3ar')

fit_glm <- glm(O3~.-1, data = as.data.frame(sel_matr1), family=binomial)
fit_MC <- logistic_MC(resp = O3, covar = sel_matr1, par0 = rep(0.5, ncol(sel_matr1) + 1))
fit_MCP <-logistic_MC(resp = O3, covar = sel_matr1, par0 = rep(0.5, ncol(sel_matr1) + 7))

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

minpfpGLM
prevGLM

