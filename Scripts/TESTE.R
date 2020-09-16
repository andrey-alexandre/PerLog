#Note ANDREY
folder  <-  "/home/andrey/Projetos/PerLog/"
# PC TRABALHO
#folder <- 'C:/Users/Andrey/Documents/Andrey/TCC/'
# PC Casa Alessandro
# folder <- '~/Dropbox/alessandro/UFES/aulas/2019_1/Orientacao_IC_Andrey/'
source(paste(folder, 'Scripts/funcoes.R', sep = ''))
arg <- commandArgs(trailingOnly = T)
# arg <- c(200, round(sin(2*pi*(1:7)/15)/2, 2), 50)
# arg <- c(280, 1, .4, .5, .2, .4, .3, .5, .45, 10, 'A')

# simulacao 
n <- as.numeric(arg[1]); k <- 2; rho <- as.numeric(arg[-c(1, 2,length(arg)-0:1)]);
beta <- c(rep(-1, k-1),  as.numeric(arg[2]));n.desc <- 7*40; m <- n+n.desc; pfn.target <- .05
# X1 <- cbind(1, matrix( runif(n = m*(k-1), min = -1.4, max = -0.6), m, (k-1)))
# X1 <- cbind(1, (.5*sin(2*pi*(1:m)/28)+runif(m,-1.25,-.75)))
# X1 <- cbind(1, matrix( seq(from = -2, to = 1, length.out = m*(k-1)), m, (k-1)))
# colnames(X1) <- c('Intercept', kronecker("X_", 1:(k-1), paste, sep=''))
# X <- X1[1:n,]; Xprev <-  X1[(n+1):m,]
REP <- as.numeric(arg[length(arg)-1])


print(beta)

#nohup Rscript --vanilla /home/andrey/Projetos/PerLog/Scripts/TESTE.R 280 .5 .4 .5 .2 .4 .3 .5 .45 1000 A > /home/andrey/Projetos/PerLog/A.txt &