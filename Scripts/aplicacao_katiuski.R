library(readxl)
require(corrplot)
require(glmnet)
library(readr)

#Note ANDREY
folder  <-  "/home/andrey/Projetos/PerLog/"
# PC CASA
#folder<-"/home/alessandro/Dropbox/alessandro/2018_2/Orientacao_Monografia_Katiuski/TCCKatiuski/Monografia/codigo/"


dados <- read_excel(paste(folder,"Dados/dados_Imputed.xlsm", sep=''))
names(dados)
attach(dados)

year<-c((2010+(1:365-1)/365), (2011+(1:365-1)/365), (2012+(1:366-1)/366),
        (2013+(1:365-1)/365), (2014+(1:365-1)/365), (2015+(1:365-1)/365),
        (2016+(1:366-1)/366))

# d<-366+365
# a<-366
v0<-which(year==2010)
vf<-which(year==2015)-1
prever<-length(year)

# vf<-length(year)

MatO3<-matrix(dados$E1_O3, ncol=24, byrow = TRUE)
# O3MaxMorn<-apply(MatO3[,1:2], MARGIN = 1, max)
#O3MaxMorn<-MatO3[,8]
HourlyMaxO3<-apply(MatO3, MARGIN = 1, which.max)

pdf(file = paste(folder,'Imagem/O3max.pdf', sep=''),width = 20, height = 8)
par(cex.axis=3, cex.lab=1, cex.main=1.2, cex.sub=1)
barplot(table(HourlyMaxO3))
dev.off()

y.h<-dados$E1_O3
mat.y<-matrix(y.h, ncol = 24, byrow = TRUE)
y.day<-apply(X = mat.y, MARGIN = 1, FUN = max)

y.f<-y.day[v0:prever]
N<-length(y.f); n<-(vf-v0+1); d<-N-n
y<-y.f[1:n];y.d<-y.f[(n+1):N]
O3<-as.numeric(y>=80)
O3p<-as.numeric(y.d>=80)

write.table(O3, paste(folder,"Dados/in_O3.txt", sep=''), col.names = FALSE, row.names = FALSE)
write.table(O3p, paste(folder,"Dados/out_O3.txt", sep=''), col.names = FALSE, row.names = FALSE)

# Morning hours
horas<-7:13
# Precursors
pol.h<-dados[,c(5:6,8)]
pol.day<-pol.day.sd<-NULL
for(v in 1:dim(pol.h)[2]){
  mat.pol<-matrix(as.vector(as.matrix(pol.h[,v])), ncol = 24, byrow = TRUE)
  pol.day<-cbind(pol.day, mat.pol[,horas])
  pol.day.sd<-cbind(pol.day.sd, apply(X = mat.pol[,horas], MARGIN = 1, FUN = sd ))
}
pol.f<-as.matrix(pol.day[v0:prever,])
#pol.sd.f<-as.matrix(pol.day.sd[v0:prever,])
pol<-as.matrix(pol.f[1:n,]);pol.d<-as.matrix(pol.f[(n+1):N,])
#pol.sd<-as.matrix(pol.sd.f[1:n,]);pol.sd.d<-as.matrix(pol.sd.f[(n+1):N,])
# colnames(pol)<-colnames(pol.f)<-colnames(pol.d)<-kronecker(colnames(pol.h), horas, FUN = paste, sep='_h')
#colnames(pol.sd)<-colnames(pol.sd.f)<-colnames(pol.sd.d)<-paste(colnames(pol.h),'sd',sep='')
pol<-as.data.frame(pol)
#pol.sd<-as.data.frame(pol.sd)
# Meteorological variables
horas<-7:13
met.h<-dados[,c(12,14)]
met.day<-met.day.sd<-NULL
for(v in 1:dim(met.h)[2]){
  mat.met<-matrix(as.vector(as.matrix(met.h[,v])), ncol = 24, byrow = TRUE)
  met.day<-cbind(met.day, apply(X = mat.met[,horas], MARGIN = 1, FUN = max))
#  met.day.sd<-cbind(met.day.sd, apply(X = mat.met[,horas], MARGIN = 1, FUN = sd))
}
met.f<-as.matrix(met.day[v0:prever,])
#met.sd.f<-as.matrix(met.day.sd[v0:prever,])
met<-as.matrix(met.f[1:n,]);met.d<-as.matrix(met.f[(n+1):N,])
#met.sd<-as.matrix(met.sd.f[1:n,]);met.sd.d<-as.matrix(met.sd.f[(n+1):N,])
# colnames(met)<-colnames(met.f)<-kronecker(colnames(met.h), horas, FUN = paste, sep='_h')
#colnames(met.sd)<-colnames(met.sd.f)<-colnames(met.sd.d)<-paste(colnames(met.h),'sd',sep='')
met<-as.data.frame(met)
#met.sd<-as.data.frame(met.sd)
# Naming the covariates
#kronecker(c('SO2', 'NO', 'NO2', 'NOx', 'CO'), horas, FUN = paste, sep='_h')
colnames(pol)<-colnames(pol.d)<-kronecker(c('NO', 'NO2', 'CO'), horas, FUN = paste, sep='_h')
#colnames(pol.sd)<-paste(c('SO2', 'NO', 'NO2', 'NOx', 'CO'), 'sd', sep='')
# colnames(met)<-kronecker(c('Prec', 'Tem', 'RH', 'WS'), horas, FUN = paste, sep='_h')
colnames(met)<-colnames(met.d)<-c('Tem', 'WS')
#colnames(met.sd)<-paste(c('Prec', 'Tem', 'RH', 'WS'), 'sd', sep='')

M<-cor(cbind(pol,met))
corrplot(M, method = 'color', tl.pos="n", add = FALSE)
corrplot(M, method = 'color', add = FALSE)

#ONI
#ONImat<- read.csv("c:/Users/alessandro/Dropbox/TCCKatiuski/Regressao_logistica/ONI.csv",
#               header = FALSE,sep = ',')
# ONImat <- read_delim("C:/Users/Katiuski  Pereira/Dropbox/TCCKatiuski/Regressao_logistica/ONI.csv",
#                     delim = ',', col_names = FALSE)
#ONImat<-read.csv("C:/Users/LEA ADM/Desktop/katiuski/Monografia/codigo/ONI.csv",sep = ',', 
#                 header = FALSE)
# ONImat <- read.csv("~/Mega Shared/TCCKatiuski/Regressao_logistica/ONI.csv", sep = ',',
#                    header = FALSE)
ONImat <- read.csv2(paste(folder,"Dados/ONI.csv", sep = ''))

ONI.f <- rep(ONImat[,5],c(c(31,28,31,30,31,30,31,31,30,31,30,31),
                         c(31,28,31,30,31,30,31,31,30,31,30,31),
                         c(31,29,31,30,31,30,31,31,30,31,30,31),
                         c(31,28,31,30,31,30,31,31,30,31,30,31),
                         c(31,28,31,30,31,30,31,31,30,31,30,31),
                         c(31,28,31,30,31,30,31,31,30,31,30,31),
                         c(31,29,31,30,31,30,31,31,30,31,30,31)))
ONI<-ONI.f[1:n];ONI.d<-ONI.f[(n+1):N]
a<- cbind(ONI,year[v0:vf],O3)

#dev.off()
pdf(paste(folder,'/Imagem/oniano.pdf',sep=''),width = 20, height = 7)
par(cex.axis=3, cex.main=.2, cex.sub=.5)
plot(ONI ~ year[v0:vf],type="l",xlab="",ylab="",ylim=c(-2,3)) # xlab="Year",ylab="ONI"
points(a[,3][O3==1]+2~a[,2][O3==1],pch=16, col="red")
lines(c(0, 2191), c(0, 0), lty = 3, lwd = 2,col=8)
legend(2014, 2.5, c(expression(Y[t]==1), "ONI"),fill = c(2,1), bty="n" , cex=1.2)
dev.off()

# #### WL7 and UL7
# dia <- c("Fri","Sat","Sun","Mon","Tue","Wed","Thu")
# week.f<-rep(dia, N)[v0:prever]
# week.f<-factor(week.f, levels = dia)
# week<-week.f[1:n];week.d<-week.f[(n+1):N]

##### WL2 and UL2
dia <- c("Wk","End","End","Wk","Wk","Wk","Wk")
week.f<-rep(dia, N)[v0:prever]
week.f<-factor(week.f, levels = dia[1:2])
week<-week.f[1:n];week.d<-week.f[(n+1):N]

Est<-c(c("Sum","Aut","Win","Spr","Sum"),c("Aut","Win","Spr","Sum"),
       c("Aut","Win","Spr","Sum"),c("Aut","Win","Spr","Sum"),
       c("Aut","Win","Spr","Sum"),c("Aut","Win","Spr","Sum"),
       c("Aut","Win","Spr","Sum"))
season.f <- rep(Est,c(c(78,93,94,89,89),c(93,94,90,89),c(93,93,90,89),
                       c(93,93,90,89),c(93,93,90,89),c(93,94,90,89),c(92,92,92,11)))
season.f <-factor(season.f, levels = c(c("Sum","Aut","Win","Spr")))
season.f <-season.f[v0:prever]
season<-season.f[1:n];season.d<-season.f[(n+1):N]


#####################
#Regressao Logistica# ALESSANDRO
#####################
# dados1<-data.frame(pol)
dados1<-data.frame(pol)
#matriz de planejamento
matr<-model.matrix(~(.)^2, dados1, contrasts.arg = NULL)
# dados1<-data.frame(met)
dados1<-data.frame(met)
matr<-cbind(matr, model.matrix(~(.)^2+ONI+ONI:season, dados1, contrasts.arg = NULL))
matr<-cbind(matr, model.matrix(~ONI*(.), dados1))

matr<-cbind(model.matrix(~week), model.matrix(~season), ONI, matr)
matr<-cbind(model.matrix(~week*., data = data.frame(scale(pol))),
             model.matrix(~season*., data = data.frame(scale(met))),
             ONI, matr)

dados1<-data.frame(pol.d)
#matriz de planejamento
matrp<-model.matrix(~(.)^2, dados1, contrasts.arg = NULL)
# dados1<-data.frame(met)
dados1<-data.frame(met.d)
matrp<-cbind(matrp, model.matrix(~(.)^2+ONI.d+ONI.d:season.d, dados1, contrasts.arg = NULL))
matrp<-cbind(matrp, model.matrix(~ONI.d*(.), dados1))

matrp<-cbind(model.matrix(~week.d), model.matrix(~season.d), ONI.d, matrp)
matrp<-cbind(model.matrix(~week.d*., data = data.frame(scale(pol.d))),
            model.matrix(~season.d*., data = data.frame(scale(met.d))),
            ONI.d, matrp)
colnames(matrp)<-colnames(matr)
matr<-matr[,!duplicated(colnames(matr))]
matr<-matr[,-1]
matrp<-matrp[,!duplicated(colnames(matrp))]
matrp<-matrp[,-1]

#dm<-model.matrix(~week+season)
#dm<-cbind(dm, model.matrix(~. + .:., data=data.frame(pol, met)))
#dm<-cbind(dm, model.matrix(~. + season:., data=data.frame(ONI)))#seasons = esta??es do ano
#matr<-dm
#smatr<-scale(matr)


#######
#pesos#
#######

grupos.f<- c(
  c(rep(1, 31), rep(2, 28), rep(3, 31), rep(4,  30), rep(5,  31), rep(6,  30),
    rep(7, 31), rep(8, 31), rep(9, 30), rep(10, 31), rep(11, 30), rep(12, 31)),
  c(rep(1, 31), rep(2, 28), rep(3, 31), rep(4,  30), rep(5,  31), rep(6,  30),
    rep(7, 31), rep(8, 31), rep(9, 30), rep(10, 31), rep(11, 30), rep(12, 31)),
  c(rep(1, 31), rep(2, 29), rep(3, 31), rep(4,  30), rep(5,  31), rep(6,  30),
    rep(7, 31), rep(8, 31), rep(9, 30), rep(10, 31), rep(11, 30), rep(12, 31)),
  c(rep(1, 31), rep(2, 28), rep(3, 31), rep(4,  30), rep(5,  31), rep(6,  30),
    rep(7, 31), rep(8, 31), rep(9, 30), rep(10, 31), rep(11, 30), rep(12, 31)),
  c(rep(1, 31), rep(2, 28), rep(3, 31), rep(4,  30), rep(5,  31), rep(6,  30),
    rep(7, 31), rep(8, 31), rep(9, 30), rep(10, 31), rep(11, 30), rep(12, 31)),
  c(rep(1, 31), rep(2, 28), rep(3, 31), rep(4,  30), rep(5,  31), rep(6,  30),
    rep(7, 31), rep(8, 31), rep(9, 30), rep(10, 31), rep(11, 30), rep(12, 31)),
  c(rep(1, 31), rep(2, 29), rep(3, 31), rep(4,  30), rep(5,  31), rep(6,  30),
    rep(7, 31), rep(8, 31), rep(9, 30), rep(10, 31), rep(11, 30), rep(12, 31)))
grupos.f<-grupos.f[v0:prever]
grupos<-grupos.f[1:n]

# grupos<- c(
#   c(rep(1,  31), rep(2,  28), rep(3,  31), rep(4,  30), rep(5, 31),  rep(6,  30),
#     rep(7,  31), rep(8,  31), rep(9,  30), rep(10, 31), rep(11, 30), rep(12, 31)), 
#   c(rep(13, 31), rep(14, 28), rep(15, 31), rep(16, 30), rep(17, 31), rep(18, 30),
#     rep(19, 31), rep(20, 31), rep(21, 30), rep(22, 31), rep(23, 30), rep(24, 31)),
#   c(rep(25, 31), rep(26, 29), rep(27, 31), rep(28, 30), rep(29, 31), rep(30, 30),
#     rep(31, 31), rep(32, 31), rep(33, 30), rep(34, 31), rep(35, 30), rep(36, 31)),
#   c(rep(37, 31), rep(38, 28), rep(39, 31), rep(40, 30), rep(41, 31), rep(42, 30),
#     rep(43, 31), rep(44, 31), rep(45, 30), rep(46, 31), rep(47, 30), rep(48, 31)),
#   c(rep(49, 31), rep(50, 28), rep(51, 31), rep(52, 30), rep(53, 31), rep(54, 30),
#     rep(55, 31), rep(56, 31), rep(57, 30), rep(58, 31), rep(59, 30), rep(60, 31)),
#   c(rep(61, 31), rep(62, 28), rep(63, 31), rep(64, 30), rep(65, 31), rep(66, 30),
#     rep(67, 31), rep(68, 31), rep(69, 30), rep(70, 31), rep(71, 30), rep(72, 31)))


## Lasso
#set.seed(555)
#lasso <- cv.glmnet(x=matr, y=O3, family='binomial', alpha=1, parallel=TRUE, standardize=TRUE, type.measure='auc')
#plot(cv.lasso)
###########################################
########### Adaptive Lasso
###### without temporal weights #un
set.seed(555)
colnames(O3)<-('O3')

ridge <- cv.glmnet(x=matr, y=O3, family = 'binomial', alpha=0, foldid = grupos,
                   parallel=FALSE, standardize=TRUE)
w3<-1/abs(matrix(coef(ridge, s=ridge$lambda.min)[,1][2:(ncol(matr)+1)]))^1 ## Using gamma = 1
w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999

lasso <- cv.glmnet(x=matr, y=O3, family = 'binomial', alpha=1,
                   parallel = FALSE, standardize = TRUE, foldid = grupos,
                   type.measure = 'auc', penalty.factor = w3)

coef.lasso<-coef(lasso, s=lasso$lambda.min)[,1][2:(ncol(matr)+1)]
nomes<-names(which(coef.lasso!=0))
## Save selected covariates (in and out of sample)
write.table(matr[,nomes], paste(folder,"dados/in_cov_UL",nlevels(week),".txt", sep=''), col.names = TRUE, row.names = FALSE)
write.table(matrp[,nomes], paste(folder,"dados/out_cov_UL",nlevels(week),".txt", sep=''), col.names = TRUE, row.names = FALSE)

###### with temporal weights
set.seed(555)
temporal.weights<-(1:length(O3)/length(O3))
ridge<- cv.glmnet(x=matr, y=O3, family = 'binomial', alpha=0, weights = temporal.weights,
                  type.measure = 'auc', foldid = grupos,
                  parallel=FALSE, standardize=TRUE)
w3<-1/abs(matrix(coef(ridge, s=ridge$lambda.min)[,1][2:(ncol(matr)+1)]))^1 ## Using gamma = 1
w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999

lasso <- cv.glmnet(x = matr, y = O3, family = 'binomial', alpha=1, weights = temporal.weights,
                   parallel = FALSE, standardize = TRUE, foldid = grupos,
                   type.measure = 'auc', penalty.factor = w3)

coef.lasso<-coef(lasso, s=lasso$lambda.min)[,1][2:(ncol(matr)+1)]
nomes<-names(which(coef.lasso!=0))
## Save selected covariates (in and out of sample)
write.table(matr[,nomes], paste(folder,"dados/in_cov_WL",nlevels(week),".txt", sep=''), col.names = TRUE, row.names = FALSE)
write.table(matrp[,nomes], paste(folder,"dados/out_cov_WL",nlevels(week),".txt", sep=''), col.names = TRUE, row.names = FALSE)
#######  End of ALASSO
###############################
