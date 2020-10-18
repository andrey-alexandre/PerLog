#Note ANDREY
folder  <-  "/home/rstudio/PerLog/"
# PC CASA
# folder<-"/home/alessandro/Dropbox/alessandro/2018_2/Orientacao_Monografia_Katiuski/TCCKatiuski/Monografia/codigo/"
# PC ALESSANDRO
# folder <- '/home/alessandro/Dropbox/alessandro/UFES/aulas/2019_1/Orientacao_IC_Andrey/'

#################
### LIBRARIES ###
#################
library(dplyr)
library(lubridate)
require(glmnet)
source(paste0(folder, 'Scripts/funcoes.R'))


replace_na <- function(x){
  x <- ifelse(is.na(x), mean(x, na.rm=T), x)
  
  return(x)
}
weekend <- function(date_){
  weekday <- weekdays(date_)
  weekend <- ifelse(weekday %in% c('Saturday', 'Sunday'), 'Weekend', 'Weekday')
  weekend <- factor(weekend, levels = c('Weekend', 'Weekday'), labels = c('Weekend', 'Weekday'))
  
  return(weekend)
}
get_season <- function(input_date){
  numeric_date <- 100 * month(input_date) + day(input_date)
  ## input Seasons upper limits in the form MMDD in the "break =" option:
  cuts <- base::cut(numeric_date, breaks = c(0,319,0620,0921,1220,1231)) 
  # rename the resulting groups (could've been done within cut(...levels=) if "Winter" wasn't double
  levels(cuts) <- c("Summer","Fall","Winter","Spring","Summer")
  return(cuts)
}


#####################
### PREPROCESSING ###
#####################
dados <- 
 read_csv2(paste0(folder,"Dados/Dados_NA.csv"))

y <- dados %>% 
  select(data_=Data, E1_O3) %>% 
  mutate(data_ = floor_date(ymd_hms(data_), 'day'),
         E1_O3 =  replace_na(E1_O3)) %>%
  group_by(data_) %>% 
  summarise(O3 = max(E1_O3)) %>%
  .$O3 %>% 
  (function(x) x>=80)
date_ <- dados %>% 
  select(data_=Data) %>% 
  mutate(data_ = floor_date(ymd_hms(data_), 'day')) %>% 
  distinct
X_pol <- dados %>% 
  select(data_=Data, no=E1_NO, no2=E1_NO2, CO=E1_CO) %>% 
  mutate(data_=ymd_hms(data_), hour_=hour(data_), 
         data_=floor_date(data_, unit='days'),
         no=replace_na(no), no2=replace_na(no2),
         CO=replace_na(CO)) %>% 
  filter(between(hour_, 7, 13)) %>% 
  group_by(data_) %>% summarise(no = max(no), no2 = max(no2), CO = max(CO)) %>% 
  select(-data_)
X_met <- dados %>% 
  select(data_=Data, Tem=E2_Temperatura, WS=E2_ScalarWindSpeed) %>% 
  mutate(data_=ymd_hms(data_), hour_=hour(data_), 
         data_=floor_date(data_, unit='days'),
         Tem=replace_na(Tem), WS=replace_na(WS)) %>% 
  filter(between(hour_, 7, 13)) %>% 
  group_by(data_) %>% 
  summarise(Tem = max(Tem), WS = max(WS)) %>% 
  select(-data_)
X_oni <- 
  readr:: read_csv2(paste0(folder,"Dados/ONI.csv")) %>% 
  mutate(data_=ymd_hms(paste(YR, MON, '01 00:00:00', sep='/')),
         ANOM=as.numeric(stringr::str_replace(ANOM, ",", "."))) %>% 
  select(data_, oni=ANOM) %>% 
  mutate(join_date=paste(month(data_), year(data_), sep='/')) %>% 
  right_join(date_ %>% 
               transmute(join_date=paste(month(data_), year(data_), sep='/')))

ONI <- X_oni$oni
season <- get_season(date_$data_)
week_ <- weekend(date_$data_)
# week_ <- weekdays(date_$data_); week_ <- factor(week_)
fl_tst <- year(date_$data_ ) > 2015

# Separa dados de treino e de tst
y_trn <- y[!fl_tst]; y_tst <- y[fl_tst]
w <- mean(y_trn); w <- ifelse(y_trn == 1, 1-w, w)
X_pol_trn <- X_pol[!fl_tst, ];X_pol_tst <- X_pol[fl_tst, ]
X_met_trn <- X_met[!fl_tst, ];X_met_tst <- X_met[fl_tst, ]
ONI_trn <- ONI[!fl_tst]; ONI_tst <- ONI[fl_tst]
season_trn <- season[!fl_tst]; season_tst <- season[fl_tst]
week_trn <- week_[!fl_tst]; week_tst <- week_[fl_tst]

write.table(y_trn, paste0(folder,"Dados/in_O3.txt"), col.names = TRUE, row.names = FALSE)
write.table(y_tst, paste0(folder,"Dados/out_O3.txt"), col.names = TRUE, row.names = FALSE)

#########################################
### CRIAÇÃO DE MATRIZ DE PLANEJAMENTO ###
#########################################
## TRAIN
X_pol_trn_i <- model.matrix(~ (.)^2, data = X_pol_trn)
X_pol_trn_i <- X_pol_trn_i[,-1]
X_pol_trn_i_s <- scale(X_pol_trn_i)
X_pol_trn_i_med <- attr(X_pol_trn_i_s, 'scaled:center'); X_pol_trn_i_sig <- attr(X_pol_trn_i_s, 'scaled:scale')

X_ONI_met_trn <- data.frame(oni=ONI_trn, X_met_trn)
X_ONI_met_trn_i <- model.matrix(~ (.)^2, data = X_ONI_met_trn)
X_ONI_met_trn_i <- X_ONI_met_trn_i[,-1]
X_ONI_met_trn_i_s <- scale(X_ONI_met_trn_i)
X_ONI_met_trn_i_med <- attr(X_ONI_met_trn_i_s, 'scaled:center'); X_ONI_met_trn_i_sig <- attr(X_ONI_met_trn_i_s, 'scaled:scale')

X_trn_df <- data.frame(week=week_trn, season=season_trn, X_pol_trn_i_s, X_ONI_met_trn_i_s)
X_trn <- model.matrix(~ . + week + week:no + week:no2 + week:CO +
                        season + season:oni + season:Tem + season:WS,
                      data = X_trn_df)
# colnames(X_trn)
# head(X_trn)

## TEST
X_pol_tst_i <- model.matrix(~ (.)^2, data = X_pol_tst)
X_pol_tst_i <- X_pol_tst_i[,-1]
X_pol_tst_i_s <- scale(X_pol_tst_i,
                       center = X_pol_trn_i_med, scale = X_pol_trn_i_sig)

X_ONI_met_tst <- data.frame(oni=ONI_tst, X_met_tst)
X_ONI_met_tst_i <- model.matrix(~ (.)^2, data = X_ONI_met_tst)
X_ONI_met_tst_i <- X_ONI_met_tst_i[,-1]
X_ONI_met_tst_i_s <- scale(X_ONI_met_tst_i,
                           center = X_ONI_met_trn_i_med, scale = X_ONI_met_trn_i_sig)

X_tst_df <- data.frame(week=week_tst, season=season_tst, X_pol_tst_i_s, X_ONI_met_tst_i_s)
X_tst <- model.matrix(~ . + week + week:no + week:no2 + week:CO +
                        season + season:oni + season:Tem + season:WS,
                      data = X_tst_df)
# colnames(X_tst)
# head(X_tst)

# # Matriz de treino
# matr <- model.matrix(~(.)^2, X_pol_trn, contrasts.arg = NULL)
# matr <- cbind(matr, model.matrix(~(.)^2+ONI_trn+ONI_trn:season_trn, X_met_trn, contrasts.arg = NULL))
# matr <- cbind(matr, model.matrix(~ONI_trn*(.), X_met_trn))
# matr <- cbind(model.matrix(~week_trn), model.matrix(~season_trn), ONI_trn, matr)
# matr <- cbind(model.matrix(~week_trn*., data = data.frame(scale(X_pol_trn))),
#               model.matrix(~season_trn*., data = data.frame(scale(X_met_trn))),
#               ONI_trn, matr)
# 
# # Matriz de tste
# matrp<-model.matrix(~(.)^2, X_pol_tst, contrasts.arg = NULL)
# matrp<-cbind(matrp, model.matrix(~(.)^2+ONI_tst+ONI_tst:season_tst, X_met_tst, contrasts.arg = NULL))
# matrp<-cbind(matrp, model.matrix(~ONI_tst*(.), X_met_tst))
# 
# matrp<-cbind(model.matrix(~week_tst), model.matrix(~season_tst), ONI_tst, matrp)
# matrp<-cbind(model.matrix(~week_tst*., data = data.frame(scale(X_pol_tst))),
#              model.matrix(~season_tst*., data = data.frame(scale(X_met_tst))),
#              ONI_tst, matrp)
# 
# colnames(matr) <- gsub('_trn', '', colnames(matr))
# colnames(matrp) <- colnames(matr)
# 
# matr<-matr[,!duplicated(colnames(matr))]
# scale_ <- caret::preProcess(matr[,-1]); matr[,-1] <- predict(scale_, matr[,-1])
# #write.table(matr, paste0(folder,"Dados/in_explan.txt"), col.names = TRUE, row.names = FALSE)
# 
# 
# matrp<-matrp[,!duplicated(colnames(matrp))]
# matrp[,-1] <- predict(scale_, matrp[,-1])
# #write.table(matrp, paste0(folder,"Dados/out_explan.txt"), col.names = TRUE, row.names = FALSE)

#############
### PESOS ###
#############
grupos <- month(date_$data_)

grupos_trn <- grupos[!fl_tst]; grupos_tst <- grupos[fl_tst]
############################
### SELEÇÃO DE VARIÁVEIS ###
############################
gam <- 1#Potência utilizada no cálculo dos w's
length_out <- 20
i <- 'UL'
for(i in c('UL', 'WL', 'PWL')){
  set.seed(555)
  if(i == 'UL'){
    ###### tradicional weights
    ols_fit <- glm(formula = y ~ .-1, family = 'binomial',
                   data = data.frame(y=as.numeric(y_trn), X_trn))
    w0 <- rep(1, length(y_trn))
    w3<-1/abs(coef(ols_fit))^gam
  }else if(i == 'WL'){
    ###### temporal weights
    w0 <- temporal.weights<-(1:length(y_trn)/length(y_trn))
    ols_fit <- glm(formula = y ~ .-1, family = 'binomial',
                   data = data.frame(y=as.numeric(y_trn), X_trn),
                   weights = w0)
    w3<-1/abs(coef(ols_fit))^gam
  }else if(i == 'PWL'){
    ###### proportional temporal weights
    weigth_aux <- 1/prop.table(table(y_trn))
    temporal.weights<-(1:length(y_trn)/length(y_trn))
    w0 <- ifelse(y_trn, temporal.weights*weigth_aux[2], temporal.weights*weigth_aux[1])
    ols_fit <- glm(formula = y ~ .-1, family = 'binomial',
                   data = data.frame(y=as.numeric(y_trn), X_trn),
                   weights = w0)
    w3<-1/abs(coef(ols_fit))^gam
  }
  fit_glm <- cv.glmnet(x=X_trn[,-1], y=y_trn, family = 'binomial', alpha=1, weights = w0,
                       parallel = FALSE, standardize = FALSE, foldid = grupos_trn,
                       type.measure = 'auc', penalty.factor = w3)
  coef_glm <- coef(fit_glm, s=fit_glm$lambda.min)[,1]
  write.table(X_trn[,names(which(coef_glm[1:ncol(X_trn)]!=0))], paste0(folder,"Dados/in_explain_RL_", i, ".txt"), col.names = TRUE, row.names = FALSE)
  write.table(X_tst[,names(which(coef_glm[1:ncol(X_trn)]!=0))], paste0(folder,"Dados/out_explain_RL_", i, ".txt"), col.names = TRUE, row.names = FALSE)

  lamb_MC <- seq(10^-3 , 2/(.4+sum(abs(coef(ols_fit)))), length.out = length_out)
  fit_MC <- pMC <- minpfpMC <- list()
  for(j in 1:length_out){
    fit_MC[[j]] <- logistic_MC(resp = as.numeric(y_trn), covar = X_trn, lamb = lamb_MC[j],
                               par0 = c(coef(ols_fit), 0), w0=w0, alpha=c(w3, 1), trace = 1)
    pMC[[j]] <- predi(beta.vet=fit_MC[[j]][1:length(fit_glm$coefficients)],
                      rho=fit_MC[[j]][length(fit_glm$coefficients)+1], covar=X_trn, resp=as.numeric(y_trn))
    minpfpMC[[j]] <- min.pfp.MC(y=as.numeric(y_trn), prob=pMC[[j]], esp.target = .9, s=1)
  }

  
  coef_MC <- fit_MC[[which.max(bind_rows(minpfpMC)$acc)]][,1]
  write.table(X_trn[,names(which(coef_MC[1:ncol(X_trn)]!=0))], paste0(folder,"Dados/in_explain_RLA_", i, ".txt"), col.names = TRUE, row.names = FALSE)
  write.table(X_tst[,names(which(coef_MC[1:ncol(X_trn)]!=0))], paste0(folder,"Dados/out_explain_RLA_", i, ".txt"), col.names = TRUE, row.names = FALSE)
  
  
  lamb_MCP <- seq(10^-3, 2/(7*0.4+sum(abs(coef(ols_fit)))), length.out = length_out)
  fit_MCP <- pMCP <- minpfpMCP <- list()
  for(j in 1:length_out){
    fit_MCP[[j]] <-logistic_MC(resp = as.numeric(y_trn), covar = X_trn, lamb = lamb_MCP[j],
                          par0 = c(coef(ols_fit), rep(0, 7)), w0=w0,
                          alpha=c(w3, rep(1, 7)), trace = 1)
    pMCP[[j]] <- predi(beta.vet=fit_MCP[[j]][1:length(fit_glm$coefficients)],
                      rho=fit_MCP[[j]][length(fit_glm$coefficients)+1:7], covar=X_trn, resp=as.numeric(y_trn))
    minpfpMCP[[j]] <- min.pfp.MC(y=as.numeric(y_trn), prob=pMCP[[j]], esp.target = .9, s=7)
  }

  coef_MCP <- fit_MCP[[which.max(rbind(lapply(minpfpMCP, function(x) x$acc)))]][,1]
  write.table(X_trn[,names(which(coef_MCP[1:ncol(X_trn)]!=0))], paste0(folder,"Dados/in_explain_RLAP_", i, ".txt"), col.names = TRUE, row.names = FALSE)
  write.table(X_tst[,names(which(coef_MCP[1:ncol(X_trn)]!=0))], paste0(folder,"Dados/out_explain_RLAP_", i, ".txt"), col.names = TRUE, row.names = FALSE)
  
  prmatrix(cbind(c(coef_glm, rep(NA, 7)), 
                 c(coef_MC, rep(NA, 6)),
                 coef_MCP))
}

