#Note ANDREY
folder  <-  "/home/andrey/Projetos/PerLog/"
# PC CASA
#folder<-"/home/alessandro/Dropbox/alessandro/2018_2/Orientacao_Monografia_Katiuski/TCCKatiuski/Monografia/codigo/"

#################
### LIBRARIES ###
#################
library(dplyr)
library(lubridate)
require(glmnet)
replace_na <- function(x){
  x <- ifelse(is.na(x), mean(x, na.rm=T), x)
  
  return(x)
}
weekend <- function(date_){
  weekday <- weekdays(date_)
  weekend <- ifelse(weekday %in% c('sábado', 'domingo'), 'Weekend', 'Weekday')
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
  readxl::read_excel(paste0(folder,"Dados/dados_NA.xlsx"))

y <- dados %>% 
  select(data_=Data, E1_O3) %>% 
  mutate(data_ = floor_date(data_, 'day'),
         E1_O3 =  replace_na(E1_O3)) %>%
  group_by(data_) %>% 
  summarise(O3 = max(E1_O3)) %>%
  .$O3 %>% 
  (function(x) x>=80)
date_ <- dados %>% 
  select(data_=Data) %>% 
  mutate(data_ = floor_date(data_, 'day')) %>% 
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
  mutate(data_=ymd_hms(paste(YR, MON, '01 00:00:00', sep='/'))) %>% 
  select(data_, oni=ANOM) %>% 
  mutate(join_date=paste(month(data_), year(data_), sep='/')) %>% 
  right_join(date_ %>% 
               transmute(join_date=paste(month(data_), year(data_), sep='/')))

ONI <- X_oni$oni
season <- get_season(date_$data_)
week_ <- weekend(date_$data_)
# week_ <- weekdays(date_$data_); week_ <- factor(week_)
fl_test <- year(date_$data_ ) > 2015

# Separa dados de treino e de test
y_train <- y[!fl_test]; y_test <- y[fl_test]
X_pol_train <- X_pol[!fl_test, ];X_pol_test <- X_pol[fl_test, ]
X_met_train <- X_met[!fl_test, ];X_met_test <- X_met[fl_test, ]
ONI_train <- ONI[!fl_test]; ONI_test <- ONI[fl_test]
season_train <- season[!fl_test]; season_test <- season[fl_test]
week_train <- week_[!fl_test]; week_test <- week_[fl_test]

write.table(y_train, paste0(folder,"Dados/in_O3.txt"), col.names = TRUE, row.names = FALSE)
write.table(y_test, paste0(folder,"Dados/out_O3.txt"), col.names = TRUE, row.names = FALSE)

#########################################
### CRIAÇÃO DE MATRIZ DE PLANEJAMENTO ###
#########################################
# Matriz de treino
matr <- model.matrix(~(.)^2, X_pol_train, contrasts.arg = NULL)
matr <- cbind(matr, model.matrix(~(.)^2+ONI_train+ONI_train:season_train, X_met_train, contrasts.arg = NULL))
matr <- cbind(matr, model.matrix(~ONI_train*(.), X_met_train))
matr <- cbind(model.matrix(~week_train), model.matrix(~season_train), ONI_train, matr)
matr <- cbind(model.matrix(~week_train*., data = data.frame(scale(X_pol_train))),
              model.matrix(~season_train*., data = data.frame(scale(X_met_train))),
              ONI_train, matr)

# Matriz de teste
matrp<-model.matrix(~(.)^2, X_pol_test, contrasts.arg = NULL)
matrp<-cbind(matrp, model.matrix(~(.)^2+ONI_test+ONI_test:season_test, X_met_test, contrasts.arg = NULL))
matrp<-cbind(matrp, model.matrix(~ONI_test*(.), X_met_test))

matrp<-cbind(model.matrix(~week_test), model.matrix(~season_test), ONI_test, matrp)
matrp<-cbind(model.matrix(~week_test*., data = data.frame(scale(X_pol_test))),
             model.matrix(~season_test*., data = data.frame(scale(X_met_test))),
             ONI_test, matrp)

colnames(matr) <- gsub('_train', '', colnames(matr))
colnames(matrp) <- colnames(matr)

matr<-matr[,!duplicated(colnames(matr))]
matr<-matr[,-1]
scale_ <- caret::preProcess(matr); matr <- predict(scale_, matr)

matrp<-matrp[,!duplicated(colnames(matrp))]
matrp<-matrp[,-1]
matrp <- predict(scale_, matrp)

saveRDS(scale_, paste(folder,"Dados/scaling",nlevels(week_),".rds", sep=''))
#############
### PESOS ###
#############
grupos <- month(date_$data_)
grupos_train <- grupos[!fl_test]; grupos_test <- grupos[fl_test]

##############
### ALASSO ###
##############

ridge <- cv.glmnet(x=matr, y=y_train, family = 'binomial', alpha=0, foldid = grupos_train,
                   parallel=FALSE, standardize=TRUE)
w3<-1/abs(matrix(coef(ridge, s=ridge$lambda.min)[,1][2:(ncol(matr)+1)]))^1 ## Using gamma = 1
w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999

lasso <- cv.glmnet(x=matr, y=y_train, family = 'binomial', alpha=1,
                   parallel = FALSE, standardize = TRUE, foldid = grupos_train,
                   type.measure = 'auc', penalty.factor = w3)

coef.lasso<-coef(lasso, s=lasso$lambda.min)[,1][2:(ncol(matr)+1)]
nomes<-names(which(coef.lasso!=0))

## Save selected covariates (in and out of sample)
write.table(matr[,nomes], paste0(folder,"Dados/in_cov_UL",nlevels(week_),".txt"), col.names = TRUE, row.names = FALSE)
write.table(matrp[,nomes], paste0(folder,"Dados/out_cov_UL",nlevels(week_),".txt"), col.names = TRUE, row.names = FALSE)

###### with temporal weights  
set.seed(555)
temporal.weights<-(1:length(y_train)/length(y_train))
ridge<- cv.glmnet(x=matr, y=y_train, family = 'binomial', alpha=0, weights = temporal.weights,
                  type.measure = 'auc', foldid = grupos_train,
                  parallel=FALSE, standardize=TRUE)
w3<-1/abs(matrix(coef(ridge, s=ridge$lambda.min)[,1][2:(ncol(matr)+1)]))^1 ## Using gamma = 1
w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999

lasso <- cv.glmnet(x = matr, y = y_train, family = 'binomial', alpha=1, weights = temporal.weights,
                   parallel = FALSE, standardize = TRUE, foldid = grupos_train,
                   type.measure = 'auc', penalty.factor = w3)

coef.lasso<-coef(lasso, s=lasso$lambda.min)[,1][2:(ncol(matr)+1)]
nomes<-names(which(coef.lasso!=0))
## Save selected covariates (in and out of sample)
write.table(matr[,nomes], paste(folder,"Dados/in_cov_WL",nlevels(week_),".txt", sep=''), col.names = TRUE, row.names = FALSE)
write.table(matrp[,nomes], paste(folder,"Dados/out_cov_WL",nlevels(week_),".txt", sep=''), col.names = TRUE, row.names = FALSE)

###### with proportional temporal weights  
set.seed(555)
weigth.aux <- 1/prop.table(table(y_train))
temporal.weights<-(1:length(y_train)/length(y_train))
temporal.weights <- ifelse(y_train, temporal.weights*weigth.aux[2], temporal.weights*weigth.aux[1])
ridge<- cv.glmnet(x=matr, y=y_train, family = 'binomial', alpha=0, weights = temporal.weights,
                  type.measure = 'auc', foldid = grupos_train,
                  parallel=FALSE, standardize=TRUE)
w3<-1/abs(matrix(coef(ridge, s=ridge$lambda.min)[,1][2:(ncol(matr)+1)]))^1 ## Using gamma = 1
w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999

lasso <- cv.glmnet(x = matr, y = y_train, family = 'binomial', alpha=1, weights = temporal.weights,
                   parallel = FALSE, standardize = TRUE, foldid = grupos_train,
                   type.measure = 'auc', penalty.factor = w3)

coef.lasso<-coef(lasso, s=lasso$lambda.min)[,1][2:(ncol(matr)+1)]
nomes<-names(which(coef.lasso!=0))
## Save selected covariates (in and out of sample)
write.table(matr[,nomes], paste(folder,"Dados/in_cov_PWL",nlevels(week_),".txt", sep=''), col.names = TRUE, row.names = FALSE)
write.table(matrp[,nomes], paste(folder,"Dados/out_cov_PWL",nlevels(week_),".txt", sep=''), col.names = TRUE, row.names = FALSE)