#Note ANDREY
folder  <-  "/home/andrey/Projetos/PerLog/"
# PC CASA
#folder<-"/home/alessandro/Dropbox/alessandro/2018_2/Orientacao_Monografia_Katiuski/TCCKatiuski/Monografia/codigo/"

#################
### LIBRARIES ###
#################
library(dplyr)
library(lubridate)
replace_na <- function(x){
  x <- ifelse(is.na(x), mean(x, na.rm=T), x)
  
  return(x)
}
weekend <- function(date_){
  weekday <- weekdays(date_)
  weekend <- ifelse(weekday %in% c('sábado', 'domingo'), 'Weekend', 'Weekday')
  
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
  readxl::read_excel(paste0(folder,"Dados/dados_NA.xls"))

y <- dados$E1_O3 %>% replace_na
X_pol <- dados %>% 
  select(data_=Data, no=E1_NO, no2=E1_NO2, CO=E1_CO) %>% 
  mutate(data_=ymd_hms(data_), hour_=hour(data_), 
         data_=floor_date(data_, unit='days'),
         no=replace_na(no), no2=replace_na(no2),
         CO=replace_na(CO)) %>% 
  filter(between(hour_, 7, 13)) %>% 
  tidyr::pivot_wider(names_from=hour_, values_from=c(no, no2, CO)) %>% 
  mutate(join_date=month(data_))
X_met <- dados %>% 
  select(data_=Data, Tem=E2_Temperatura, WS=E2_ScalarWindSpeed) %>% 
  mutate(data_=ymd_hms(data_), hour_=hour(data_), 
         data_=floor_date(data_, unit='days'),
         Tem=replace_na(Tem), WS=replace_na(WS)) %>% 
  filter(between(hour_, 7, 13)) %>% 
  group_by(data_) %>% 
  summarise(Tem = max(Tem), WS = max(WS)) %>% 
  mutate(join_date=month(data_))
X_oni <- 
  readr:: read_csv2(paste0(folder,"Dados/ONI.csv")) %>% 
  mutate(data_=ymd_hms(paste(YR, MON, '01 00:00:00', sep='/'))) %>% 
  select(data_, oni=ANOM)
season <- get_season(X_pol$data_)
weekend_ <- weekend(X_pol$data_)

model.matrix(~(.)^2, X_pol, contrasts.arg = NULL)

matr<-cbind(matr, model.matrix(~(.)^2+ONI+ONI:season, X_met, contrasts.arg = NULL))
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