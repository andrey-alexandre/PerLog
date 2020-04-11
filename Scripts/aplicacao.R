#Note ANDREY
folder  <-  "/home/andrey/Projetos/PerLog/"
# PC CASA
#folder<-"/home/alessandro/Dropbox/alessandro/2018_2/Orientacao_Monografia_Katiuski/TCCKatiuski/Monografia/codigo/"
library(dplyr)
library(lubridate)
replace_na <- function(x){
  x <- ifelse(is.na(x), mean(x, na.rm = T), x)
  
  return(x)
}

dados <- 
  readxl::read_excel(paste0(folder,"Dados/dados_NA.xls")) %>% 
  mutate_all(replace_na)

y <- dados$E1_O3
X_pol <- dados %>% 
  select(Data, NO = E1_NO, NO2 = E1_NO2, CO = E1_CO) %>% 
  mutate(Data = ymd_hms(Data), Hora = hour(Data), 
         Data = floor_date(Data, unit = 'days')) %>% 
  filter(between(Hora, 7, 13)) %>% 
  tidyr::pivot_wider(names_from = Hora, values_from = c(NO, NO2, CO))

X_met <- dados %>% 
  select(Data, Tem = E2_Temperatura, WS = E2_ScalarWindSpeed) %>% 
  mutate(Data = ymd_hms(Data), Hora = hour(Data), 
         Data = floor_date(Data, unit = 'days')) %>% 
  filter(between(Hora, 7, 13)) %>% 
  group_by(Data) %>% 
  summarise(Tem = max(Tem), WS = max(WS))

X_oni <- 
  readr:: read_csv2(paste0(folder,"Dados/ONI.csv")) %>% 
    mutate(Data = ymd_hms(paste(YR, MON, '01 00:00:00', sep='/'))) %>% 
    select(Data, ANOM)


data.frame(row = rep(c(1, 51), each = 3),
           var = c("Sepal.Length", "Species", "Species_num"),
           value = c(5.1, "setosa", 1, 7.0, "versicolor", 2)) %>% tidyr::spread(var, value)
