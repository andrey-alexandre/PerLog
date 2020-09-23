#Note ANDREY
folder  <-  "/home/andrey/Projetos/PerLog/"
# PC CASA
#folder<-"/home/alessandro/Dropbox/alessandro/2018_2/Orientacao_Monografia_Katiuski/TCCKatiuski/Monografia/codigo/"

library(dplyr)
library(tidyr)
library(ggplot2)
turn_greek <- function(string){
  stringr::str_replace(string, "beta", "\U03B2") %>% 
    stringr::str_replace("rho", "\U03C1") %>% 
    stringr::str_replace("_0", "\U2080") %>% 
    stringr::str_replace("_1", "\U2081") %>% 
    stringr::str_replace("_2", "\U2082") %>% 
    stringr::str_replace("_3", "\U2083") %>% 
    stringr::str_replace("_4", "\U2084") %>% 
    stringr::str_replace("_5", "\U2085") %>% 
    stringr::str_replace("_6", "\U2086") %>% 
    stringr::str_replace("_7", "\U2087")
}

# Boxplots ----
parameters <- tibble()
for(i in LETTERS[5:10]){
  parameters_RL <- readr::read_csv(paste0(folder, "Dados/Parciais/", i, "/parcial_est_glm.csv"), col_names = paste0('beta_', 0:1)) %>% mutate(Sim=i, Modelo="RL", rep=1:1000)
  parameters_RLA <- readr::read_csv(paste0(folder, "Dados/Parciais/", i, "/parcial_est_MC.csv"), col_names =  c(paste0('beta_', 0:1), 'rho_1')) %>% mutate(Sim=i, Modelo="RLA", rep=1:1000)
  parameters_RLAP <- readr::read_csv(paste0(folder, "Dados/Parciais/", i, "/parcial_est_MCP.csv"), col_names =  c(paste0('beta_', 0:1), paste0('rho_', 1:7))) %>% mutate(Sim=i, Modelo="RLAP", rep=1:1000)
  
  parameters <-  bind_rows(parameters, parameters_RLA, parameters_RL, parameters_RLAP)
}

# Beta ----
parameters  %>% 
  select(-(starts_with("rho"))) %>% 
  mutate_all(~tidyr::replace_na(.x, 0)) %>% 
  pivot_longer(-c(Sim, Modelo, rep)) %>% 
  mutate(Tamanho=case_when(Sim %in% c('E', 'G','I') ~ "n=280",
                           Sim %in% c('F', 'H','J') ~ "n=840"),
         Casos=case_when(Sim %in% c('E', 'F') ~ "Sem correlação",
                         Sim %in% c('G', 'H') ~ "Com correlação",
                         Sim %in% c('I', 'J') ~ "Com correlação periódica") %>% 
           factor(levels = c("Sem correlação", "Com correlação", "Com correlação periódica")),
         real_value = case_when(name == "beta_0" ~ -1, 
                                name == "beta_1" ~ .5,
                                Casos == "Sem correlação" ~ 0,
                                Casos =="Com correlação"~ .5,
                                name=="rho_1" ~ 0.6,
                                name=="rho_2" ~ 0.3, 
                                name=="rho_3" ~ 0.1,
                                name=="rho_4" ~ -0.2, 
                                name=="rho_5" ~ -0.4,
                                name=="rho_6" ~ -0.2, 
                                name=="rho_7" ~ 0.2)) %>% 
  ggplot(aes(x=Tamanho, y=value, col=Modelo))+
  geom_boxplot() + 
  geom_hline(aes(yintercept=real_value))+
  facet_grid(Casos~turn_greek(name), scales = 'free')

# Rho ----
parameters  %>% 
  select(-(starts_with("beta"))) %>% 
  mutate_at(vars(matches("rho")), ~ifelse(is.na(.x) & Modelo =="RLA", rho_1, .x)) %>% 
  mutate_all(~tidyr::replace_na(.x, 0)) %>% 
  pivot_longer(-c(Sim, Modelo, rep)) %>% 
  mutate(Tamanho=case_when(Sim %in% c('E', 'G','I') ~ "n=280",
                           Sim %in% c('F', 'H','J') ~ "n=840"),
         Casos=case_when(Sim %in% c('E', 'F') ~ "Sem correlação",
                         Sim %in% c('G', 'H') ~ "Com correlação",
                         Sim %in% c('I', 'J') ~ "Com correlação periódica") %>% 
           factor(levels = c("Sem correlação", "Com correlação", "Com correlação periódica")),
         real_value = case_when(name == "beta_0" ~ -1, 
                                name == "beta_1" ~ .5,
                                Casos == "Sem correlação" ~ 0,
                                Casos =="Com correlação"~ .5,
                                name=="rho_1" ~ 0.6,
                                name=="rho_2" ~ 0.3, 
                                name=="rho_3" ~ 0.1,
                                name=="rho_4" ~ -0.2, 
                                name=="rho_5" ~ -0.4,
                                name=="rho_6" ~ -0.2, 
                                name=="rho_7" ~ 0.2)) %>% 
  ggplot(aes(x=Tamanho, y=value, col=Modelo))+
  geom_boxplot() + 
  geom_hline(aes(yintercept=real_value))+
  facet_grid(Casos~turn_greek(name), scales = 'free')


# Barplots ----
# RMSE - Beta ----
parameters  %>% 
  select(-(starts_with("rho"))) %>% 
  mutate_all(~tidyr::replace_na(.x, 0)) %>% 
  pivot_longer(-c(Sim, Modelo, rep)) %>% 
  mutate(Tamanho=case_when(Sim %in% c('E', 'G','I') ~ "n=280",
                           Sim %in% c('F', 'H','J') ~ "n=840"),
         Casos=case_when(Sim %in% c('E', 'F') ~ "Sem correlação",
                         Sim %in% c('G', 'H') ~ "Com correlação",
                         Sim %in% c('I', 'J') ~ "Com correlação periódica") %>% 
           factor(levels = c("Sem correlação", "Com correlação", "Com correlação periódica")),
         real_value = case_when(name == "beta_0" ~ -1, 
                                name == "beta_1" ~ .5)) %>% 
  mutate(SE = (value - real_value)^2) %>% 
  group_by(Sim, Modelo, name, Tamanho, Casos) %>% 
  summarise(RMSE = sqrt(mean(SE))) %>% 
  ggplot(aes(x=Tamanho, y=RMSE, fill=Modelo))+
  geom_col(position="dodge") + 
  facet_grid(Casos~turn_greek(name), scales = 'free')

# RMSE - Rho ----
parameters  %>% 
  select(-(starts_with("beta"))) %>% 
  mutate_at(vars(matches("rho")), ~ifelse(is.na(.x) & Modelo =="RLA", rho_1, .x)) %>% 
  mutate_all(~tidyr::replace_na(.x, 0)) %>% 
  pivot_longer(-c(Sim, Modelo, rep)) %>% 
  mutate(Tamanho=case_when(Sim %in% c('E', 'G','I') ~ "n=280",
                           Sim %in% c('F', 'H','J') ~ "n=840"),
         Casos=case_when(Sim %in% c('E', 'F') ~ "Sem correlação",
                         Sim %in% c('G', 'H') ~ "Com correlação",
                         Sim %in% c('I', 'J') ~ "Com correlação periódica") %>% 
           factor(levels = c("Sem correlação", "Com correlação", "Com correlação periódica")),
         real_value = case_when(Casos == "Sem correlação" ~ 0,
                                Casos =="Com correlação"~ .5,
                                name=="rho_1" ~ 0.6,
                                name=="rho_2" ~ 0.3, 
                                name=="rho_3" ~ 0.1,
                                name=="rho_4" ~ -0.2, 
                                name=="rho_5" ~ -0.4,
                                name=="rho_6" ~ -0.2, 
                                name=="rho_7" ~ 0.2)) %>%
  mutate(SE = (value - real_value)^2) %>% 
  group_by(Sim, Modelo, name, Tamanho, Casos) %>% 
  summarise(RMSE = sqrt(mean(SE))) %>% 
  ggplot(aes(x=Tamanho, y=RMSE, fill=Modelo))+
  geom_col(position="dodge") + 
  facet_grid(Casos~turn_greek(name), scales = 'free')

# Especificity, Sensibility and Accuracy ----
metrics <- tibble()
for(i in LETTERS[5:10]){
  acc_RL <- readr::read_csv(paste0(folder, "Dados/Parciais/", i, "/parcial_acc.glm.csv"), col_names = 'acc')
  pfn_RL <- readr::read_csv(paste0(folder, "Dados/Parciais/", i, "/parcial_pfn.glm.csv"), col_names = 'pfn')
  pfp_RL <- readr::read_csv(paste0(folder, "Dados/Parciais/", i, "/parcial_pfp.glm.csv"), col_names = 'pfp')
  metrics_RL <- bind_cols(acc_RL, pfp_RL, pfn_RL) %>% mutate(Sim=i, Modelo="RL", rep=1:1000)
  
  acc_RLA <- readr::read_csv(paste0(folder, "Dados/Parciais/", i, "/parcial_acc.MC.csv"), col_names = 'acc')
  pfn_RLA <- readr::read_csv(paste0(folder, "Dados/Parciais/", i, "/parcial_pfn.MC.csv"), col_names = 'pfn')
  pfp_RLA <- readr::read_csv(paste0(folder, "Dados/Parciais/", i, "/parcial_pfp.MC.csv"), col_names = 'pfp')
  metrics_RLA <- bind_cols(acc_RLA, pfp_RLA, pfn_RLA) %>% mutate(Sim=i, Modelo="RLA", rep=1:1000)
  
  acc_RLAP <- readr::read_csv(paste0(folder, "Dados/Parciais/", i, "/parcial_acc.MCP.csv"), col_names = 'acc')
  pfn_RLAP <- readr::read_csv(paste0(folder, "Dados/Parciais/", i, "/parcial_pfn.MCP.csv"), col_names = 'pfn')
  pfp_RLAP <- readr::read_csv(paste0(folder, "Dados/Parciais/", i, "/parcial_pfp.MCP.csv"), col_names = 'pfp')
  metrics_RLAP <- bind_cols(acc_RLAP, pfp_RLAP, pfn_RLAP) %>% mutate(Sim=i, Modelo="RLAP", rep=1:1000)
  
  metrics <-  bind_rows(metrics, metrics_RL, metrics_RLA, metrics_RLAP)
}

metrics  %>% 
  mutate_all(~tidyr::replace_na(.x, 0)) %>% 
  pivot_longer(-c(Sim, Modelo, rep)) %>% 
  mutate(Tamanho=case_when(Sim %in% c('E', 'G','I') ~ "n=280",
                           Sim %in% c('F', 'H','J') ~ "n=840"),
         Casos=case_when(Sim %in% c('E', 'F') ~ "Sem correlação",
                         Sim %in% c('G', 'H') ~ "Com correlação",
                         Sim %in% c('I', 'J') ~ "Com correlação periódica") %>% 
           factor(levels = c("Sem correlação", "Com correlação", "Com correlação periódica"))) %>% 
  group_by(Sim, Modelo, Tamanho, Casos, name) %>% 
  summarise(value=mean(value)) %>% 
  ggplot(aes(x=Tamanho, y=value, fill=Modelo))+
  geom_col(position="dodge") + 
  facet_grid(Casos~turn_greek(name))
