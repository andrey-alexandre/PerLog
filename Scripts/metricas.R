library(dplyr)
library(tidyr)

for(i in LETTERS[5:10]){
  load(paste0('/home/andrey/Projetos/PerLog/Dados/Metricas/', i, '/resultado.rds', sep = ''))
  
  pf_mat %>% 
    round(3) %>% as.data.frame %>% tibble::rownames_to_column() %>% 
    separate(col = 'rowname', into = c('Metrica', 'Amostra'), ' - ')%>% 
    arrange(Amostra) %>% 
    select(Amostra, Metrica, contains('M')) %>% 
    write.table(file = '/home/andrey/Projetos/PerLog/Dados/Metricas/Geral/pf_mat2.txt', 
                sep = '\t & \t', append = T, row.names = F, col.names = F, quote = F, , na = '-')
  
  mat.EP %>% cbind(mat_Bias_RMSE) %>% t() %>% 
    round(3) %>% as.data.frame %>% tibble::rownames_to_column() %>% 
    separate(col = 'rowname', into = c('Metrica', 'Modelo'), ' - ') %>% 
    arrange(Modelo) %>% 
    select(Modelo, Metrica, starts_with('V')) %>% 
    write.table(file = '/home/andrey/Projetos/PerLog/Dados/Metricas/Geral/mat_EP2.txt', 
                sep = '\t & \t', append = T, row.names = F, col.names = F, quote = F, na = '-')
}
