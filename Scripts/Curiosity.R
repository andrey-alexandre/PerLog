library(dplyr)
library(tidyr)
library(ggplot2)

z_test_p <- function(x, y){
  H_0 <- !between(x, quantile(y, probs = .025), quantile(y, probs = .975))
  
  return(H_0)
}

z_test <- function(x){
  ux <- mean(x); sx <- sd(x); n <- length(x)
  upper <- ux + qnorm(.975) * sx/sqrt(n)
  lower <- ux - qnorm(.975) * sx/sqrt(n)
  
  
  H_0 <- !between(x, lower, upper)
  
  return(H_0)
}
sim <- c(LETTERS, paste0('A', LETTERS[1:5]))

Z <- data.frame(beta = seq(0, 1, along.with = sim), glm = 0, MC = 0, MCP = 0)
for(i in seq_along(sim)){
  std_glm <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/Fraca_', sim[i] ,'/parcial_std_glm.csv'), stringsAsFactors = F, header = F)[,2]
  std_MC <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/Fraca_', sim[i] ,'/parcial_std_MC.csv'), stringsAsFactors = F, header = F)[,2]
  std_MCP <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/Fraca_', sim[i] ,'/parcial_std_MCP.csv'), stringsAsFactors = F, header = F)[,2]
  est_glm <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/Fraca_', sim[i] ,'/parcial_est_glm.csv'), stringsAsFactors = F, header = F)[,2]
  est_MC <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/Fraca_', sim[i] ,'/parcial_est_MC.csv'), stringsAsFactors = F, header = F)[,2]
  est_MCP <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Dados/Parciais/Fraca_', sim[i] ,'/parcial_est_MCP.csv'), stringsAsFactors = F, header = F)[,2]
  
  if(sim[i] == 'E'){
    z_glm_0 <- est_glm/std_glm
    Z[i,2] <- mean(abs(z_glm_0) > qnorm(.975)) 
  }else{
    z_glm <- est_glm/std_glm
    Z[i,2] <- mean(z_test_p(z_glm, z_glm_0)) 
  }

  if(sim[i] == 'E'){
    z_MC_0 <- est_MC/std_MC
    Z[i,3] <- mean(abs(z_MC_0) > qnorm(.975))
  }else{
    z_MC <- est_MC/std_MC
    Z[i,3] <- mean(z_test_p(z_MC, z_MC_0)) 
  }
  
  if(sim[i] == 'E'){
    z_MCP_0 <- est_MCP/std_MCP
    Z[i,4] <- mean(abs(z_MCP_0) > qnorm(.975))
  }else{
    z_MCP <- est_MCP/std_MCP
    Z[i,4] <- mean(z_test_p(z_MCP, z_MCP_0)) 
  }
}

h0_rejection <- c(paste('bold("GLM"):', round(Z[1,2], 3)), paste('bold("MC"):', round(Z[1,3], 3)), paste('bold("MCP"):', round(Z[1,4], 3)))

Z %>% 
  gather(key = 'Modelo', value = 'TestPower', -beta) %>% 
  filter(beta != 0) %>% 
  ggplot(aes(x = beta, y = TestPower, col = Modelo)) +
  geom_line() +
  geom_point() +
  annotate("text", x = .1, y = .95, label = h0_rejection[1], parse = T) +
  annotate("text", x = .1, y = .90, label = h0_rejection[2], parse = T) +
  annotate("text", x = .1, y = .85, label = h0_rejection[3], parse = T)