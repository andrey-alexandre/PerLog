z_confidence_interval <- function(x){
  ux <- mean(x); sdx <- sd(x); n <- length(x)
  norm <- qnorm(p = .025, mean = ux, sd = sdx)
  bound <- norm * sdx * sqrt(n)
  
  output <- list(lower = ux - bound, upper = ux + bound)
  return(output)
}
z_test <- function(x){
  bound <- z_confidence_interval(x)
  
  y <- ifelse(between(x, bound$lower, bound$upper), 0, 1)
  return(y)
}

normal_bound <- function(x){
  ux <- mean(x); sdx <- sd(x); n <- length(x)
  norm <- qnorm(p = .025)
  
  output <- list(lower = norm, upper = -norm)
  return(output)
}
bound_test <- function(x){
  bound <- normal_bound(x)
  
  y <- ifelse(between(x, bound$lower, bound$upper), 0, 1)
  return(y)
}

std_glm <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Scripts/Resultados/E/parcial_std_glm.csv'), stringsAsFactors = F, header = F)
std_MC <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Scripts/Resultados/E/parcial_std_MC.csv'), stringsAsFactors = F, header = F)
std_MCP <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Scripts/Resultados/E/parcial_std_MCP.csv'), stringsAsFactors = F, header = F)

std_glm
std_glm %>% apply(2, z_confidence_interval)
std_glm %>% apply(2, z_test) %>% colMeans
std_glm %>% apply(2, normal_bound)
std_glm %>% apply(2, bound_test) %>% colMeans

std_MC
std_MC %>% apply(2, z_confidence_interval)
std_MC %>% apply(2, z_test) %>% colMeans
std_MC %>% apply(2, normal_bound)
std_MC %>% apply(2, bound_test) %>% colMeans

std_MCP
std_MCP %>% apply(2, z_confidence_interval)
std_MCP %>% apply(2, z_test) %>% colMeans
std_MCP %>% apply(2, normal_bound)
std_MCP %>% apply(2, bound_test) %>% colMeans
      