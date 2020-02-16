z <- function(x){
  ux <- mean(x); sdx <- sd(x); n <- length(x)
  norm <- qnorm(p = .025, mean = ux, sd = sdx)
  bound <- norm * sdx * sqrt(n)
  
  output <- list(lower = ux - bound, upper = ux + bound)
  return(output)
}

test <- function(x){
  bound <- z(x)
  
  y <- ifelse(between(x, bound$lower, bound$upper), 0, 1)
  return(y)
}

std_glm <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Scripts/Resultados/', arg[length(arg)], '/parcial_std_glm.csv'), stringsAsFactors = F, header = F)
std_MC <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Scripts/Resultados/', arg[length(arg)], '/parcial_std_MC.csv'), stringsAsFactors = F, header = F)
std_MCP <- read.csv(file = paste0('/home/andrey/Projetos/PerLog/Scripts/Resultados/', arg[length(arg)], '/parcial_std_MCP.csv'), stringsAsFactors = F, header = F)

std_glm
std_glm %>% apply(2, z)
std_glm %>% apply(2, test) %>% colMeans

std_MC
std_MC %>% apply(2, z)
std_MC %>% apply(2, test) %>% colMeans

std_MCP
std_MCP %>% apply(2, z)
std_MCP %>% apply(2, test) %>% colMeans
