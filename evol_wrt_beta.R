library(tidyverse)
sigma_beta <- function(psi, beta=1/2){
  c_beta <- 1/log((1+exp(beta)) / (1+exp(-beta)))
  out <- c_beta * log((1+exp(beta*psi))/(1+exp(-beta)))
  return(out)
}

tibble(
  psi=seq(-1,1,length.out=10)
) %>% 
  mutate(
    sigma_beta_1=sigma_beta(psi,beta=1),
    sigma_beta_2=sigma_beta(psi,beta=2),
    sigma_beta_0.5=sigma_beta(psi,beta=0.5)
) %>% 
  pivot_longer(
    cols=!psi,
    names_to="beta",
    names_prefix="sigma_beta_",
    names_transform = as.factor,
    values_to="sigma_beta"
  ) %>% 
  ggplot()+
  geom_line(aes(x=psi,y=sigma_beta,color=beta))
