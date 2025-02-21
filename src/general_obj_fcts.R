setwd("~/Documents/PhD/Project 1 - Policy learning - Constraints - Multiple outcome/simulations_new_approach_AC")
library(tidyverse)

# Policy definition function
sigma_beta <- function(psi, beta=1/2){
  c_beta <- 1/log((1+exp(beta)) / (1+exp(-beta)))
  out <- c_beta * log((1+exp(beta*psi))/(1+exp(-beta)))
  return(out)
}

# Risk function for CATE
R_p0 <- function(psi,X,option){
  out <- mean(psi^2 -2*psi*delta_Y(X,option))
  return(out)
}

# Constraint function
S_p0 <- function(psi,X,option, beta=1/2,alpha=0.1){
  out<- mean(sigma_beta(psi, beta)*delta_Z(X,option))-alpha
  return(out)
}

# Define the objective function
L <- function(psi,X, lambda,option,beta,alpha) {
  term1 <- R_p0(psi,X,option)  # Summing over all dimensions
  term2 <- S_p0(psi,X,option,beta,alpha)  # Summing over all dimensions
  out <- term1 + lambda*term2
  return(out)  # Adjusting for n dimensions
}

policy_values <- function(psi, counterfactual_outcomes, beta, alpha=0.1){
  sigma_psi <- sigma_beta(psi, beta = beta)
  policy <- rbinom(n,1,sigma_psi)
  y1<- counterfactual_outcomes[1]
  y0 <- counterfactual_outcomes[2]
  out <- mean(policy*y1 + (1-policy)*y0)
  return(out)
}
