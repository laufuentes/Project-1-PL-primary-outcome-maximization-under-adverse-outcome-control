setwd("~/Documents/PhD/Project 1 - Policy learning - Constraints - Multiple outcome/Project1_PL_risk_benefit")
library(tidyverse)

# Policy definition function
sigma_beta <- function(psi, beta=1/2, centered=FALSE){
  c_beta <- 1/log((1+exp(beta)) / (1+exp(-beta)))
  if (centered==TRUE){
    cent <- 0.5 - c_beta * log((1+exp(beta*0))/(1+exp(-beta)))
  }else{
    cent <- 0
  }
  out <- c_beta * log((1+exp(beta*psi))/(1+exp(-beta))) + cent
  return(out)
}

# Risk function for CATE
R_p0 <- function(psi,X,option){
  out <- mean(psi^2 -2*psi*delta_Y(X,option))
  return(out)
}

# Constraint function
S_p0 <- function(psi,X,option, beta=1/2,alpha=0.1, centered=FALSE){
  out<- mean(sigma_beta(psi, beta, centered)*delta_Z(X,option))-alpha
  return(out)
}

# Define the objective function
L <- function(psi,X, lambda,option,beta,alpha, centered=FALSE) {
  term1 <- R_p0(psi,X,option)  # Summing over all dimensions
  term2 <- S_p0(psi,X,option,beta,alpha,centered)  # Summing over all dimensions
  out <- term1 + lambda*term2
  return(out)  # Adjusting for n dimensions
}

policy_values <- function(psi, counterfactual_outcomes, beta, centered=FALSE, alpha=0.1){
  sigma_psi <- sigma_beta(psi, beta = beta, centered)
  policy <- rbinom(n,1,sigma_psi)
  y1<- counterfactual_outcomes[1]
  y0 <- counterfactual_outcomes[2]
  out <- mean(policy*y1 + (1-policy)*y0)
  return(out)
}
