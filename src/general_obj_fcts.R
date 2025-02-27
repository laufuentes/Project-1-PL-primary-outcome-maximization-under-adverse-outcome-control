library(tidyverse)
library(parallel)

# Policy definition function
sigma_beta <- function(psi, beta = 1 / 2, centered = FALSE) {
  c_beta <- 1 / log((1 + exp(beta)) / (1 + exp(-beta)))
  if (centered == TRUE) {
    cent <- 0.5 - c_beta * log((1 + exp(beta * 0)) / (1 + exp(-beta)))
  } else {
    cent <- 0
  }
  out <- c_beta * log((1 + exp(beta * psi)) / (1 + exp(-beta))) + cent
  return(out)
}

# Risk function for CATE
R_p0 <- function(psi, X, option) {
  out <- mean(psi^2 - 2 * psi * delta_Y(X, option))
  return(out)
}

# Constraint function
S_p0 <- function(psi, X, option, beta = 1 / 2, alpha = 0.1, centered = FALSE) {
  out <- mean(sigma_beta(psi, beta, centered) * delta_Z(X, option)) - alpha
  return(out)
}

# Define the objective function
L <- function(psi, comb) {
  beta <- comb$beta
  lambda <- comb$lambda
  term1 <- R_p0(psi, covariates, option) # Summing over all dimensions
  term2 <- S_p0(psi, covariates, option, beta, alpha, centered) # Summing over all dimensions
  out <- term1 + lambda * term2
  return(out) # Adjusting for n dimensions
}

policy_values <- function(
  psi,
  counterfactual_outcomes,
  beta,
  centered = FALSE,
  alpha = 0.1
) {
  sigma_psi <- sigma_beta(psi, beta = beta, centered)
  policy <- rbinom(n, 1, sigma_psi)
  y1 <- counterfactual_outcomes[1]
  y0 <- counterfactual_outcomes[2]
  out <- mean(policy * y1 + (1 - policy) * y0)
  return(out)
}

optimize_grid <- function(
  beta,
  lambda,
  covariates,
  option,
  alpha,
  centered,
  epsilon
) {
  optim_result <- optim(
    par = delta_Y(covariates, option),
    fn = function(psi)
      L(psi, covariates, lambda, option, beta, alpha, centered = centered),
    method = "L-BFGS-B",
    lower = -1,
    upper = 1
  )

  return(optimal_x = optim_result$par) # Return only optimal_x
}
