set.seed(2025)
setwd(
  "~/Documents/PhD/Project 1 - Policy learning - Constraints - Multiple outcome/Project1_PL_risk_benefit"
)

library(tidyverse)
library(dplyr)
library(kernlab)
library(lbfgs)

library(tidyverse)
library(dplyr)
library(optimx)

source("src/plot_fcts.R")
source("src/tool_box.R")
source("src/synthetic_data.R")
source("src/general_obj_fcts.R")

############################
##### Data generation #####
############################

#### General parameters ####

# Number of individuals
n <- 1e4
# Constraint-tolerance parameter
alpha <- 0.1 #0.5 #0.1 #0.05
epsilon <- 0.03
centered <- TRUE

# Setting
setting <- "Other_1"
#"IVF_1" #"IVF_2" #"Other_2" # "Other_1"

#### Generate synthetic data ####
option <- option_det(setting, "_")
exp <- data_gen(n, option)
df_complete <- exp[[1]]

###
#true_cate_Y <- (df_complete$y1-df_complete$y0)
#true_cate_Z <- df_complete$p1-df_complete$p0

### Plot synthetic_setting scenario
synthetic_setting_plot(df_complete, option)


############################
#### Optimization part ####
############################
# Grid search candidates
## Beta
beta_values <- seq(0.05, 2, 0.05)
## Lambda
lambda_values <- seq(0, 15, 0.1)


#pi_opt <- ifelse(true_cate_Y>0,1,0)
#mean(pi_opt*df_complete$y1 + (1-pi_opt)*df_complete$y0)

# Optimization: Find the x that minimizes the objective for each lambda

results <- expand.grid(
  lambda = lambda_values,
  beta = beta_values,
  optimal_x = NA,
  risk = NA,
  constraint = NA,
  obj = NA,
  policy_value = NA
)
covariates <- df_complete %>% select(starts_with("X"))
results$optimal_x <- vector("list", nrow(results))

for (beta in beta_values) {
  for (lambda in lambda_values) {
    # Get the current row index for this (lambda, beta) pair
    idx <- which(results$lambda == lambda & results$beta == beta)

    optim_result <- optim(
      par = delta_Y(covariates, option),
      fn = function(psi)
        L(psi, covariates, lambda, option, beta, alpha, centered = centered),
      method = "L-BFGS-B",
      lower = -1,
      upper = 1
    )

    # Store results
    results$optimal_x[[idx]] <- optim_result$par
    results$risk[idx] <- R_p0(optim_result$par, covariates, option)
    results$constraint[idx] <- S_p0(
      optim_result$par,
      covariates,
      option,
      beta,
      centered = centered
    )
    results$obj[idx] <- L(
      optim_result$par,
      covariates,
      lambda,
      option,
      beta,
      alpha,
      centered = centered
    )
    results$policy_value[idx] <- policy_values(
      optim_result$par,
      c(df_complete$y1, df_complete$y0),
      beta
    )

    # Stop increasing lambda if constraint is met
    if (results$constraint[idx] < -epsilon) {
      break # Exit the loop for this beta, move to the next beta
    }
  }
}

results <- na.omit(results)

write.csv(
  results %>% select(-optimal_x),
  paste0("opt_results/", setting, ".csv")
)


################################
#### Plot results ####
################################

idx_opt_pol <- which(
  results$policy_value == max(results$policy_value[results$constraint <= 0])
)
idx_opt_obj <- which(results$obj == max(results$obj[results$constraint <= 0]))


idx_opt <- idx_opt_pol[1]

lambda_evol(
  results %>% select(-optimal_x),
  option,
  results$lambda[idx_opt],
  results$beta[idx_opt]
)
geom_points_fct(results, idx_opt, df_complete, option)

gamma_lambda_plot(results, option, df_complete)

animation_plot(
  results[which(results$beta == results$beta[idx_opt]), ],
  results$lambda[which(results$beta == results$beta[idx_opt])],
  option,
  df_complete
)
