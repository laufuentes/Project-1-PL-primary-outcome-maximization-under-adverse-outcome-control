library(parallel)
library(foreach)
library(snow)
library(tidyverse)

source("src/optim_functions.R")
source("src/tool_box.R")
source("src/general_obj_fcts.R")
source("src/synthetic_data.R")
source("src/plot_fcts.R")

# Create dummy data
n <- 1e4
# Setting
setting <- "Other_1"
#"IVF_1" #"IVF_2" #"Other_2" # "Other_1"

#### Generate synthetic data ####
option <- option_det(setting, "_")
exp <- data_gen(n, option)
df_complete <- exp[[1]]


alpha <- 0.1
centered <- TRUE
epsilon <- 0.03


beta_values <- seq(0.05, 2, 0.05)
## Lambda
lambda_values <- seq(0, 15, 0.1)

# Create a dummy results data frame
covariates <- as.data.frame(df_complete %>% select(starts_with("X"))) # Ensure covariates is a data frame
initial_guess <- delta_Y(covariates, option)
param_combinations <- expand.grid(lambda = lambda_values, beta = beta_values)


policies <- mclapply(
  1:nrow(param_combinations),
  function(i) {
    if (i %% 100 == 0) {
      print(i)
    }
    optimize_combination(i, initial_guess, param_combinations)
  },
  mc.cores = detectCores(),
  mc.preschedule = FALSE
)

res <- mclapply(
  1:nrow(param_combinations),
  function(i) {
    process_policy(
      i,
      param_combinations,
      policies,
      covariates,
      option,
      df_complete,
      centered,
      alpha
    )
  },
  mc.cores = detectCores(),
  mc.preschedule = FALSE
)

results <- as.data.frame(do.call(rbind, res))

# Save the results to a CSV file
write.csv(
  results %>% select(-optimal_x),
  paste0("opt_results/", setting, ".csv")
)

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

# animation_plot(
#   results[which(results$beta == results$beta[idx_opt]), ],
#   results$lambda[which(results$beta == results$beta[idx_opt])],
#   option,
#   df_complete
# )
