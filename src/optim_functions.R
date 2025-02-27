library(tidyverse)
library(optimx)
library(lbfgs)

source("src/general_obj_fcts.R")

optimize_combination <- function(idx, initial_guess, param_combinations) {
    params <- param_combinations[idx, ]
    optim_result <- optim(
        par = initial_guess,
        fn = L,
        comb = params,
        method = "L-BFGS-B",
        lower = rep(-1, n),
        upper = rep(1, n)
    )

    policy <- optim_result$par
    return(optim_result$par)
}

process_policy <- function(
    idx,
    param_combinations,
    policies,
    covariates,
    option,
    df_complete,
    centered,
    alpha
) {
    # Extract the policy for the current index
    policy <- policies[[idx]]
    results <- data.frame(
        lambda = param_combinations$lambda[idx],
        beta = param_combinations$beta[idx],
        optimal_x = I(list(policy)), # I() wraps the list to avoid issues with data frames
        risk = R_p0(policy, covariates, option),
        constraint = S_p0(
            policy,
            covariates,
            option,
            param_combinations$beta[idx],
            centered = centered
        ),
        obj = L(psi = policy, comb = param_combinations[idx, ]),
        policy_value = policy_values(
            psi = policy,
            counterfactual_outcomes = c(df_complete$y1, df_complete$y0),
            beta = param_combinations$beta[idx],
            centered = centered,
            alpha = alpha
        )
    )
    colnames(results) <- c(
        "lambda",
        "beta",
        "optimal_x",
        "risk",
        "constraint",
        "obj",
        "policy_value"
    )
    return(results) # Return the updated results for this index
}
