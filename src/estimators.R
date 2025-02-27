nuissance_params <- function(n.folds, df) {
}

dr_learner <- function(counterfacts_test, df_test, prop_score, n.folds = 5) {
  n <- nrow(df_test)
  tau_DR <- rep(0, n)
  indices <- split(seq(n), sort(seq(n) %% n.folds))
  for (idx in indices) {
    train <- df_test[-idx, ]
    test <- df_test[idx, ]
    pseudo_outcome <- ((train$W - prop_score[-idx]) /
      (prop_score[-idx] * (1 - prop_score[-idx]))) *
      (train$Y -
        ifelse(
          train$W == 1,
          counterfacts_test[-idx, ]$mu1,
          counterfacts_test[-idx, ]$mu0
        )) +
      (counterfacts_test[-idx, ]$mu1 - counterfacts_test[-idx, ]$mu0)

    fit_DR <- lm(
      DR ~ .,
      data = data.frame(train %>% dplyr::select(-Y, -W), DR = pseudo_outcome)
    )
    tau_DR[idx] <- predict(fit_DR, test)
  }
  return(tau_DR)
}

train_cond_mean <- function(train, test) {
  T.learner.1 <- regression_forest(
    X = train %>%
      dplyr::filter(A = 1) %>%
      dplyr::select(-A, -Y) %>%
      as.matrix(),
    Y = train %>%
      dplyr::filter(A == 1) %>%
      dplyr::select(Y) %>%
      as.matrix()
  )

  T.learner.0 <- regression_forest(
    X = train %>%
      dplyr::filter(A == 0) %>%
      dplyr::select(-A, -Y) %>%
      as.matrix(),
    Y = train %>%
      dplyr::filter(A == 0) %>%
      dplyr::select(Y) %>%
      as.matrix()
  )

  mu1 <- predict(
    T.learner.1,
    newdata = test %>%
      dplyr::select(-Y, -A)
  )$predictions %>%
    as.vector()
  mu0 <- predict(
    T.learner.0,
    newdata = test %>%
      dplyr::select(-Y, -A)
  )$predictions %>%
    as.vector()

  return(list(mu1, mu0))
}

compute_propensity <- function(train, test) {
  n_obs <- dim(train)[1]
  logit <- probability_forest(
    train %>% dplyr::select(-Y, -W),
    as.factor(train$W)
  )
  prop.score <- predict(
    logit,
    newdata = test %>% dplyr::select(-Y, -W),
    type = "response"
  )$pred[, 2]
  return(prop.score)
}
