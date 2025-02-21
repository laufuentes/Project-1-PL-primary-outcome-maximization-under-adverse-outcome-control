set.seed(2025)
setwd("~/Documents/PhD/Project 1 - Policy learning - Constraints - Multiple outcome/simulations_new_approach_AC")
library(tidyverse)
library(dplyr)
library(kernlab)
library(lbfgs)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)


############################
##### Data generation #####
############################
n <- 5*1e3
option <-1

h_Y<- function(X,A, option){
  if(option==3){
    2*(1- X[,1]-X[,2])*A 
  }else{
    8*(1-X[,1]^2 -X[,2]*2)*A
  }
}

h_R<-function(X,A, option){
  if(option==3){
    (1+X[,1]-X[,2])*A
  }else{
    (X[,1]+X[,2]-0.3)*A
  }
}

data_gen <- function(n, option){
  if(option==1){
    u <- 10
    beta <- 3
    beta.0 <- 1
    
    theta_Y0 <- c(-0.01,0.2,-0.05, 1.5/4,2.5/4, 1.5/4,1.5/4 )
    theta_Y1<- c(0,0,0, 1, 1, 1, 1)
    
    theta_Z <- c(-0.009,0,-0.09,0,0,0,0)
    # basic characteristics 
    height <- rnorm(n, mean=163, sd=7)
    age<- as.integer(runif(n, min=18,max=60))
    bmi <- rnorm(n,21.75, sd=3)
    
    #diseases
    diabetes<- rbinom(n,1,p=ifelse(bmi>25, 0.2,4.7*1e-2))
    cancer <- rbinom(n,1, p=1e-2)
    hyperthrd <- rbinom(n,1, p=2*1e-2)
    hypothrd <- rbinom(n,1,p=ifelse(age>60,1*1e-2, 3*1e-2))
    
    X <-data.frame(height,age,bmi,diabetes, cancer, hyperthrd, hypothrd)
    Treatment <- rbinom(n,1,p=0.5) # Random treatment allocation
    
    epsilon <- rnorm(n,mean=0, sd = 1)
    
    Y_0 <- as.integer(pmax(u - (as.matrix(X)%*%theta_Y0)+rnorm(n, mean=0, sd=1),0))
    Y_1 <- as.integer(pmax(Y_0 + beta*(1 - as.matrix(X)%*%theta_Y1) + rnorm(n, mean=0, sd=1),0))
    Y.0 <- (Y_0-min(Y_0,Y_1))/(max(Y_0,Y_1)-min(Y_0,Y_1))
    Y.1 <- (Y_1-min(Y_0,Y_1))/(max(Y_0,Y_1)-min(Y_0,Y_1))

    logit_prob_Z.0 <- beta.0 + as.matrix(X)%*%theta_Z 
    logit_prob_Z.1 <- logit_prob_Z.0 + 2
    
    p0 <- 1/(1+exp(-logit_prob_Z.0))
    p1 <- 1/(1+exp(-logit_prob_Z.1))
    Z.1<- rbinom(n,1,p1)
    Z.0<- rbinom(n,1,p0)
  }else if(option==2){
    u <- 10
    beta <- 3
    beta.0 <- 1
    
    theta_Y0 <- c(-0.01,0.2,-0.05, 1.5/4,2.5/4, 1.5/4,1.5/4 )
    theta_Y1<- c(0,0,0, 1, 1, 1, 1)
    
    # basic characteristics
    height <- rnorm(n, mean=163, sd=7)
    age<- as.integer(runif(n, min=18,max=60))
    bmi <- rnorm(n,21.75, sd=3)
    
    #diseases
    diabetes<- rbinom(n,1,p=ifelse(bmi>25, 0.2,4.7*1e-2))
    cancer <- rbinom(n,1, p=1e-2)
    hyperthrd <- rbinom(n,1, p=2*1e-2)
    hypothrd <- rbinom(n,1,p=ifelse(age>60,1*1e-2, 3*1e-2))
    
    X <-data.frame(height,age,bmi,diabetes, cancer, hyperthrd, hypothrd)
    Treatment <- rbinom(n,1,p=0.5) # Random treatment allocation
    
    epsilon <- rnorm(n,mean=0, sd = 1)
    Y_0 <- as.integer(pmax(u - (as.matrix(X)%*%theta_Y0)+rnorm(n, mean=0, sd=1),0))
    Y_1 <- as.integer(pmax(Y_0 + beta*(1 - as.matrix(X)%*%theta_Y1) + rnorm(n, mean=0, sd=1),0))
    Y.0 <- (Y_0-min(Y_0,Y_1))/(max(Y_0,Y_1)-min(Y_0,Y_1))
    Y.1 <- (Y_1-min(Y_0,Y_1))/(max(Y_0,Y_1)-min(Y_0,Y_1))
    
    Y <- ifelse(Treatment==1,Y.1, Y.0)
    
    p0 <- ifelse(
      bmi < 20 & height < 150,
      runif(length(bmi), 0.4, 0.6), # Generate random values for true condition
      runif(length(bmi), 0.01, 0.1) # Generate random values for false condition
    )
    
    p1 <- ifelse(
      bmi < 25 & height < 155,
      p0+ runif(length(bmi),0.3,0.4),
      p0+ runif(length(bmi),0.2,0.3)
    )
    Z.1<- rbinom(n,1,p1)
    Z.0<- rbinom(n,1,p0)
  }else{
    Treatment <- rbinom(n,1,0.5)
    X <- matrix(runif(n*10,0,1),n,10)
    epsilon_Y <- rnorm(n,0,1)
    epsilon_R <- pmin(rnorm(n,0,1),1)
    
    Y_1 <- 1 - 2*X[,1] + X[,2] - X[,3] + h_Y(X,rep(1,n),option) + epsilon_Y
    Y_0 <- 1 - 2*X[,1] + X[,2] - X[,3] + h_Y(X,rep(-1,n),option) + epsilon_Y
    
    Y.1 <- (Y_1-min(Y_0,Y_1))/(max(Y_0,Y_1)-min(Y_0,Y_1))
    Y.0 <- (Y_0-min(Y_0,Y_1))/(max(Y_0,Y_1)-min(Y_0,Y_1))

    p1<- 2+ X[,1] + h_R(X,rep(1,n),option)+epsilon_R
    p0 <- 2+ X[,1] + h_R(X,rep(-1,n),option)+epsilon_R
  }
  df_complete <- data.frame(X,Treatment,y1=Y.1,y0=Y.0,p1=p1,p0=p0, CATE.sign_Y = as.factor(ifelse(Y.1-Y.0>0,1,-1)))
  df_obs<- data.frame(X,Treatment,Y=ifelse(Treatment==1,Y.1,Y.0),Z=ifelse(Treatment==1, p1, p0))
  return(list(df_complete, df_obs))
}

exp <- data_gen(n,option)

df_complete<- exp[[1]]
#### REGARDER ESP COND != CAUSALES
delta_Y <- (df_complete$y1-df_complete$y0)
delta_R <- df_complete$p1-df_complete$p0

if(option<3){
  plot_Y_sign<- ggplot(df_complete, aes(x=height, y=bmi, color=CATE.sign_Y))+
    geom_point(alpha = 0.5)
  p0_plot<- ggplot(df_complete, aes(x=height, y=bmi, color=(p0)))+
    geom_point(alpha = 0.5)+
    scale_color_gradient(low = "yellow", high = "red", limits = c(0, max(df_complete$p0)))
  p1_plot<- ggplot(df_complete, aes(x=height, y=bmi, color=(p1)))+
    geom_point(alpha = 0.5)+
    scale_color_gradient(low = "yellow", high = "red", limits = c(0, max(df_complete$p1)))
  
  p_plot<- ggplot(df_complete, aes(x=height, y=bmi, color=(delta_R)))+
    geom_point(alpha = 0.5)+
    scale_color_gradient(low = "blue", high = "green", limits = c(0.1, 0.6))
  
  if (option == 2) {
    p_plot <- p_plot +
      geom_rect(aes(xmin = min(height), xmax = 150, ymin = min(bmi), ymax = 20), 
                color = "black", fill = NA, linewidth = 1, inherit.aes = FALSE) +
      geom_rect(aes(xmin = min(height), xmax = 155, ymin = min(bmi), ymax = 25), 
                color = "black", fill = NA, linewidth = 1, inherit.aes = FALSE)
  }
  combined_plot <- grid.arrange(
    plot_Y_sign, p0_plot, p1_plot, p_plot,
    ncol = 1  
  )
}else{
  plot_Y_sign<- ggplot(df_complete, aes(x=X1, y=X2, color=CATE.sign_Y))+
    geom_point(alpha = 0.5)
  p0_plot<- ggplot(df_complete, aes(x=X1, y=X2, color=(p0)))+
    geom_point(alpha = 0.5)+
    scale_color_gradient(low = "yellow", high = "red", limits = c(min(df_complete$p0), max(df_complete$p0)))
  p1_plot<- ggplot(df_complete, aes(x=X1, y=X2, color=(p1)))+
    geom_point(alpha = 0.5)+
    scale_color_gradient(low = "yellow", high = "red", limits = c(min(df_complete$p1), max(df_complete$p1)))
  
  p_plot<- ggplot(df_complete, aes(x=X1, y=X2, color=(delta_R)))+
    geom_point(alpha = 0.5)+
    scale_color_gradient(low = "blue", high = "green")
  combined_plot <- grid.arrange(
    plot_Y_sign, p0_plot, p1_plot, p_plot,
    ncol = 1  
  )
}


ggsave(paste0("images/synthetic_setting_",option,".pdf"),combined_plot)

alpha<- 0.1
#0.5 
#0.1
#0.05
############################
##### Optimization pb #####
############################
sigma_beta <- function(psi, beta=1/2){
  c_beta <- 1/log((1+exp(beta)) / (1+exp(-beta)))
  out <- c_beta * log((1+exp(beta*psi))/(1+exp(-beta)))
  return(out)
}

pi_opt <- ifelse(delta_Y>0,1,0)
mean(pi_opt*df_complete$y1 + (1-pi_opt)*df_complete$y0)


R_p0 <- function(psi){
  mean(psi^2 -2*psi*delta_Y)
}

S_p0 <- function(psi){
  mean(sigma_beta(psi)*delta_R)-alpha
}

# Define the objective function
objective <- function(x, lambda) {
  term1 <- sum(x^2 - 2 * x * delta_Y)  # Summing over all dimensions
  term2 <- lambda * sum(sigma_beta(x) * delta_R)  # Summing over all dimensions
  return(term1 + term2 - lambda * alpha * n)  # Adjusting for n dimensions
}
# Define lambda values
lambda_values <- seq(0, 15, 0.05)  

# Optimization: Find the x that minimizes the objective for each lambda
library(optimx)

results <- data.frame(lambda = lambda_values, optimal_x = NA, r_lambda = NA, s_lambda=NA, obj=NA)

for (i in seq_along(lambda_values)) {
  lambda <- lambda_values[i]
  optim_result <- optim(par = delta_Y, fn = function(x) objective(x, lambda), method = "L-BFGS-B", lower=-1,upper=1)
  results$optimal_x[i] <- list(optim_result$par)  # Store optimal x as a list
  results$r_lambda[i] <- R_p0(optim_result$par)
  results$s_lambda[i] <- S_p0(optim_result$par)
  results$obj[i]<-R_p0(optim_result$par) +lambda*S_p0(optim_result$par)
}


# Plot lambda vs optimal objective value
library(ggplot2)
library(ggpubr)

lambda_evol<- ggplot(results, aes(x = lambda)) +
  geom_point(aes(y = r_lambda, color = "R_lambda")) +
  geom_point(aes(y = s_lambda, color = "S_lambda")) +
  geom_point(aes(y=r_lambda+lambda*s_lambda, color="U_lambda"))+
  labs(title = "R_lambda and S_lambda vs Lambda", x = "Lambda", y = "Value") +
  scale_color_manual(name = "Functions", values = c("R_lambda" = "blue", "S_lambda" = "red", "U_lambda"="gray")) +
  theme_minimal()

ggsave(paste0("images/lambda_evol_",option,".pdf"),lambda_evol)


gamma_plot_funct <- function(results, idx, option) {
  policy <- sigma_beta(results$optimal_x[[idx]])
  
  # Initialize base plot
  if(option<3){
    p <- ggplot(
      cbind(df_complete, treat_proba = policy),
      aes(x = height, y = bmi, color = treat_proba)
    ) +
      geom_point(alpha = 0.5) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                            midpoint = 0.5, limits = c(0, 1), 
                            oob = scales::squish) +
      labs(title = bquote(lambda == .(lambda_values[[idx]]))) +
      theme_minimal() +
      theme(legend.position = "right")
    
    # Add rectangles **only if option == 1**
    if (option == 2) {
      p <- p +
        geom_rect(aes(xmin = min(height), xmax = 150, ymin = min(bmi), ymax = 20), 
                  color = "black", fill = NA, linewidth = 1, inherit.aes = FALSE) +
        geom_rect(aes(xmin = min(height), xmax = 155, ymin = min(bmi), ymax = 25), 
                  color = "black", fill = NA, linewidth = 1, inherit.aes = FALSE)
    }
  }else{
    p <- ggplot(
      cbind(df_complete, treat_proba = policy),
      aes(x = X1, y = X2, color = treat_proba)
    ) +
      geom_point(alpha = 0.5) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                            midpoint = 0.5, limits = c(0, 1), 
                            oob = scales::squish) +
      labs(title = bquote(lambda == .(lambda_values[[idx]]))) +
      theme_minimal() +
      theme(legend.position = "right")
  }
  
  return(p)
}

lambda_discr <- as.integer(seq(1, length(lambda_values), length.out = 10))

plots <- lapply(lambda_discr, function(x) gamma_plot_funct(results, x, option))
plots_no_legend <- lapply(plots, function(p) p + theme(legend.position = "none"))

legend <- get_legend(plots[[1]])
combined_plots <- plot_grid(plotlist = plots_no_legend, ncol = 5, nrow = 2, align = "hv")
final_plot <- plot_grid(combined_plots, legend, ncol=1, rel_heights = c(5,2))
# Display the final plot
print(final_plot)

ggsave(paste0("images/geom_points_lambda_",option,".pdf"),final_plot)


idx_min <- which.min(results$s_lambda[results$s_lambda>0])
plot_none<-gamma_plot_funct(results, 1, option )
plot_min <- gamma_plot_funct(results, idx_min, option )
plot_max <-  gamma_plot_funct(results, idx_min+1, option)
opt_plots <- plot_grid(plot_none, plot_min,plot_max, ncol=3)
ggsave(paste0("images/geom_points_lambda_opt_",option,".pdf"),opt_plots)


###############################
library(ggplot2)
library(gganimate)

if(option<3){
  df_covariates <- df_complete %>% 
    select("bmi","height") %>% 
    mutate(id=1:n)
  
  df_policy <- as.data.frame(sapply(results$optimal_x, sigma_beta))
  colnames(df_policy) <- paste0("lambda_",lambda_values)  # Adjust step size if needed
  df_policy$id <- 1:n
  
  df_long <- df_policy %>%
    pivot_longer(cols = starts_with("lambda_"),
                 names_to = "lambda",
                 values_to = "policy") %>%
    mutate(lambda = as.numeric(sub("lambda_", "", lambda)))
  
  df_final <- df_long %>%
    left_join(df_covariates, by = "id")
  
  
  p <- ggplot(df_final, aes(x = height, y = bmi, color = policy)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                          midpoint = 0.5, limits = c(0, 1), 
                          oob = scales::squish)+
    labs(title = "Treatment Probability for λ = {closest_state}",
         x = "Height", y = "BMI", color = "Treatment Probability") +
    transition_states(lambda) +
    ease_aes('linear')
  
  if (option == 2) {
    p <- p +
      geom_rect(aes(xmin = min(height), xmax = 150, ymin = min(bmi), ymax = 20), 
                color = "black", fill = NA, linewidth = 1, inherit.aes = FALSE) +
      geom_rect(aes(xmin = min(height), xmax = 155, ymin = min(bmi), ymax = 25), 
                color = "black", fill = NA, linewidth = 1, inherit.aes = FALSE)
  }
}else{
  df_covariates <- df_complete %>% 
    select("X1","X2") %>% 
    mutate(id=1:n)
  
  df_policy <- as.data.frame(sapply(results$optimal_x, sigma_beta))
  colnames(df_policy) <- paste0("lambda_",lambda_values)  # Adjust step size if needed
  df_policy$id <- 1:n
  
  df_long <- df_policy %>%
    pivot_longer(cols = starts_with("lambda_"),
                 names_to = "lambda",
                 values_to = "policy") %>%
    mutate(lambda = as.numeric(sub("lambda_", "", lambda)))
  
  df_final <- df_long %>%
    left_join(df_covariates, by = "id")
  
  
  p <- ggplot(df_final, aes(x = X1, y = X2, color = policy)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                          midpoint = 0.5, limits = c(0, 1), 
                          oob = scales::squish)+
    labs(title = "Treatment Probability for λ = {closest_state}",
         x = "X1", y = "X2", color = "Treatment Probability") +
    transition_states(lambda) +
    ease_aes('linear')
  
}

animate(p, fps = 20, duration = 10, width = 800, height = 600)
anim_save(paste0("images/policy_animation_",option,".gif"))

