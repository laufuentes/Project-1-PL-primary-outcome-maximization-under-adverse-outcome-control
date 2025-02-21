setwd("~/Documents/PhD/Project 1 - Policy learning - Constraints - Multiple outcome/simulations_new_approach_AC")
set.seed(2025)

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
    Y.0 <- as.integer(pmax(u - (as.matrix(X)%*%theta_Y0)+rnorm(n, mean=0, sd=1),0))
    Y.1 <- as.integer(pmax(Y.0 + beta*(1 - as.matrix(X)%*%theta_Y1) + rnorm(n, mean=0, sd=1),0))
    
    logit_prob_Z.0 <- beta.0 + as.matrix(X)%*%theta_Z 
    logit_prob_Z.1 <- logit_prob_Z.0 + 2
    
    prob.0 <- 1/(1+exp(-logit_prob_Z.0))
    prob.1 <- 1/(1+exp(-logit_prob_Z.1))
    Z.1<- rbinom(n,1,prob.1)
    Z.0<- rbinom(n,1,prob.0)
  }else{
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
    Y.0 <- as.integer(pmax(u - (as.matrix(X)%*%theta_Y0)+rnorm(n, mean=0, sd=1),0))
    Y.1 <- as.integer(pmax(Y.0 + beta*(1 - as.matrix(X)%*%theta_Y1) + rnorm(n, mean=0, sd=1),0))
    
    Y <- ifelse(Treatment==1,Y.1, Y.0)
    
    prob.0 <- ifelse(
      bmi < 20 & height < 150,
      runif(length(bmi), 0.4, 0.6), # Generate random values for true condition
      runif(length(bmi), 0.01, 0.1) # Generate random values for false condition
    )
    
    prob.1 <- ifelse(
      bmi < 25 & height < 155,
      prob.0+ runif(length(bmi),0.3,0.4),
      prob.0+ runif(length(bmi),0.2,0.3)
    )
    Z.1<- rbinom(n,1,prob.1)
    Z.0<- rbinom(n,1,prob.0)
  }
  df_complete <- data.frame(X,Treatment,y1=Y.1,y0=Y.0,p1=prob.1,p0=prob.0, CATE.sign_Y = as.factor(ifelse(Y.1-Y.0>0,1,-1)))
  df_obs<- data.frame(X,Treatment,Y=ifelse(Treatment==1,Y.1,Y.0),Z=ifelse(Treatment==1, Z.1, Z.0))
  return(list(df_complete, df_obs))
}

n <- 1e6
dfs <- data_gen(n, option=2)
df_complete <- dfs[[1]]
df_obs <- dfs[[2]]

alpha <- 0.1
beta <- 1

sigma_beta <- function(psi){
  c_beta <- log((1+exp(beta)) / (1+exp(-beta)))
  return(c_beta * log((1+exp(beta*psi))/(1+exp(-beta))))
}

psi<- ((df_complete$y1 - min(df_complete$y1))/(max(df_complete$y1)-min(df_complete$y1)))-(df_complete$y0 -min(df_complete$y0))/(max(df_complete$y0) - min(df_complete$y0))

pi_X_etoile <- sapply(sigma_beta(psi),function(x){rbinom(1,1,p=x)})

pi_opt <- ifelse(df_complete$y1>df_complete$y0,1,0)

pi_opt_Z <- ifelse(df_complete$p1-df_complete$p0<alpha,1,0)


mean(pi_opt*(df_complete$y1-df_complete$y0))
mean(pi_opt*(df_complete$p1-df_complete$p0))-alpha

