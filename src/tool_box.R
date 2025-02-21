set.seed(2025)
setwd("~/Documents/PhD/Project 1 - Policy learning - Constraints - Multiple outcome/simulations_new_approach_AC")


option_det <- function(string, split_char){
  res<- strsplit(string, split = split_char)
  return(res[[1]])
}