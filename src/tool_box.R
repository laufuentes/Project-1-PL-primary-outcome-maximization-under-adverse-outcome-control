set.seed(2025)

option_det <- function(string, split_char) {
  res <- strsplit(string, split = split_char)
  return(res[[1]])
}
