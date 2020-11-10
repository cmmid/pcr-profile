pos_test <- function(x, qu){
 if(x < 0){
  return(0)
 }else{
  return(p_tab[diff == x, quantile(value, prob = qu)])
 }
}


detpr <- function(freqX, detect_within = 5, delay_to_result = 2, qu = 0.5){
 
 pos_sampling_times <- list()
 
 if(freqX >= detect_within){
  for(i in 0:(detect_within))
   pos_sampling_times[[i + 1]] <- i
 }else if(freqX == 1){
  pos_sampling_times[[1]] <- 0:(detect_within)
 }else{
  for(i in 0:(freqX - 1)) {
   pos_sampling_times[[i + 1]] <- seq(i, detect_within, by = freqX)
  }
 }
 # pos_sampling_times
 pj <- 0
 for(j in 1:length(pos_sampling_times)) {
  pk <- 0
  for(k in 1:length(pos_sampling_times[[j]])) {
   pk <- pos_test(pos_sampling_times[[j]][k] - delay_to_result, qu = qu) * (1 - inc_cumulative(pos_sampling_times[[j]][k])) *
    ifelse(k > 1, prod(vapply(pos_sampling_times[[j]][1:(k - 1)]  - delay_to_result, FUN = function(x){(1 - pos_test(x, qu = qu)) * (1 - inc_cumulative(x))}, FUN.VALUE = 1)), 1)
   pj <- pj + pk
  }
 }
 return(pj / freqX)
}


# Probability of testing positive
detpr2 <- function(freqX, detect_within = 5, delay_to_result = 2, qu = 0.5){
 
 pos_sampling_times <- list()
 
 if(freqX >= detect_within){
  for(i in 0:(detect_within - delay_to_result))
   pos_sampling_times[[i + 1]] <- i
 }else if(freqX == 1){
  pos_sampling_times[[1]] <- 0:(detect_within - delay_to_result)
 }else{
  for(i in 0:(freqX - 1)) {
   pos_sampling_times[[i + 1]] <- seq(i, detect_within - delay_to_result, by = freqX)
  }
 }
 # pos_sampling_times
 pj <- 0
 for(j in 1:length(pos_sampling_times)) {
  pk <- 0
  for(k in 1:length(pos_sampling_times[[j]])) {
   pk <- pos_test(pos_sampling_times[[j]][k], qu = qu) * 
    ifelse(k > 1, prod(vapply(pos_sampling_times[[j]][1:(k - 1)], FUN = function(x){(1 - pos_test(x, qu = qu))}, FUN.VALUE = 1)), 1)
   pj <- pj + pk
  }
 }
 return(pj / freqX)
} 