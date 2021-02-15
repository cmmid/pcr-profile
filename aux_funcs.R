
# This is the function that gets run over the table of different parameters for delay and testing frequency
# It runs the function "detpr" below on every sampled trajectory in p_tab and
# returns the median + intervals of estimated probability
prob_det_before_symp <- function(freqX, delay) {
  print(paste0("Frequency: ", freqX))
  print(paste0("Delay:", delay))
  
  out <- p_tab[, detpr(freqX = freqX, ptab = .SD, detect_within = 30, delay_to_result = delay), by = "iter"
  ][, .(median = median(V1), 
        bottom = quantile(V1, 0.025),
        top = quantile(V1, 0.975))]
  return(out)
}


## This function calculates the probability of detecting an infection 
# prior to symptom onset for a given delay and testing frequency
detpr <- function(freqX, detect_within = 5, delay_to_result = 2, ptab){
 
 # First step is to generate all the different testing regimes that
 # could occur for a given testing frequency
 pos_sampling_times <- list()
 
 if(freqX >= detect_within){
  for(i in 0:(detect_within))
   # If testing frequency > maximum number of days then
   # can only be tested once in the window
   pos_sampling_times[[i + 1]] <- i
 }else if(freqX == 1){
   # If the testing frequency is every day then
   # there is one regime
  pos_sampling_times[[1]] <- 0:(detect_within)
 }else{
  for(i in 0:(freqX - 1)) {
    # Generating different testing regimes depending on timing
    # of first test after infection
   pos_sampling_times[[i + 1]] <- seq(i, detect_within, by = freqX)
  }
 }

 # Second step is to calculate the probability of detection prior
 # to symptom onset for each testing regime.
 pj <- 0
 for(j in 1:length(pos_sampling_times)) {
  pk <- 0
  for(k in 1:length(pos_sampling_times[[j]])) {
   pk <- pos_test(pos_sampling_times[[j]][k] - delay_to_result, ptab = ptab) * (1 - inc_cumulative(pos_sampling_times[[j]][k])) *
    ifelse(k > 1, prod(vapply(pos_sampling_times[[j]][1:(k - 1)]  - delay_to_result, FUN = function(x){(1 - pos_test(x, ptab = ptab)) * (1 - inc_cumulative(x))}, FUN.VALUE = 1)), 1)
   pj <- pj + pk
  }
 }
 return(pj / freqX)
}

# A way of safely accessing the probability of detecting
# an infection from the posterior samples
pos_test <- function(x, ptab){
  if(x < 0){
    return(0)
  }else{
    return(ptab[diff == x, value])
  }
}

# This is the function that gets run over the table of different parameters for delay and testing frequency
# It runs the function "detpr2" below on every sampled trajectory in p_tab and
# returns the median + intervals of estimated probability
prob_det_before_X_asymp <- function(freqX, delay, within) {
  print(paste0("Frequency: ", freqX))
  print(paste0("Delay:", delay))
  print(paste0("Within:", within))
  
  out <- p_tab[, detpr2(freqX = freqX, ptab = .SD, detect_within = within, delay_to_result = delay), by = "iter"
  ][, .(median = median(V1), 
        bottom = quantile(V1, 0.025),
        top = quantile(V1, 0.975))]
  return(out)
}

## This function calculates the probability of detecting an asymptomatic infection 
# within a certain number of days for a given delay and testing frequency
detpr2 <- function(freqX, detect_within = 5, delay_to_result = 2, ptab){
 
 pos_sampling_times <- list()
 # Again we find all the different possible testing regimes for 
 # a given testing frequency
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
 # Slightly simpler probabilities of detection because there
 # is no need to worry about symptom onset occurring beforehand
 pj <- 0
 for(j in 1:length(pos_sampling_times)) {
  pk <- 0
  for(k in 1:length(pos_sampling_times[[j]])) {
   pk <- pos_test(pos_sampling_times[[j]][k], ptab = ptab) * 
    ifelse(k > 1, prod(vapply(pos_sampling_times[[j]][1:(k - 1)], FUN = function(x){(1 - pos_test(x, ptab = ptab))}, FUN.VALUE = 1)), 1)
   pj <- pj + pk
  }
 }
 return(pj / freqX)
} 

# Run the main analysis for a different ct value
fit_different_ct <- function(ct_threshold) {
  
  # Create new synthetic PCR results with different ct threshold
  test_lft <- test_final[, pcr_result := ifelse(ct == 0 | is.na(ct), 0, ifelse(ct <= ct_threshold, TRUE, FALSE))]
  
  # Update data
  dat_lft <- dat
  
  dat_lft$test_result <- test_lft$pcr_result %>% as.numeric()
  dat_lft$te_upper_bound <- test_lft[, te_upper_bound := ifelse(
    any(day[pcr_result == TRUE] < first_symp_day[pcr_result == TRUE]),
    min(day[pcr_result == TRUE & day < first_symp_day]),
    first_symp_day), by = num_id
  ][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]
  
  fit_lft <- rstan::sampling(mod, chains = 4, 
                             iter = 5000,
                             warmup = 2000,
                             data = dat_lft,
                             seed = seedx,
                             control = list(adapt_delta = 0.9, 
                                            stepsize = 0.75,
                                            max_treedepth = 13))
  res_lft <- rstan::extract(fit_lft)
  
  # Samples from LFT positive curve at different times since infection
  p_tab_lft <- as.data.table(res_lft$p)
  p_tab_lft <- data.table::melt(p_tab_lft)
  p_tab_lft$diff <- rep(p_vals, rep(12000, length(p_vals)))
  p_tab_lft$iter <- rep(1:12000, length(p_vals))
  p_tab_lft[, variable := NULL]
  
  # Figure S3A
  figS3a <- figure3a(ct_plot_dt = ct_plot_dt, ct_threshold = ct_threshold)
  
  # Figure S3B
  
  figS3b <- figure3b(res = res_lft, test_final = test_lft, ribbon_col = "firebrick2")
  
  # Figure S3C
  # Evaluated testing frequencies
  day_list <- c(2, 3, 4)
  
  tab <- data.table(every = rep(rep(day_list, rep(1, length(day_list))), 1),
                    within = 30,
                    delay = 0)
  
  tab <- tab[, prob_det_before_symp(freqX = every, delay = delay), by = c("every", "delay")]
  
  figS3c <- figureS3cd(tab, symp = TRUE)
  
  # Figure S3D
  tab2 <- data.table(every = rep(rep(day_list, rep(length(1), length(day_list))), 1),
                     within = 7,
                     delay = 0)
  
  tab2 <- tab2[, prob_det_before_X_asymp(freqX = every, delay = delay, within = within), by = c("every", "delay")]

  figS3d <- figureS3cd(tab2, symp = FALSE)
  

  bot_panel <- (figS3c + figS3d) + patchwork::plot_layout(guides = "collect")
  figureS3 <- (figS3a + figS3b) / bot_panel + plot_annotation(tag_levels = "A") 
  
  return(list("plot" = figureS3, "p_tab" = p_tab_lft, "fit" = fit_lft))
}
