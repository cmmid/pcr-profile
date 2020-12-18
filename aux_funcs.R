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
  
  # LFT sensitivity plot
  p_lft <- data.frame(top = apply(res_lft$p, 2, quantile, prob = 0.975), 
                      bottom = apply(res_lft$p, 2, quantile, prob = 0.025),
                      y = apply(res_lft$p, 2, median),
                      days = seq(0, 30, 0.1)) %>%
    ggplot2::ggplot(ggplot2::aes(x = days, y=y,  ymin = bottom, ymax = top, fill = "Posterior distribution")) +
    ggplot2::geom_ribbon(alpha = 0.75) + 
    cowplot::theme_cowplot() + 
    ggplot2::geom_line(aes(lty = "Posterior median")) +
    ggplot2::labs(y = "Probability of positive LFT (%)", x = "Days since infection") +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), labels = paste0(seq(0, 100, 20))) +
    ggplot2::scale_x_continuous(breaks = c(0, seq(5, 30, 5)))
  
  # Generate empirical distribution of PCR curve from posterior samples 
  # of infection times
  
  res_day_lft <- as.data.table(t(res_lft$T_e)) %>%
    melt(value.name = "inf_day", variable.name = "iter")
  
  res_day_lft[ , num_id := 1:.N, iter]
  res_day_lft <- res_day_lft[, .(iter = 1:.N, inf_day), num_id]
  
  res_day_lft <- merge(test_lft, res_day_lft, by = "num_id", allow.cartesian = TRUE)
  
  res_day_lft[, x := day - round(inf_day)]
  
  res_day_lft <- res_day_lft[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
  ][, .(top = quantile(pos, 0.975), 
        mean = mean(pos),
        bottom = quantile(pos, 0.025)), by = x][x >= 0]
  
  # Add empirical distribution to posterior distribution plot
  pcols_lft <- c("Posterior distribution" = "firebrick2","Empirical Distribution" = "black")
  figS3 <- p_lft + 
    geom_ribbon(data = res_day_lft, inherit.aes = FALSE, aes(x = x, y = mean, ymin = bottom, ymax = top, fill = "Empirical Distribution"), alpha = 0.25) +
    geom_line(data = res_day_lft, inherit.aes = FALSE, aes(x = x, y = mean, lty = "Empirical mean")) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
    coord_cartesian(xlim = c(0, 30)) +
    ggtitle(paste0("LFT positivity over the course of infection (ct threshold = ", ct_threshold,")")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = pcols_lft, name = "") +
    scale_linetype_discrete(name = "")
  
  return(figS3)
}
