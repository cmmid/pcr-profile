
testing_func <- function(freqX = NULL, detect_within = NULL, 
                         delay_to_result = NULL, symp = NULL, 
                         ptab = NULL) {
  
  # For symptomatic cases, we put the delay in probability of onset by X
  # For asymptomatic cases, we put the delay in detect_within
  det <- ifelse(symp == TRUE, detect_within, detect_within - delay_to_result)
  
  # Probability of being tested on day X given testing frequency
  p_at_X <- 1 / freqX
  
  # Output matrix
  p_out <- matrix(ncol = det + 1, nrow = 4000)
  
  # Special case for testing at X = 0
  # P(tested on day 0) * P(test positive on day 0) * P(no onset by day 0)
  p_out[, 1] <- p_at_X * (ptab[diff == 0, value] * (1 - inc_cumulative(delay_to_result))) 
  
  # Loop over time since infection X
  for(X in 1:det) { # Need special case for X = 0
    
    # Probability that you have not onset by time X + delay to result from being tested now
    # Set to 1 at all times for asymptomatic infections
    p_no_onset <- ifelse(symp == TRUE, 1 - inc_cumulative(X + delay_to_result), 1)
    
    # Probability positive at X
    p_pos_X <- ptab[diff == X, value]
    
    # You could have been tested before,
    # depending on frequency
    if(X >= freqX){
      # Times that you could have been tested from X
      # to zero
      prev_t <- seq(X - freqX, 0, by = -freqX)
      
      # Probability that you tested negative at these previous times
      pt <- matrix(ncol = length(prev_t), nrow = 4000)
      for(i in 1:length(prev_t)) {
        pt[, i] <- ptab[diff == prev_t[i], 1 - value]
      }
      
      # Multiply together and add to matrix of output
      # P(tested on X) * P(no onset by the time you get the result) * 
      # P(test positive at X) * P(test negative at all times previous)
      p_out[, X + 1] <- p_at_X * p_no_onset * p_pos_X * apply(X = pt, MARGIN = 1, FUN = prod)
      
    }else{
      # If you haven't been tested before
      p_out[, X + 1] <- p_at_X * p_no_onset * p_pos_X
    }
    
  }
  
  # Summing over all times since infection for each posterior sample and taking quantiles
  out <- data.table(median = quantile(rowSums(p_out), probs = 0.5),
                    bottom = quantile(rowSums(p_out), probs = 0.025),
                    top = quantile(rowSums(p_out), probs = 0.975))
  
  return(out)
}

# Run the main analysis for a different ct value
fit_different_ct <- function(ct_threshold = NULL, test_final = NULL, mod = NULL, seedx = NULL) {
  
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
  p_vals <- seq(0, 30, 0.1)
  p_tab_lft <- as.data.table(res_lft$p)
  p_tab_lft <- data.table::melt(p_tab_lft)
  p_tab_lft$diff <- rep(p_vals, rep(12000, length(p_vals)))
  p_tab_lft$iter <- rep(1:12000, length(p_vals))
  p_tab_lft[, variable := NULL]
  p_tab_lft <- p_tab_lft[iter > 8000]
  
  # Figure S3A
  figS3a <- figure3a(ct_plot_dt = ct_plot_dt, ct_threshold = ct_threshold)
  
  # Figure S3B
  
  figS3b <- figure3b(res = res_lft, test_final = test_lft, ribbon_col = "firebrick2")
  
  # Figure S3C
  # Evaluated testing frequencies
  day_list_lft <- c(2, 3, 4)
  
  tab <- data.table(every = rep(rep(day_list_lft, rep(1, length(day_list_lft))), 1),
                    within = 30,
                    delay = 0)
  
  tab <- tab[, testing_func(freqX = every, 
                     detect_within = 30, 
                     delay_to_result = delay,
                     symp = TRUE,
                     ptab = p_tab_lft), 
      by = c("every", "delay")]
  
  figS3c <- figureS3cd(tab, symp = TRUE)
  
  # Figure S3D
  tab2 <- data.table(every = rep(rep(day_list_lft, rep(length(1), length(day_list_lft))), 1),
                     within = 7,
                     delay = 0)
  
  tab2 <- tab2[, testing_func(freqX = every, 
                              detect_within = 7, 
                              delay_to_result = delay,
                              symp = FALSE,
                              ptab = p_tab_lft), 
               by = c("every", "delay")]

  figS3d <- figureS3cd(tab2, symp = FALSE)
  

  bot_panel <- (figS3c + figS3d) + patchwork::plot_layout(guides = "collect")
  figureS3 <- (figS3a + figS3b) / bot_panel + plot_annotation(tag_levels = "A") 
  
  return(list("plot" = figureS3, "p_tab" = p_tab_lft, "fit" = fit_lft))
}
