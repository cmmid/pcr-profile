
# A function to provide estimates of the detection probability for any continuous value of time since infection
# Have this data table in your environment
dt <- data.table::fread("fitted_params.csv")

prob_cont <- function(time_since_inf, q = 0.5) {
  dt[, p := boot::inv.logit(beta1 + beta2 * (time_since_inf - cutpoint) + (time_since_inf - cutpoint) * beta3 * beta2 * fifelse(time_since_inf - cutpoint > 0, 1, 0))]
  out <- quantile(dt$p, probs = q)
  return(out)
}

# Can pass a vector of quantiles
prob_cont(time_since_inf = 7, q = c(0.025, 0.5, 0.975))


# A nice speedy way of calculating number of detectable cases from
# a time series of daily infections

pvec <- vapply(X = 0:365, FUN = prob_cont, FUN.VALUE = 1, q = 0.5)

detectable_cases <- function(x) {
  
  out <- rep(0, length(x))
  for(i in 1:length(x)) {
    for(j in i:1) {
      out[i] <- out[i] + x[i] * pvec[i - j + 1]
    }
  }
  
  return(out)
}