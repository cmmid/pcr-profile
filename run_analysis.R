# Load required libraries
library(data.table)
library(magrittr)
library(ggplot2)
library(dplyr)
library(rstan)
library(ggplot2)
library(patchwork)
library(ggridges)
library(ggnewscale)

# Load auxillary functions
source(here::here("aux_funcs.R"))
source(here::here("figure_functions.R"))

# Load in participant data
test_final <- data.table::fread("test_data.csv")
symp_final <- data.table::fread("symptom_data.csv")
symp_final$date <- as.Date(symp_final$date)

## Set days relative to a chosen start date
start_date <- as.Date("2020-01-01")
test_final$date <- as.Date(test_final$date)
test_final$serology_date <- as.Date(test_final$serology_date)
test_final[, day := as.integer(date - start_date)]
test_final[, serology_day := as.integer(serology_date - start_date)]

## Add initial asymptomatic reports on enrollment day
symp_final <- rbind(symp_final, test_final[, .(date = min(date), symptom = FALSE), by = num_id])

## Find first symptomatic report dates and last asymptomatic report dates
symp_final[, first_symp := min(date[symptom == TRUE]), by = num_id]
symp_final[, last_asym := max(date[symptom == FALSE & date < first_symp]), by = num_id]

first_last_df <- symp_final[, .(first_symp_day = unique(as.integer(first_symp - start_date)),
               last_asym_day = unique(as.integer(last_asym - start_date))), by = num_id]

symp_final[, day := as.integer(date - start_date)]
symp_final[, last_asym_day := as.integer(last_asym - start_date)]
symp_final[, first_symp_day := as.integer(first_symp - start_date)]

# first_last_df[, num_id := 1:.N]

setkey(first_last_df, num_id)
setkey(test_final, num_id)

test_final <- merge(test_final, first_last_df)


######################################
# Figure 1: Symptom and testing data #
######################################

# Merge testing and symptom data into one data table
dfy <- merge(symp_final, test_final, by = c("num_id","day"), all.y = TRUE)

fig1 <- figure1(dfy) 

# Save Figure 1
ggsave(fig1, filename = "figure1.pdf", height = 20, width = 40, units = "cm")


#################
# Model fitting #
#################

# Generate list of data for stan
dat <- list()
dat$P <- first_last_df[, .N]
dat$N <- test_final[, .N]
dat$day_of_test <- test_final$day
dat$test_result <- test_final$pcr_result %>% as.numeric()
dat$patient_ID <- test_final$num_id
dat$time_first_symptom <- first_last_df$first_symp_day
dat$time_last_asym <- first_last_df$last_asym_day
dat$lmean <- EpiNow2::incubation_periods$mean
dat$lsd <- EpiNow2::incubation_periods$sd

# Upper bound on time of infection, infection must occur before
# first symptomatic report or first positive PCR, whichever
# is first
dat$te_upper_bound <- test_final[, te_upper_bound := ifelse(
  any(day[pcr_result == TRUE] < first_symp_day[pcr_result == TRUE]),
  min(day[pcr_result == TRUE & day < first_symp_day]),
  first_symp_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

# Compile and run stan model
options(mc.cores = parallel::detectCores())
mod <- rstan::stan_model("pcr_breakpoint.stan")
seedx <- 16782497
fit <- rstan::sampling(mod, chains = 4, 
                       iter = 2000,
                       warmup = 1000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.9, 
                                      stepsize = 0.75,
                                      max_treedepth = 13))
res <- rstan::extract(fit)


##########################################
# Figure 2: Posterior of infection times #
##########################################

# Create data table of all infection time samples
allsamp <- data.table::melt(as.data.table(res$T_e)[, iterations := 1:.N], 
                            id.vars = c("iterations"), 
                            value.name = "sample")[, .(iterations, num_id = variable, smp = sample)
                                                   ][, num_id := as.numeric(num_id)]
# Bin into discrete days
allsamp[, smp := round(smp) + start_date]

# Generate figure 2
fig2 <- figure2(allsamp)

# Save figure 2
ggsave(fig2, filename = "figure2.pdf", width = 25, height = 20, unit = "cm")

########################
# Figure 3A: Ct values #
########################

ct_plot_dt <- merge(test_final, allsamp[, .(infection_date = median(smp)), num_id], by = "num_id")
ct_plot_dt[, ct := ifelse(is.na(ct), 40, ifelse(ct == 0, 40, ct))]

ct_plot_dt[, x_date := date - infection_date]

fig3a <- figure3a(ct_plot_dt = ct_plot_dt, ct_threshold = 37)


#####################################
# Figure 3B: Posterior of PCR curve #
######################################

fig3b <- figure3b(res = res, test_final = test_final, ribbon_col = "dodgerblue")


#############################################################
# Figure 3C: Probability of detecting symptomatic infection #
#############################################################

# Evaluated testing frequencies
day_list <- c(1,2,4,7,14)

# Samples from PCR positive curve at different times since infection
p_tab <- as.data.table(res$p)
p_vals <- seq(0, 30, 0.1)
p_tab <- data.table::melt(p_tab)
p_tab$diff <- rep(p_vals, rep(4000, length(p_vals)))
p_tab$iter <- rep(1:4000, length(p_vals))
p_tab[, variable := NULL]

# Write PCR curve to csv
fwrite(p_tab[, .(days_since_infection = diff, sample = iter, value)], file = "PCR_curve.csv")

# Write PCR curve summary to csv
fwrite(p_tab[, .(median = median(value), 
                 lower_95 = quantile(value, 0.025), 
                 upper_95 = quantile(value, 0.975)), 
             by = list(days_since_infection = diff)],
       file = "PCR_curve_summary.csv")

# Incubation period distribution functions
inc_cumulative <- function(x){plnorm(x, meanlog = dat$lmean, sdlog = dat$lsd) }
inc_probability <- function(x){dlnorm(x, meanlog = dat$lmean, sdlog = dat$lsd) }

# Create a table of parameter value combinations to run over
tab <- data.table(every = rep(rep(day_list, rep(1, length(day_list))), 2),
                  within = 30,
                  delay = rep(c(1, 2), c(5, 5)))

# Calculate summary statistics for each parameter combination
tab <- tab[, testing_func(freqX = every, 
                          detect_within = 30, 
                                  delay_to_result = delay,
                                  symp = TRUE,
                                  ptab = p_tab), 
           by = c("every", "delay")]

# Generate figure 3c
fig3c <- figure3c(tab)


##############################################################
# Figure 3D: Probability of detecting asymptomatic infection #
##############################################################

# Same procedure as for 3A but with different probability
# of detection function
tab2 <- data.table(every = rep(rep(day_list, rep(length(1), length(day_list))), 2),
                  within = 7,
                  delay = rep(c(1, 2), c(5, 5)))

# Calculate summary statistics for each parameter combination
tab2 <- tab2[, testing_func(freqX = every, 
                          detect_within = 7, 
                          delay_to_result = delay,
                          symp = FALSE,
                          ptab = p_tab), 
           by = c("every", "delay")]

# Generate figure 3c
fig3d <- figure3d(tab2, wth = 7)

# Patchwork the three plots together
bot_panel <- (fig3c + fig3d) + patchwork::plot_layout(guides = "collect")
figure3 <- (fig3a + fig3b) / bot_panel + plot_annotation(tag_levels = "A") 

# Save figure 3
ggsave(figure3, filename = "figure3.pdf", height = 30, width = 40, units = "cm")

# Alternative figure 3 with different ct thresholds:
# THESE WILL TAKE A WHILE TO RUN AS THEY INVOLVE GENERATING 
# NEW COPIES OF FIGURE 3C + D EACH
test_final_lft_28 <- data.table::copy(test_final)
test_final_lft_25 <- data.table::copy(test_final)
figure3_28 <- fit_different_ct(ct_threshold = 28, 
                               test_final = test_final_lft_28, 
                               mod = mod, seedx = seedx)
ggsave(figure3_28$plot, filename = "figure3_28.pdf", height = 30, width = 40, units = "cm")
figure3_25 <- fit_different_ct(ct_threshold = 25, 
                               test_final = test_final_lft_25, 
                               mod = mod, seedx = seedx)
ggsave(figure3_25$plot, filename = "figure3_25.pdf", height = 30, width = 40, units = "cm")

##############################################
# Figure S1A: Comparison with other findings #
##############################################
library(tidyverse)
library(doMC)

# Directory pointing to the github repository provided by Kucirka et. al.
# from https://github.com/HopkinsIDD/covidRTPCR
kdir <- "~/repos/covidRTPCR"

n_iter <- 1500
n_warmup <- 250
p_adapt_delta <- 0.99
n_max_treedepth <- 20
## the max number of days after exposure to estimate
T_max <- 30
exposed_n <- 686
exposed_pos <- 77

source(paste0(kdir,"/R/utils.R"))

raw_data <- read_csv(paste0(kdir,"/data/antibody-test-data.csv")) %>% 
    filter(grepl("RT_PCR", test),
           study != "Danis_no_4")

pcr_dat <- raw_data %>% 
    ## add non-quantified positives to other positives for Danis et al.
    mutate(n_adj=n+nqp,
           test_pos_adj=test_pos+nqp) %>% 
    ## remove estimates without observations
    filter(n_adj > 0,
           ## days needs to be above -5
           day > -5,
           ## only use the nasal swabs from Kujawski, not throat swabs
           !(study == "Kujawski" & test == "RT_PCR_oro")) %>% 
    mutate(study_idx=paste(study, test, sep="_") %>% as.factor() %>% as.numeric(),
           pct_pos=test_pos_adj/n_adj)

day_poly <- poly(log(pcr_dat$day+5), degree=3)
poly_predict <- predict(day_poly, log(1:T_max))

npv_onset_model <- stan_model(paste0(kdir,"/Stan/npv-fixed-onset.stan"))

main_analysis <- make_analysis_data(stan_model=npv_onset_model,
                                    dat=pcr_dat,
                                    T_max=T_max,
                                    poly_est=as.matrix(day_poly),
                                    poly_pred=poly_predict,
                                    exposed_n=exposed_n,
                                    exposed_pos=exposed_pos,
                                    spec=1,
                                    iter=n_iter,
                                    warmup=n_warmup,
                                    control=list(adapt_delta=p_adapt_delta,
                                                 max_treedepth=n_max_treedepth),
                                    save_warmup=F,
                                    save_stan=F)

infections <- data.table(med = apply(res$T_e, 2, median) - first_last_df$last_asym_day,
                         bottom = apply(res$T_e, 2, quantile, prob = 0.025) - first_last_df$last_asym_day,
                         top = apply(res$T_e, 2, quantile, prob = 0.975) - first_last_df$last_asym_day,
                         num_id = 1:nrow(first_last_df)) 

infections <- merge(infections, first_last_df[, .(num_id, id)], by = "num_id")

dt <- merge(test_final, infections, by = "num_id")[, x := round((day - last_asym_day) - med)
][, xmin := round((day - last_asym_day) - bottom)
][, xmax := round((day - last_asym_day) - top)]

pcr_dat_ext <- dt[, .(test_pos = sum(pcr_result), n = .N), by = x
][order(x)
][, `:=`(study = "SAFER", test = "RT_PCR", day = x, day_min = x, day_max = x,
         days_since_exposure = x,
         inconclusive = 0, nqp = 0, pct_pos = test_pos / n, notes = NA, n_adj = n,
         test_pos_adj = test_pos, study_idx = 9, model = "Hellewell/Russell")][, x := NULL]

pcr_dat %<>%
  mutate(model = "Kucirka", days_since_exposure = day + 5)

pcr_dat_ext <- rbindlist(list(pcr_dat, pcr_dat_ext), use.names = TRUE, fill = TRUE)


kucirka_res <- as.data.table(main_analysis$plot_dat)[, .(days_since_exposure, fnr_med, fnr_ub, fnr_lb)][, model := "Kucirka"]

safer_res <- data.table(fnr_ub = 1 - apply(res$p, 2, quantile, prob = 0.975), 
                 fnr_lb = 1 - apply(res$p, 2, quantile, prob = 0.025),
                 fnr_med = 1 - apply(res$p, 2, median),
                 days_since_exposure = seq(0, 30, 0.1),
                 model = "Hellewell/Russell")

hay_res <- fread("Hay_S1_data.csv")
hay_res <- hay_res[, .(days_since_exposure = t, fnr_med = 1 - median, fnr_lb = 1 - upper, fnr_ub = 1 - lower, model = "Hay/Kennedy-Shaffer")][days_since_exposure <= 30]

borremans <- fread("borremans_swab_pos.csv")
borremans <- borremans[, .(days_since_exposure = days_since_onset + 5, pct_pos = (pct_pos / 100), n, model = "Hay/Kennedy-Shaffer")]
figS1a_points <- rbindlist(list(pcr_dat_ext, borremans), use.names = TRUE, fill = TRUE)[days_since_exposure >= 0 & days_since_exposure <= 30]

figS1a <- rbindlist(list(kucirka_res,hay_res, safer_res), use.names = TRUE)%>%
  ggplot(aes(x = days_since_exposure, y = 1 - fnr_med, ymin = 1 - fnr_lb, ymax = 1 - fnr_ub)) +
  geom_vline(xintercept = 5, lty = 2) +
  geom_line(aes(col = model)) +
  geom_ribbon(aes(fill = model), alpha = 0.45) +
  geom_point(inherit.aes = FALSE, data = figS1a_points, size = 2,
             aes(x = days_since_exposure, y = pct_pos, col = model), alpha = 0.6) +
  facet_wrap(~ model) +
  scale_x_continuous(breaks = seq(1, 30, 2)) + 
  cowplot::theme_cowplot() +
  scale_size_continuous(name = "Number of tests") +
  scale_color_brewer(guide = "none", palette = "Set2") +
  scale_fill_brewer(guide = "none", palette = "Set2") +
  labs(x = "Days since infection", y = "Probability of positive PCR test (%)") +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25))# +
  # geom_vline(data = data.frame(model = c("Hellewell/Russell", "Kucirka", "Hay/Kennedy-Shaffer"), xintercept = c(4, 8, 4)),
             # aes(xintercept = xintercept), alpha = 0.5)


##################################################
# Figure S1B: Re-fitting Kucirka with SAFER data #
##################################################

# Need to convert SAFER study days since exposure into
# days since onset
pcr_dat_ext[study == "SAFER", day := day - 5]
# Remove tests further than 4 days prior to onset
# Should only be from SAFER, negative tests
pcr_dat_ext <- pcr_dat_ext[day > -5]

day_poly_ext <- poly(log(pcr_dat_ext$day+5), degree=3)
poly_predict_ext <- predict(day_poly, log(1:T_max))

main_analysis_ext <- make_analysis_data(stan_model=npv_onset_model,
                                        dat=pcr_dat_ext,
                                        T_max=T_max,
                                        poly_est=as.matrix(day_poly_ext),
                                        poly_pred=poly_predict_ext,
                                        exposed_n=exposed_n + 200,
                                        exposed_pos=exposed_pos + 27,
                                        spec=1,
                                        iter=n_iter,
                                        warmup=n_warmup,
                                        control=list(adapt_delta=p_adapt_delta,
                                                     max_treedepth=n_max_treedepth),
                                        save_warmup=F,
                                        save_stan=F)


# Format results
kucirka_ext_res <- as.data.table(main_analysis_ext$plot_dat)[, .(days_since_exposure, fnr_med, fnr_ub, fnr_lb)
                                                ][, model := "Kucirka (original + SAFER data)"]

kucirka_res[, model := "Kucirka (original data)"]

# Format our results
safer_res[, model := "Hellewell/Russell (SAFER data only)"]


figS1b <- rbindlist(list(kucirka_res, kucirka_ext_res, safer_res), use.names = TRUE) %>%
  ggplot(aes(x = days_since_exposure, y = 1 - fnr_med, ymin = 1 - fnr_lb, ymax = 1 - fnr_ub)) +
  # geom_point(data = pcr_dat_ext, aes(x = days_since_exposure, y = pct_pos, size = n), inherit.aes = FALSE, alpha = 0.5) +
  geom_vline(xintercept = 5, lty = 2) +
  geom_line(aes(col = model)) +
  geom_ribbon(aes(fill = model), alpha = 0.25) +
  scale_x_continuous(breaks = seq(0, 30, 1)) + 
  cowplot::theme_cowplot() +
  scale_shape_discrete(name = "Number of tests") +
  # scale_fill_brewer(name = "Model", palette = "Set1") +
  # scale_colour_brewer(name = "Model", palette = "Set1") +
  scale_fill_manual(values = c("#D95F02", "dimgrey", "#7570B3")) +
  scale_colour_manual(values = c("#D95F02", "dimgrey", "#7570B3")) +
  labs(x = "Days since infection", y = "Probability of positive PCR test (%)") +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25)) +
  # geom_vline(xintercept = c(4, 6), alpha = 0.5) +
  theme(legend.position = "bottom")

# Put figure S1 together
figS1 <- figS1a / figS1b + patchwork::plot_annotation(tag_levels = "A")

# Save figure S1
ggsave(figS1, filename = "figureS1.pdf", height = 30, width = 40, units = "cm")

###########################################
# S3: LFT sensitivity curve approximation #
###########################################

# Create new synthetic PCR results with different ct threshold
test_lft <- test_final[, pcr_result := ifelse(ct == 0 | is.na(ct), 0, ifelse(ct <= 28, TRUE, FALSE))]

# Update data
dat_lft <- dat

dat_lft$test_result <- test_lft$pcr_result %>% as.numeric()
dat_lft$te_upper_bound <- test_lft[, te_upper_bound := ifelse(
  any(day[pcr_result == TRUE] < first_symp_day[pcr_result == TRUE]),
  min(day[pcr_result == TRUE & day < first_symp_day]),
  first_symp_day), by = num_id
][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]

fit_lft <- rstan::sampling(mod, chains = 4, 
                       iter = 2000,
                       warmup = 1000,
                       data = dat_lft,
                       seed = seedx,
                       control = list(adapt_delta = 0.9, 
                                      stepsize = 0.75,
                                      max_treedepth = 13))
res_lft <- rstan::extract(fit_lft)

# Samples from LFT positive curve at different times since infection
p_tab_lft <- as.data.table(res_lft$p)
p_tab_lft <- data.table::melt(p_tab_lft)
p_tab_lft$diff <- rep(p_vals, rep(4000, length(p_vals)))
p_tab_lft$iter <- rep(1:4000, length(p_vals))
p_tab_lft[, variable := NULL]

# Write LFT curve to csv
fwrite(p_tab_lft[, .(days_since_infection = diff, sample = iter, value)], file = "LFT_curve.csv")

# Write PCR curve summary to csv
fwrite(p_tab_lft[, .(median = median(value), 
                 lower_95 = quantile(value, 0.025), 
                 upper_95 = quantile(value, 0.975)), 
             by = list(days_since_infection = diff)],
       file = "LFT_curve_summary.csv")

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
pcols_lft <- c("Posterior Distribution" = "firebrick2","Empirical Distribution" = "black")
figS3 <- p_lft + 
  geom_ribbon(data = res_day_lft, inherit.aes = FALSE, aes(x = x, y = mean, ymin = bottom, ymax = top, fill = "Empirical Distribution"), alpha = 0.25) +
  geom_line(data = res_day_lft, inherit.aes = FALSE, aes(x = x, y = mean, lty = "Empirical mean")) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  coord_cartesian(xlim = c(0, 30)) +
  ggtitle("LFT positivity over the course of infection") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = pcols_lft, name = "") +
  scale_linetype_discrete(name = "")

# Save figure S3
ggsave(figS3, filename = "figureS3.pdf", height = 15, width = 20, units = "cm")


########################
# Sensitivity analysis #
########################

###############
## FIGURE S2 ##
###############

# Loop over each participant and leave them out when fitting model
oos_fit <- list()
for(i in 1:dat$P) {
  print(i)
  dat <- list()
  dat$P <- first_last_df[, .N] 
  dat$N <- test_final[num_id != i, .N] 
  dat$day_of_test <- test_final[num_id != i, day]
  dat$test_result <- test_final[num_id != i, as.numeric(pcr_result)]
  dat$patient_ID <- test_final[num_id != i, num_id]
  dat$time_first_symptom <- first_last_df[, first_symp_day]
  dat$time_last_asym <- first_last_df[, last_asym_day]
  dat$lmean <- EpiNow2::incubation_periods$mean
  dat$lsd <- EpiNow2::incubation_periods$sd
  
  dat$te_upper_bound <- test_final[, te_upper_bound := ifelse(
    any(day[pcr_result == TRUE] < first_symp_day[pcr_result == TRUE]),
    min(day[pcr_result == TRUE & day < first_symp_day]),
    first_symp_day), by = num_id
  ][, .(te_upper_bound = unique(te_upper_bound)), num_id][,te_upper_bound]
  
  oos_fit[[i]] <- rstan::sampling(mod, chains = 4, 
                                  iter = 2000,
                                  warmup = 1000,
                                  data = dat,
                                  seed = seedx,
                                  control = list(adapt_delta = 0.9, 
                                                 stepsize = 0.75,
                                                 max_treedepth = 13))
}

# Bind results together and format
res_oos <- data.table::rbindlist(lapply(oos_fit, FUN = function(x){as.data.table(rstan::extract(x)$p)}), idcol = TRUE)
p_vals <- seq(0, 30, 0.1)
p_oos <- data.table::melt(res_oos, id.vars = ".id")
p_oos$diff = (as.integer(substr(p_oos$variable, start = 2, stop = 4))) * 0.1 - 0.1
p_oos <- p_oos[, .(med = median(value), top = quantile(value, 0.975), bottom = quantile(value, 0.025)), by = c(".id", "diff")]

# Plot output
figS2 <- figureS2(data = p_oos)
ggsave(figS2, filename = "figureS2.pdf", height = 15, width = 20, units = "cm" )


##################################
## FIGURE S4 - TRAJECTORY + MAP ##
##################################

figS4 <- fig_MAP(p_tab = p_tab, seedx = seedx)
ggsave(figS4, filename = "figureS4.pdf", height = 15, width = 20, units = "cm" )

####################################
## FIGURE S5 - 5 DAY ASYMPTOMATIC ##
####################################

tab5 <- data.table(every = rep(rep(day_list, rep(length(1), length(day_list))), 2),
                   within = 5,
                   delay = rep(c(1, 2), c(5, 5)))

# Calculate summary statistics for each parameter combination
tab5 <- tab5[, testing_func(freqX = every, 
                            detect_within = within, 
                            delay_to_result = delay,
                            symp = FALSE,
                            ptab = p_tab), 
             by = c("every", "delay")]

# Generate figure 3c
figS5 <- figure3d(tab5, wth = 5)
figS5
ggsave(figS5, filename = "figureS5.pdf", height = 15, width = 20, units = "cm" )

############################################
## FIGURE S6 - FLEXIBLE INCUBATION PERIOD ##
############################################

dat$lmean <- 1.63
dat$lsd <- 0.5

fit_flex <- rstan::sampling(mod, chains = 4, 
                       iter = 4000,
                       warmup = 2000,
                       data = dat,
                       seed = seedx,
                       control = list(adapt_delta = 0.99, 
                                      stepsize = 0.75,
                                      max_treedepth = 13))
res_flex <- rstan::extract(fit_flex)

figS6 <- figureS6(res_lauer = res, res_mcaloon = res_flex)
ggsave(figS6, filename = "figureS6.pdf", height = 15, width = 20, units = "cm")
