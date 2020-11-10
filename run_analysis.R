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

dfy <- merge(symp_final, test_final, by = c("num_id","day"))
cols1 <- scales::viridis_pal(option = "D")(10)
cols2 <- scales::viridis_pal(option = "A")(10)

# Generate figure 1
fig1 <- dfy %>%
  ggplot(aes(x = day  - last_asym_day.x, y = as.factor(num_id))) +
  geom_point(aes(fill = symptom), size = 3, shape = 21, stroke = 1.5, col = "white") +
  geom_point(aes(col = pcr_result), size = 3, shape = 21, stroke = 1.5) +
  cowplot::theme_minimal_grid() +
  scale_color_manual(values = cols1[c(2,7)], name = "PCR result", labels = c("Negative", "Positive")) +
  scale_fill_manual(values = c("white","red"), name = "Symptoms") +
  labs(x = "Day", y = "Participant ID") +
  theme(axis.text.x=element_blank()) + 
  scale_x_continuous(breaks =  seq(-19, 38, 1)) +
  coord_cartesian(xlim = c(-19, 38)) +
  new_scale_color() +
  geom_point(data = dfy[!is.na(serology_day), .(serology_day = serology_day[1] - last_asym_day.x), num_id][, sero := TRUE],
             inherit.aes = FALSE,
             aes(x = serology_day, y = as.factor(num_id), col = sero), pch = 4, size = 4) +
  scale_color_manual(name = "Serology result", values = "black", labels = "Negative") 
 

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
mod <- rstan::stan_model(here::here("pcr_breakpoint.stan"))
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
allsamp <- as.data.table(data.table::melt(res$T_e, value.name = "sample"))[, .(iterations, num_id = Var2, smp = sample)]
# Bin into discrete days
allsamp[, smp := round(smp) + start_date]
cols <- c("Posterior infection date" = "darkorchid4", "Censored onset interval" = "green4")

# Generate figure 2
fig2 <- allsamp %>%
  ggplot(aes(x = smp, y = as.factor(num_id), col = "Posterior infection date")) + 
  geom_density_ridges(stat = "binline", binwidth = 1, rel_min_height = 0.025, scale = 0.7, 
                       fill = "darkorchid4", alpha = 0.5) +
  geom_errorbarh(data = first_last_df, inherit.aes = FALSE, position = position_dodge(width = 5),
                 aes(y = as.factor(num_id), xmin = last_asym_day + start_date, xmax = first_symp_day + start_date, color = "Censored onset interval"), 
                 lty = 2, size = 1, height = 0.5
  ) +
  theme_bw() +
  labs(y = "Participant ID", x = "Date") +
  scale_color_manual(name = "", values = cols) +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14)) +
  guides(color = guide_legend(override.aes = list(fill = "white"))) +
  coord_cartesian(xlim = c(as.Date("2020-03-13"), as.Date("2020-04-17"))) +
  scale_x_date(minor_breaks = "1 day", breaks = "4 days", date_labels = "%d %B")

# Save figure 2
ggsave(fig2, filename = "figure2.pdf", width = 25, height = 20, unit = "cm")


#####################################
# Figure 3A: Posterior of PCR curve #
######################################

# Generate curve from posterior fit of model
pcols <- c("Posterior distribution" = "dodgerblue","Empirical Distribution" = "black")
pshp <- c("Posterior distribution" = 1, "Empirical Distribution" = 3)
p <- data.frame(top = apply(res$p, 2, quantile, prob = 0.975), 
                bottom = apply(res$p, 2, quantile, prob = 0.025),
                y = apply(res$p, 2, median),
                days = seq(0, 30, 0.1)) %>%
  ggplot2::ggplot(ggplot2::aes(x = days, y=y,  ymin = bottom, ymax = top, fill = "Posterior distribution")) +
  ggplot2::geom_ribbon(alpha = 0.75) + 
  cowplot::theme_cowplot() + 
  ggplot2::geom_line(aes(lty = "Posterior median")) +
  ggplot2::labs(y = "Probability of positive PCR (%)", x = "Days since infection") +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), labels = paste0(seq(0, 100, 20))) +
  ggplot2::scale_x_continuous(breaks = c(0, seq(5, 30, 5)))

# Generate empirical distribution of PCR curve from posterior samples 
# of infection times

res_day <- as.data.table(t(res$T_e)) %>%
  melt(value.name = "inf_day", variable.name = "iter")

res_day[ , num_id := 1:.N, iter]
res_day <- res_day[, .(iter = 1:.N, inf_day), num_id]

res_day <- merge(test_final, res_day, by = "num_id", allow.cartesian = TRUE)

res_day[, x := day - round(inf_day)]

res_day <- res_day[, .(pos = sum(pcr_result) / .N), by = list(iter, x)
][, .(top = quantile(pos, 0.975), 
      mean = mean(pos),
      bottom = quantile(pos, 0.025)), by = x][x >= 0]

# Add empirical distribution to posterior distribution plot
fig3a <- p + 
  geom_ribbon(data = res_day, inherit.aes = FALSE, aes(x = x, y = mean, ymin = bottom, ymax = top, fill = "Empirical Distribution"), alpha = 0.25) +
  geom_line(data = res_day, inherit.aes = FALSE, aes(x = x, y = mean, lty = "Empirical mean")) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  coord_cartesian(xlim = c(0, 30)) +
  ggtitle("PCR positivity over the course of infection") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = pcols, name = "") +
  scale_linetype_discrete(name = "")


#############################################################
# Figure 3B: Probability of detecting symptomatic infection #
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

# Incubation period distribution functions
inc_cumulative <- function(x){plnorm(x, meanlog = dat$lmean, sdlog = dat$lsd) }
inc_probability <- function(x){dlnorm(x, meanlog = dat$lmean, sdlog = dat$lsd) }

# Create data table that calculate detection probabilities for different values
# the function detpr is in the file aux_funcs.R
tab <- data.table(every = rep(rep(day_list, rep(1, length(day_list))), 2),
                  within = 30,
                  delay = rep(c(1, 2), c(5, 5)))
tab[, med := detpr(freqX = every, qu = 0.5, detect_within = within, delay_to_result = delay),
    by = c("every", "within", "delay") ]
tab[, top :=  detpr(freqX = every, qu = 0.975, detect_within = within, delay_to_result = delay),
    by = c("every", "within", "delay") ]
tab[, bottom :=  detpr(freqX = every, qu = 0.025, detect_within = within, delay_to_result = delay),
    by = c("every", "within", "delay") ]

tab[, every_lab := paste0("every ", every, " day(s)")]
tab$every_lab <- factor(tab$every_lab, levels = paste0("every ", day_list, " day(s)"))
tab$delay <- factor(tab$delay, labels = c("1 day", "2 days"))

fig3b <- tab %>%
  ggplot(aes(x = as.factor(every), y = med, col = delay, ymin = bottom, ymax = top)) +
  geom_errorbar(position = position_dodge(0.5), width = 0.5) +
  geom_point(position = position_dodge(0.5)) +
  scale_color_brewer(palette = "Set1", name = "Results delay") +
  cowplot::theme_minimal_hgrid() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 15, face = "bold")) +
  labs(x = "Testing frequency", y = "Probability (%)", title = "Probability of detecting symptomatic case before onset") +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25)) +
  coord_cartesian(ylim = c(0 ,1)) +
  geom_rect(mapping = aes(xmin = 1.99, xmax = 3.01, ymax = 0.6, ymin = 0.2), fill = NA, col = "black", lty = 2)


##############################################################
# Figure 3B: Probability of detecting asymptomatic infection #
##############################################################

# Same procedure as for 3A but with different probability
# of detection function
tab2 <- data.table(every = rep(rep(day_list, rep(length(1), length(day_list))), 2),
                  within = 7,
                  delay = rep(c(1, 2), c(5, 5)))
tab2[, med := detpr2(freqX = every, qu = 0.5, detect_within = within, delay_to_result = delay), 
    by = c("every", "within", "delay") ]
tab2[, top :=  detpr2(freqX = every, qu = 0.975, detect_within = within, delay_to_result = delay),
    by = c("every", "within", "delay") ]
tab2[, bottom :=  detpr2(freqX = every, qu = 0.025, detect_within = within, delay_to_result = delay),
    by = c("every", "within", "delay") ]

tab2[, every_lab := paste0("every ", every, " day(s)")]
tab2$every_lab <- factor(tab2$every_lab, levels = paste0("every ", day_list, " day(s)"))

# tab$within <- factor(tab$within, labels = c("Detection within 5 days", "Detection within 7 days"))
tab2$delay <- factor(tab2$delay, labels = c("1 day", "2 days"))

fig3c <- tab2 %>%
  ggplot(aes(x = every_lab, y = med, col = delay, ymin = bottom, ymax = top)) +
  geom_errorbar(position = position_dodge(0.5), width = 0.5) +
  geom_point(position = position_dodge(0.5)) +
  scale_color_brewer(palette = "Set1", name = "Results delay") + 
  cowplot::theme_minimal_hgrid() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 15, face = "bold")) +
  # facet_wrap(~ within) +
  labs(x = "Testing frequency", y = "Probability (%)", title = "Probability of detecting asymptomatic case within 7 days") +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25)) +
  coord_cartesian(ylim = c(0 ,1))

# Patchwork the three plots together
bot_panel <- (fig3b + fig3c) + patchwork::plot_layout(guides = "collect")
figure3 <- fig3a / bot_panel + plot_annotation(tag_levels = "A") 

# Save figure 3
ggsave(figure3, filename = "figure3.pdf", height = 30, width = 40, units = "cm")

########################
# Sensitivity analysis #
########################

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
    first_symp_day), by = id
  ][, .(te_upper_bound = unique(te_upper_bound)), id][,te_upper_bound]
  
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
res <- data.table::rbindlist(lapply(oos_fit, FUN = function(x){as.data.table(extract(x)$p)}), idcol = TRUE)
p_vals <- seq(0, 30, 0.1)
p_tab <- data.table::melt(res, id.vars = ".id")
p_tab$diff = (as.integer(substr(p_tab$variable, start = 2, stop = 4))) * 0.1 - 0.1
p_tab[, variable := NULL]
datt <- p_tab[, .(med = median(value), top = quantile(value, 0.975), bottom = quantile(value, 0.025)), by = c(".id", "diff")]

# Plot output
datt %>%
  ggplot(aes(x = diff, y = med, ymin = bottom, ymax = top, group = .id)) +
  geom_ribbon(fill = "blue", alpha = 1/27) +
  geom_line() +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  cowplot::theme_cowplot() +
  labs(y = "Probability of positive PCR", x = "Days since infection")
