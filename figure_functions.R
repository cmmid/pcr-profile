figure1 <- function(dfy = NULL) {
  
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

  return(fig1)  
}

figure2 <- function(allsamp = NULL) {
  
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
  
  return(fig2)
}

figure3a <- function(ct_plot_dt = NULL, ct_threshold = NULL) {
  
  fig3a <- ct_plot_dt %>%
    ggplot(aes(x = x_date, y = ct, group = num_id)) + 
    geom_point(aes(col = ifelse(ct >= ct_threshold, "Negative", "Positive"))) +
    scale_y_reverse() +
    cowplot::theme_cowplot() +
    labs(y = "Ct value", x = "Days since infection") +
    geom_vline(xintercept = 0, lty = 2) +
    geom_hline(yintercept = ct_threshold, lty = 2) +
    scale_color_manual(name = "Test result", values = c("black", "red")) +
    ggtitle(paste0("Test results with ct value threshold of ", ct_threshold)) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5))
  
  return(fig3a)
}

figure3b <- function(res = NULL, test_final = NULL, ribbon_col = NULL) {
  
  # Generate curve from posterior fit of model
  pcols <- c("Posterior Distribution" = ribbon_col,"Empirical Distribution" = "black")
  pshp <- c("Posterior Distribution" = 1, "Empirical Distribution" = 3)
  p <- data.frame(top = apply(res$p, 2, quantile, prob = 0.975), 
                  bottom = apply(res$p, 2, quantile, prob = 0.025),
                  y = apply(res$p, 2, median),
                  days = seq(0, 30, 0.1)) %>%
    ggplot2::ggplot(ggplot2::aes(x = days, y=y,  ymin = bottom, ymax = top, fill = "Posterior Distribution")) +
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
  fig3b <- p + 
    geom_ribbon(data = res_day, inherit.aes = FALSE, aes(x = x, y = mean, ymin = bottom, ymax = top, fill = "Empirical Distribution"), alpha = 0.25) +
    geom_line(data = res_day, inherit.aes = FALSE, aes(x = x, y = mean, lty = "Empirical mean")) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
    coord_cartesian(xlim = c(0, 30)) +
    ggtitle("PCR positivity over the course of infection") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = pcols, name = "") +
    scale_linetype_discrete(name = "")
  
  return(fig3b)
}

figure3c <- function(tab = NULL) {
  tab[, every_lab := paste0("every ", every, " day(s)")]
  tab$every_lab <- factor(tab$every_lab, levels = paste0("every ", day_list, " day(s)"))
  tab$delay <- factor(tab$delay, labels = c("1 day", "2 days"))
  
  fig3c <- tab %>%
    ggplot(aes(x = every_lab, y = median, col = delay, ymin = bottom, ymax = top)) +
    geom_errorbar(position = position_dodge(0.5), width = 0.5) +
    geom_point(position = position_dodge(0.5)) +
    scale_color_brewer(palette = "Set1", name = "Results delay") +
    cowplot::theme_minimal_hgrid() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "Testing frequency", y = "Probability (%)", title = "Probability of detecting symptomatic case before onset") +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25)) +
    coord_cartesian(ylim = c(0 ,1)) +
    geom_rect(mapping = aes(xmin = 1.99, xmax = 3.01, ymax = 0.6, ymin = 0.2), fill = NA, col = "black", lty = 2)
  
  return(fig3c)
}

figure3d <- function(tab2 = NULL) {
  tab2[, every_lab := paste0("every ", every, " day(s)")]
  tab2$every_lab <- factor(tab2$every_lab, levels = paste0("every ", day_list, " day(s)"))
  tab2$delay <- factor(tab2$delay, labels = c("1 day", "2 days"))
  
  fig3d <- tab2 %>%
    ggplot(aes(x = every_lab, y = median, col = delay, ymin = bottom, ymax = top)) +
    geom_errorbar(position = position_dodge(0.5), width = 0.5) +
    geom_point(position = position_dodge(0.5)) +
    scale_color_brewer(palette = "Set1", name = "Results delay") + 
    cowplot::theme_minimal_hgrid() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "Testing frequency", y = "Probability (%)", title = "Probability of detecting asymptomatic case within 7 days") +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25)) +
    coord_cartesian(ylim = c(0 ,1))
  
  return(fig3d)
}

figureS3cd <- function(tab = NULL, symp = NULL) {
  
  tt <- ifelse(symp == TRUE, 
               "Probability of detecting symptomatic case before onset",
               "Probability of detecting asymptomatic case within 7 days")
  
  tab[, every_lab := paste0("every ", every, " day(s)")]
  tab$every_lab <- factor(tab$every_lab, levels = paste0("every ", day_list, " day(s)"))
  tab$delay <- factor(tab$delay, labels = c("0 days"))
  
  figS3c <- tab %>%
    ggplot(aes(x = every_lab, y = median, ymin = bottom, ymax = top)) +
    geom_errorbar(position = position_dodge(0.5), width = 0.5) +
    geom_point(position = position_dodge(0.5)) +
    scale_color_brewer(palette = "Set1", name = "Results delay") +
    cowplot::theme_minimal_hgrid() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "Testing frequency", y = "Probability (%)", title = tt) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25)) +
    coord_cartesian(ylim = c(0 ,1))
  
  return(figS3c)
}

fig_MAP <- function(p_tab = NULL, seedx = NULL) {
  
  set.seed(seedx)
  out <- p_tab[iter %in% sample(1:4000, size = 150, replace = FALSE)] %>%
    ggplot() +
    geom_line(aes(x = diff, y = value, group = iter), col = "dodgerblue", alpha = 0.3) +
    geom_line(data = p_tab[iter == which(res$lp__ == min(res$lp__))],
              aes(x = diff, y = value), size = 1, lty = 2) +
    cowplot::theme_cowplot() +
    scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
    labs(y = "Probability of detecting infection (%)", 
         x = "Time since infection (days)")
  
  return(out)
}
