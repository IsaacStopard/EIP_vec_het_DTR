rm(list = ls())

library(tidyverse); library(readxl); library(lubridate); library(mgcv)
library(rstanarm); library(patchwork)

df <- read_excel(path = "data/sporozoite_prev_BF.xlsx",
           sheet = "Data") %>% 
  subset(Village == "Tiefora" & Method != "MET") 

m_d <- min(df$Date)

df <- df %>% rowwise() %>% mutate(month = month(Date),
              day_y = yday(Date),
              year = year(Date),
              day = difftime(Date, m_d, units = "days")[[1]],
              spz = ifelse(Sporozoite_infection == "Negative", 0, 1)) 

df_plot <- df %>% group_by(day, Date, day_y, year) %>% summarise(tot_p = sum(spz),
                                              tot = n())

b_m <- rstanarm::stan_gamm4(formula = cbind(tot_p, tot - tot_p) ~ s(day_y, bs = "cc", k = 10) + year,
                            knots = list(day_y = c(1, 365)),
                            data = df_plot, 
                            family = binomial("logit"),
                            prior = normal(location = 0, scale = 5),
                            prior_intercept = normal(location = 0, scale = 5),
                            prior_smooth = exponential(rate = 1, autoscale = TRUE), 
                            iter = 4000,
                            cores = 4,
                            adapt_delta = 0.99,
                            seed = 123)

b_m_day <- rstanarm::stan_gamm4(formula = cbind(tot_p, tot - tot_p) ~ s(day_y, bs = "cc", k = 10),
                            knots = list(day_y = c(1, 365)),
                            data = df_plot, 
                            family = binomial("logit"),
                            prior = normal(location = 0, scale = 5),
                            prior_intercept = normal(location = 0, scale = 5),
                            prior_smooth = exponential(rate = 1, autoscale = TRUE), 
                            iter = 4000,
                            cores = 4,
                            adapt_delta = 0.999,
                            seed = 123)

new_data <- expand.grid("day_y" = (seq(1, 365)), "year" = seq(2016, 2018))

new_data_day <- data.frame("day_y" =seq(1, 365))

p_b_m <- posterior_epred(b_m, newdata = new_data)
p_b_m_day <- posterior_epred(b_m_day, newdata = new_data_day)

pred_df <- cbind(t(apply(p_b_m, 2, quantile, probs = c(0.025, 0.5, 0.975))) %>% as.data.frame(), new_data)
pred_df_day <- cbind(t(apply(p_b_m_day, 2, quantile, probs = c(0.025, 0.5, 0.975))) %>% as.data.frame(), new_data_day)

colnames(pred_df) <- c("lower", "p", "upper", "day", "year")
colnames(pred_df_day) <- c("lower", "p", "upper", "day")

t_m <- rstanarm::stan_gamm4(formula = tot ~ s(day_y, bs = "cc", k = 10) + year,
                            knots=list(day_y=c(1,365)),
                            data = df_plot,
                            family = neg_binomial_2,
                            prior = normal(location = 0, scale = 5),
                            prior_intercept = normal(location = 0, scale = 5),
                            prior_smooth = exponential(rate = 1, autoscale = TRUE), 
                            iter = 4000,
                            adapt_delta = 0.99,
                            cores = 4, 
                            seed = 123)

t_m_day <- rstanarm::stan_gamm4(formula = tot ~ s(day_y, bs = "cc", k = 10),
                            knots=list(day_y=c(1,365)),
                            data = df_plot,
                            family = neg_binomial_2,
                            prior = normal(location = 0, scale = 5),
                            prior_intercept = normal(location = 0, scale = 5),
                            prior_smooth = exponential(rate = 1, autoscale = TRUE), 
                            iter = 4000,
                            adapt_delta = 0.99,
                            cores = 4, 
                            seed = 123)

p_t_m <- posterior_epred(t_m, newdata = new_data)
p_t_m_day <- posterior_epred(t_m_day, newdata = new_data_day)

pred_df_t <- cbind(t(apply(p_t_m, 2, quantile, probs = c(0.025, 0.5, 0.975))) %>% as.data.frame(), new_data)
pred_df_t_day <- cbind(t(apply(p_t_m_day, 2, quantile, probs = c(0.025, 0.5, 0.975))) %>% as.data.frame(), new_data_day)

colnames(pred_df_t) <- c("lower", "p", "upper", "day", "year")
colnames(pred_df_t_day) <- c("lower", "p", "upper", "day")

# EIR
t_EIR <- rstanarm::stan_gamm4(formula = tot_p ~ s(day_y, bs = "cc", k = 10) + year,
                            knots=list(day_y=c(1,365)),
                            data = df_plot,
                            family = neg_binomial_2,
                            prior = normal(location = 0, scale = 5),
                            prior_intercept = normal(location = 0, scale = 5),
                            prior_smooth = exponential(rate = 1, autoscale = TRUE), 
                            iter = 4000,
                            adapt_delta = 0.99,
                            cores = 4, 
                            seed = 123)

t_EIR_day <- rstanarm::stan_gamm4(formula = tot_p ~ s(day_y, bs = "cc", k = 10),
                              knots=list(day_y=c(1,365)),
                              data = df_plot,
                              family = neg_binomial_2,
                              prior = normal(location = 0, scale = 5),
                              prior_intercept = normal(location = 0, scale = 5),
                              prior_smooth = exponential(rate = 1, autoscale = TRUE), 
                              iter = 4000,
                              adapt_delta = 0.99,
                              cores = 4, 
                              seed = 123)

p_t_EIR <- posterior_epred(t_EIR, newdata = new_data)
p_t_EIR_day <- posterior_epred(t_EIR_day, newdata = new_data_day)

pred_EIR <- cbind(t(apply(p_t_EIR, 2, quantile, probs = c(0.025, 0.5, 0.975))) %>% as.data.frame(), new_data)
pred_EIR_day <- cbind(t(apply(p_t_EIR_day, 2, quantile, probs = c(0.025, 0.5, 0.975))) %>% as.data.frame(), new_data_day)

colnames(pred_EIR) <- c("lower", "EIR", "upper", "day", "year")
colnames(pred_EIR_day) <- c("lower", "EIR", "upper", "day")

# t_p_t_m <- cbind(t(p_t_m), new_data)

# t_p_t_m <- t_p_t_m %>% pivot_longer(1:8000, names_to = "iteration",
#                          values_to = "m") %>% 
#   group_by(day_y, iteration) %>% summarise(n_m = n(), m = sum(m) / n_m) %>% ungroup() %>% 
#   group_by(day_y) %>% summarise(m_m = median(m),
#                                 m_l = quantile(m, probs = 0.025)[[1]],
#                                 m_u = quantile(m, probs = 0.975)[[1]]) %>% 
#   rename(day = day_y)
# 
# pred_EIR <- cbind(t(apply(p_b_m * p_t_m, 2, quantile, probs = c(0.025, 0.5, 0.975))) %>% as.data.frame(), new_data)
# 
# colnames(pred_EIR) <- c("lower", "EIR", "upper", "day", "year")

saveRDS(file = "data/format_spz_data.rds", list("raw" = df_plot, 
                                                "fit_s" = b_m,
                                                "pred_s" = pred_df,
                                                "fit_m" = t_m,
                                                "pred_m" = pred_df_t,
                                                #"t_p_t_m" = t_p_t_m,
                                                "t_EIR" = t_EIR,
                                                "pred_EIR" = pred_EIR,
                                                "fit_s_day" = b_m_day,
                                                "pred_s_day" = pred_df_day,
                                                "fit_m_day" = t_m_day,
                                                "pred_m_day" = pred_df_t_day,
                                                #"t_p_t_m" = t_p_t_m,
                                                "t_EIR_day" = t_EIR_day,
                                                "pred_EIR_day" = pred_EIR_day
                                                ))

plot(b_m, "trace")

plot(t_m, "trace")

plot(t_EIR, "trace")

mos_plot <- ggplot() + 
  geom_point(data = spz_data$raw, aes(x = date_p, y = tot, col = factor(year)), size = 3.5) +
  xlab("Day of the year") +
  ylab("Human biting rate\nper person day") +
  geom_line(data = spz_data$pred_m, aes(x = date_p, y = p, col = factor(year)), linewidth = 1) +
  geom_ribbon(data = spz_data$pred_m, aes(x = date_p, ymin = lower, ymax = upper, fill = factor(year)), alpha = 0.1) +
  theme_bw() + theme(text = element_text(size = 18)) +
  scale_colour_manual(name = "year", values = c("#000000", "#0072B2", "#009E73")) +
  scale_fill_manual(name = "year", values = c("#000000", "#0072B2", "#009E73")) +
  scale_x_date(date_labels = "%d %b",
               date_breaks = "1 month"#,
               #limits = as.Date(c("01/01/2017", "01/01/2018"), format = "%d/%m/%Y")
  )

spz_plot <- ggplot() + 
  geom_point(data = spz_data$raw, aes(x = date_p, y = tot_p/tot, col = factor(year)), size = 3.5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.05)) +
  xlab("Day of the year") +
  ylab("Sporozoite prevalence") +
  geom_line(data = spz_data$pred_s, aes(x = date_p, y = p, col = factor(year)), linewidth = 1) +
  geom_ribbon(data = spz_data$pred_s, aes(x = date_p, ymin = lower, ymax = upper, fill = factor(year)), alpha = 0.1) +
  theme_bw() + theme(text = element_text(size = 18)) +
  #scale_size_continuous(name = "sample size") +
  scale_colour_manual(name = "year", values = c("#000000", "#0072B2", "#009E73")) +
  scale_fill_manual(name = "year", values = c("#000000", "#0072B2", "#009E73")) +
  scale_x_date(date_labels = "%d %b",
               date_breaks = "1 month",
               limits = as.Date(c("01/01/2017", "01/01/2018"), format = "%d/%m/%Y")) +
  coord_cartesian(ylim = c(0, 0.2))

ggplot() + 
  geom_point(data = df_plot, aes(x = day_y, y = tot_p/tot, size = tot, col = factor(year))) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.05)) +
  xlab("Day of the year") +
  ylab("Sporozoite prevalence") +
  geom_line(data = pred_df, aes(x = day, y = p, col = factor(year))) +
  geom_ribbon(data = pred_df, aes(x = day, ymin = lower, ymax = upper, fill = factor(year)), alpha = 0.15) +
  theme_bw() + theme(text = element_text(size = 18)) +
  scale_size_continuous(name = "sample size") +
  geom_line(data = pred_df_day, aes(x = day, y = p)) +
  geom_ribbon(data = pred_df_day, aes(x = day, ymin = lower, ymax = upper), alpha = 0.15)

ggplot() + 
  geom_point(data = df_plot, aes(x = day, y = tot_p/tot, size = tot)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.05)) +
  xlab("Day of the year") +
  ylab("Sporozoite prevalence") +
  geom_line(data = pred_df, aes(x = year, y = p)) +
  geom_ribbon(data = pred_df, aes(x = year, ymin = lower, ymax = upper), alpha = 0.15) +
  theme_bw() + theme(text = element_text(size = 18)) +
  scale_size_continuous(name = "sample size")


ggplot() + 
  geom_point(data = df_plot, aes(x = day_y, y = tot, col = factor(year))) +
  xlab("Day of the year") +
  ylab("Mosquito biting rate") +
  geom_line(data = pred_df_t, aes(x = day, y = p, col = factor(year))) +
  geom_ribbon(data = pred_df_t, aes(x = day, ymin = lower, ymax = upper, fill = factor(year)), alpha = 0.15) +
  theme_bw() + theme(text = element_text(size = 18))

ggplot() + 
  geom_point(data = df_plot, aes(x = day_y, y = tot_p, col = factor(year))) +
  xlab("Day of the year") +
  ylab("EIR") +
  geom_line(data = pred_EIR, aes(x = day, y = EIR, col = factor(year))) +
  geom_ribbon(data = pred_EIR, aes(x = day, ymin = lower, ymax = upper, fill = factor(year)), alpha = 0.15) +
  theme_bw() + theme(text = element_text(size = 18))

