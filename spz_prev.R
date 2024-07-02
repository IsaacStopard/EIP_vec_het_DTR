rm(list = ls())

library(tidyverse); library(readxl); library(lubridate); library(mgcv)
library(rstanarm); library(patchwork)

# sporozoite prevalence
# ignore the household and location when calculating the sporozoite prevalence
df <- read_excel(path = "data/sporozoite_prev_BF.xlsx",
           sheet = "Data") %>% 
  subset(Village == "Tiefora" & Method != "MET") 

m_d <- min(df$Date)

df <- df %>% rowwise() %>% mutate(month = month(Date),
              day_y = yday(Date),
              year = year(Date),
              day = difftime(Date, m_d, units = "days")[[1]],
              spz = ifelse(Sporozoite_infection == "Negative", 0, 1))

df_plot <- df %>% group_by(Date, day, day_y, month, year) %>% summarise(tot_p = sum(spz),
                                                                 tot = n())

# human biting rate
df_c <- read_excel(path = "data/field and sporozoites_MiRA.xlsx",
                   sheet = "F1_field_.data16-19_mira") %>% subset(Village == "Tiefora" & Method != "MET") %>% 
  subset(Date <= max(df$Date))

df_c[which(df_c$Date == "2017-11-10"),"Date"] <- as.Date("2017-11-09", format = "%Y-%m-%d")
  
df_c <- df_c %>% rowwise() %>% mutate(month = month(Date),
                         day_y = yday(Date),
                         year = year(Date),
                         day = difftime(Date, m_d, units = "days")[[1]])

u_h <- unique(df_c$Hour)

# checking each household has the all the hours
df_c <- df_c %>% group_by(Date, day, day_y, month, year, Household, Location) %>% mutate(all_hours = ifelse(n() == length(u_h), 1, 0)) %>% filter(all_hours == 1) %>% 
  summarise(tot_hlc = sum(tot.gamb)) %>% ungroup()

# sporozoite prevalence

b_m <- rstanarm::stan_gamm4(formula = cbind(tot_p, tot - tot_p) ~ s(day_y, bs = "cc", k = 10),
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

new_data <- data.frame("day_y" = seq(1, 365))

p_b_m <- posterior_epred(b_m, newdata = new_data)
pred_df <- cbind(t(apply(p_b_m, 2, quantile, probs = c(0.025, 0.5, 0.975))) %>% as.data.frame(), new_data)
colnames(pred_df)[1:3] <- c("lower", "p", "upper")

ggplot() +
  geom_line(data = pred_df, aes(x = day_y, y = p)) +
  geom_point(data = df_plot, aes(x = day_y, y = tot_p/tot))

# human biting rate estimates
# constant number of births

t_m <- rstanarm::stan_gamm4(formula = tot_hlc ~ s(day_y, bs = "cc", k = 10) + Location,
                            knots=list(day_y=c(1,365)),
                            data = df_c,
                            family = neg_binomial_2,
                            prior = normal(location = 0, scale = 5),
                            prior_intercept = normal(location = 0, scale = 5),
                            prior_smooth = exponential(rate = 1, autoscale = TRUE), 
                            iter = 4000,
                            adapt_delta = 0.999,
                            cores = 4, 
                            seed = 123)

new_data_t <- rbind(new_data %>% mutate(Location = "IN"), new_data %>% mutate(Location = "OUT"))

p_t_m <- posterior_epred(t_m, newdata = new_data_t)
pred_df_t <- cbind(t(apply(p_t_m, 2, quantile, probs = c(0.025, 0.5, 0.975))) %>% as.data.frame(), new_data_t)
colnames(pred_df_t)[1:3] <- c("lower", "p", "upper")

ggplot() +
  geom_line(data = pred_df_t, aes(x = day_y, y = p, col = Location)) +
  geom_point(data = df_c, aes(x = day_y, y = tot_hlc, col = Location))

# EIR

# t_EIR <- rstanarm::stan_gamm4(formula = tot_p ~ s(day_y, bs = "cc", k = 10) + year,
#                             knots=list(day_y=c(1,365)),
#                             data = df_plot,
#                             family = neg_binomial_2,
#                             prior = normal(location = 0, scale = 5),
#                             prior_intercept = normal(location = 0, scale = 5),
#                             prior_smooth = exponential(rate = 1, autoscale = TRUE), 
#                             iter = 4000,
#                             adapt_delta = 0.99,
#                             cores = 4, 
#                             seed = 123)
# 
# t_EIR_day <- rstanarm::stan_gamm4(formula = tot_p ~ s(day_y, bs = "cc", k = 10),
#                               knots=list(day_y=c(1,365)),
#                               data = df_plot,
#                               family = neg_binomial_2,
#                               prior = normal(location = 0, scale = 5),
#                               prior_intercept = normal(location = 0, scale = 5),
#                               prior_smooth = exponential(rate = 1, autoscale = TRUE), 
#                               iter = 4000,
#                               adapt_delta = 0.99,
#                               cores = 4, 
#                               seed = 123)
# 
# p_t_EIR <- posterior_epred(t_EIR, newdata = new_data)
# p_t_EIR_day <- posterior_epred(t_EIR_day, newdata = new_data_day)

in_i <- which(new_data_t$Location == "IN")
p_t_EIR <- p_t_m[,in_i] * p_b_m

#pred_EIR <- cbind(t(apply(p_t_EIR, 2, quantile, probs = c(0.025, 0.5, 0.975))) %>% as.data.frame(), new_data)
pred_EIR <- cbind(t(apply(p_t_EIR, 2, quantile, probs = c(0.025, 0.5, 0.975))) %>% as.data.frame(), new_data)

#colnames(pred_EIR) <- c("lower", "EIR", "upper", "day", "year")
colnames(pred_EIR)[c(1:3)] <- c("lower", "EIR", "upper")

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
                                                "raw_m" = df_c,
                                                "fit_s" = b_m,
                                                "pred_s" = pred_df,
                                                "fit_m" = t_m,
                                                "pred_m" = pred_df_t,
                                                "pred_EIR" = pred_EIR
                                                ))

plot(b_m, "trace")

plot(t_m, "trace")



