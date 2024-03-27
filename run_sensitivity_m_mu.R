# Sensitivity analysis of the impact of EIP (mean and variance) on sporozoite prevalence
# Author: Isaac J Stopard
# Version: 0.01 
# Last updated: 14/06/2021
# Notes: 

rm(list = ls())
library(deSolve); library(tidyverse); library(zipfR); library(GoFKernel); library(ggnewscale); library(cowplot)
library(lubridate); library(zoo); library(viridis); library(viridisLite)

source(file = "utils/functions_temp_only.R"); source(file = "utils/model_functions.R")

###########
### EIP ###
###########
# getting the EIP values for feeds 2 and 3 only
fit <- readRDS(file = "data/fit_mSOS_temp_only_f2_f3.rds")
fit_ <- rstan::extract(fit)
comp <- 20 # number of compartments to discretize

# the experimental temps at which the EIP model was fit
# index 1 is the lowest temperature - 17
# index 11 is the highest temperature - 30
temps <- c(17, 18, 19, 20, 21, 23, 25, 27, 28, 29, 30) 
n_t <- length(temps)

n_iter <- 10000 # number iterations

i_1 <- which(temps == 20)
i_2 <- which(temps == 30)

PDF_quantiles_out <- get_PDF_multi(fit_, gt = seq(1, n_t), temp = temps)$PDF_quantiles %>%
  mutate(temp_label = paste0("Temperature: ", temp, "°C"))

# discretising the EIP PDF daily
EIP_vals_sum_t <- discretize_EIP_multi_t(diff = 1, fit_$shape_total_S[,1],
                                           fit_$rate_total_S[,1],
                                           fit_$mu[,1],
                                           fit_$k)$EIP_vals_sum %>% mutate(
                                             "temp" = temps[1],
                                             diff = 1
                                           )
for(i in 2:n_t){
  EIP_vals_sum_t <- rbind(EIP_vals_sum_t,
                        discretize_EIP_multi_t(diff = 1, fit_$shape_total_S[,i],
                                                   fit_$rate_total_S[,i],
                                                   fit_$mu[,i],
                                                   fit_$k)$EIP_vals_sum %>% # single k value
                             mutate("temp" = temps[i],
                                    diff = 1))
}

EIP_vals_sum_t <- EIP_vals_sum_t %>% mutate(temp_label = paste0("Temperature: ", temp, "°C"))

png(file = "report/discrete_PDF.png", height = 450, width = 700)
ggplot(subset(EIP_vals_sum_t, n < 50 & temp %in% c(20, 30))) +
  geom_rect(aes(xmin = lower, xmax = upper, ymin = 0, ymax = p), alpha = 0.75, fill = "steelblue1", col = "steelblue3", size = 0.25) +
  geom_line(data = subset(PDF_quantiles_out, time < 50 & temp %in% c(20, 30)),
            aes(x = time, y = mean)) + facet_wrap(~temp_label, scales = "free_x") +
  theme_bw() + xlab("EIP") + ylab("Probability density") + theme(text = element_text(size = 15))
dev.off()

###
# creating a dataframe with the MCMC parameter values for 20 and 30 degrees
params_all <- data.frame("shape" = c(fit_$shape_total_S[,i_1], fit_$shape_total_S[,i_2]),
                         "rate" = c(fit_$rate_total_S[,i_1], fit_$rate_total_S[,i_2]),
                         "mu" =  c(fit_$mu[,i_1], fit_$mu[,i_2]),
                         "k" = c(fit_$k, fit_$k),
                         "temp" = c(rep(20, n_iter), rep(30, n_iter))) %>%
  mutate(mean_ = v.int_mean_func(a = shape,
                                 b = rate,
                                 mu = mu,
                                 k = k),
         var_ = v.int_var_func(a = shape,
                               b = rate,
                               mu = mu,
                               k = k))

# checking the EIP values
params_all %>% na.omit() %>% group_by(temp) %>% summarise(mean_mean = mean(mean_),
                                                          lower_mean = quantile(mean_, 0.025),
                                                          upper_mean = quantile(mean_, 0.975),
                                                          mean_var = mean(var_),
                                                          lower_var = quantile(var_, 0.025),
                                                          upper_var = quantile(var_, 0.975)
)

# getting the mean posterior parameter values
params_df <- params_all %>% group_by(temp) %>% summarise(shape = mean(shape),
                                                         rate = mean(rate),
                                                         mu_PL = mean(mu),
                                                         k = mean(k)) %>% rowwise() %>%
  mutate(mean = int_mean_func(a = shape, b = rate,
                              mu = mu_PL, k = k),
         var = int_var_func(a = shape, b = rate,
                            mu = mu_PL, k = k))

################################
### EIP sensitivity analysis ###
################################
# possible parameter values from https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-10-202/tables/4
# a: 1/3, b: 0.35, c: 0.5, m (mosquito density per person): 0.5:100, r (recovery rate): 0.02
threshold <- 10E-7
a <- 1/3
b <- 0.35
c <- 0.35
r <- 1/100#1/7

s_params <- expand.grid("a_" = a, "b_" = b, "c_" = c, "r_" = r,
                        "m_" = seq(0.5, 40, 0.1), "mu_" = seq(5, 20, 0.1))

s_params_full <- rbind(cbind(s_params, params_df[rep(1, nrow(s_params)),]),
                       cbind(s_params, params_df[rep(2, nrow(s_params)),])) %>%
  mutate(
    #e_p_surv = v.prop_surv_int(mu = 1/mu_, shape = shape, rate = rate, mu_PL = mu_PL, k = k),
    I_m_M1 = v.eq_I_m_full(a = a_, b = b_, c = c_, mu = 1/mu_, m = m_, n = mean, r = r_),
    I_h_M1 = v.eq_I_h(a = a_, b = b_, c = c_, mu = 1/mu_, m = m_, n = mean, r = r_),
    I_m_M2 = v.continuous_eq_I_m_full(a = a_, b = b_, c = c_, mu = 1/mu_, m = m_, r = r_, shape = shape, rate = rate, mu_PL = mu_PL, k = k),
    I_h_M2 = v.continuous_eq_I_h_full(a = a_, b = b_, c = c_, mu = 1/mu_, m = m_, r = r_, shape = shape, rate = rate, mu_PL = mu_PL, k = k),
    diff_I_m = (I_m_M2 - I_m_M1),
    EIR_M1 = I_m_M1 * m_ * a_, # so per day
    EIR_M2 = I_m_M2 * m_ * a_, # so per day
    diff_EIR = EIR_M2 - EIR_M1,
    R0_M1 = v.R0(a = a_, b = b_, c = c_, mu = 1/mu_, m = m_, r = r_, n = mean),
    R0_M2 = v.continuous_R0(a = a_ , b = b_, c = c_, mu = 1/mu_, m = m_, r = r_, shape = shape, rate = rate, mu_PL = mu_PL, k = k),
    diff_R0 = R0_M2 - R0_M1,
    I_m_M1_not_full = v.eq_I_m(a = a_, c = c_, I_h = I_h_M1, mu = 1/mu_, n = mean),
    temp_facet = ifelse(temp == 20, "Temperature: 20°C", "Temperature: 30°C")
  )

s_params_relative <- s_params_full %>% mutate(diff_I_m = diff_I_m / I_m_M2,
                                            diff_EIR = diff_EIR / EIR_M2,
                                            diff_R0 = diff_R0 / R0_M2)

s_params_full <- s_params_full[-which(s_params_full$I_h_M1 < 0 & s_params_full$I_m_M1 < 0 & s_params_full$I_m_M2 < 0 & s_params_full$I_h_M2 < 0),]

s_params_full <- s_params_full %>% mutate(I_h_M1 = ifelse(I_h_M1 < threshold | I_m_M1 < threshold, NA, I_h_M1),
                         I_m_M1 = ifelse(I_h_M1 < threshold | I_m_M1 < threshold, NA, I_m_M1),
                         diff_I_m = ifelse(I_h_M1 < threshold | I_m_M1 < threshold, NA, diff_I_m),
                         EIR_M1 = ifelse(I_h_M1 < threshold | I_m_M1 < threshold, NA, EIR_M1),
                         diff_EIR = ifelse(I_h_M1 < threshold | I_m_M1 < threshold, NA, diff_EIR),
                         R0_M1 = ifelse(I_h_M1 < threshold | I_m_M1 < threshold, NA, R0_M1),
                         diff_R0 = ifelse(I_h_M1 < threshold | I_m_M1 < threshold, NA, diff_R0))

s_params_relative <- s_params_relative[-which(s_params_relative$I_h_M1 < threshold & 
                                                s_params_relative$I_m_M1 < threshold & 
                                                s_params_relative$I_m_M2 < threshold & 
                                                s_params_relative$I_h_M2 < threshold),]

s_params_relative <- s_params_relative %>% mutate(I_h_M1 = ifelse(I_h_M1 < threshold | I_m_M1 < threshold, NA, I_h_M1),
                                          I_m_M1 = ifelse(I_h_M1 < threshold | I_m_M1 < threshold, NA, I_m_M1),
                                          diff_I_m = ifelse(I_h_M1 < threshold | I_m_M1 < threshold, NA, diff_I_m),
                                          EIR_M1 = ifelse(I_h_M1 < threshold | I_m_M1 < threshold, NA, EIR_M1),
                                          diff_EIR = ifelse(I_h_M1 < threshold | I_m_M1 < threshold, NA, diff_EIR),
                                          R0_M1 = ifelse(I_h_M1 < threshold | I_m_M1 < threshold, NA, R0_M1),
                                          diff_R0 = ifelse(I_h_M1 < threshold | I_m_M1 < threshold, NA, diff_R0))

########################
### model comparison ###
########################

# controlling the human malaria prevalence
m_params <- expand.grid("a_" = a, "b_" = b, "c_" = c,
                        "mu_" = seq(3, 20, 0.1), "m_" = c(1, 20, 40), "r_" = r)

m_params_ <- rbind(cbind(m_params, params_df[rep(1, nrow(m_params)),]),
                   cbind(m_params, params_df[rep(2, nrow(m_params)),])) %>%
  
  mutate(
         classic = v.eq_I_m_full(a = a_, 
                               b = b_,
                               c = c_,
                               mu = 1/mu_, 
                               m = m_,
                               n = mean,
                               r = r_),
         multiple = v.continuous_eq_I_m_full(a = a_,
                                           b = b_,
                                           c = c_, 
                                           mu = 1/mu_,
                                           m = m_,
                                           r = r_,
                                           shape = shape, 
                                           rate = rate,
                                           mu_PL = mu_PL, 
                                           k = k)) %>% rowwise() %>%
  mutate(temp_ = temp, 
         discrete_ep = discrete_eq_I_m_full(a = a_, 
                                         b = b_,
                                         c = c_, 
                                         mu = 1/mu_,
                                         m = m_,
                                         r = r_,
                                         n = subset(EIP_vals_sum_t, temp == temp_)$n,
                                         p = subset(EIP_vals_sum_t, temp == temp_)$p),
         discrete = discrete_eq_I_m_full(a = a_, 
                                            b = b_,
                                            c = c_, 
                                            mu = 1/mu_,
                                            m = m_,
                                            r = r_,
                                            n = subset(EIP_vals_sum_t, temp == temp_)$n,
                                            p = subset(EIP_vals_sum_t, temp == temp_)$p),
         R0_discrete = discrete_R0(a = a_, 
                                      b = b_,
                                      c = c_, 
                                      mu = 1/mu_,
                                      m = m_,
                                      r = r_,
                                      n = subset(EIP_vals_sum_t, temp == temp_)$n,
                                      p = subset(EIP_vals_sum_t, temp == temp_)$p))  

### running the Erlang distribution model to equilibrium
for(i in 1:nrow(m_params_)){
  print(i)
  params <- c(beta = 0,
              mu = 1/m_params_[i, "mu_"][[1]], # assuming there is no mosquito emergence
              a = a,
              b = b,
              c = c,
              m = m_params_[i, "m_"][[1]],
              shape = 47,
              rate = 47/m_params_[i, "mean"][[1]],
              r = r
              )
  
  state <- c(#S = 1,
             E = rep(0, params[["shape"]]),
             I = 0,
             I_h = 0.01
  )
  
  out_mc <- as.data.frame(ode(y = state, times = seq(0, 1500, 1), func = multi_comp_model_cb,
                              parms = params)) %>% rowwise() %>%
    mutate(#M = S + sum(c_across("E1":paste0("E",params[["shape"]]))) + I,
           s_prev = I)
  
  m_params_[i, "model_3"] <- out_mc[1500, "s_prev"][[1]]
}

ggplot(data = out_mc, aes(x = time, y = s_prev)) + geom_line()

saveRDS(m_params_, file = "data/m_params_20.rds")
m_params_ <- readRDS(file = "data/m_params_20.rds")

m_params_ <- m_params_ %>% mutate(classic = ifelse(classic < threshold, 0, classic),
                     discrete = ifelse(discrete < threshold, 0, discrete),
                     multiple = ifelse(multiple < threshold, 0, multiple),
                     R0_classic =  v.R0(a = a_, b = b_, c = c_, mu = 1/mu_, m = m_, r = r_, n = mean),
                     R0_M2 = v.continuous_R0(a = a_ , b = b_, c = c_, mu = 1/mu_, m = m_, r = r_, shape = shape, rate = rate, mu_PL = mu_PL, k = k),
                     R0_classic = ifelse(R0_classic < threshold, 0, R0_classic),
                     R0_M2 = ifelse(R0_M2 < threshold, 0, R0_M2),
                     R0_discrete = ifelse(R0_discrete < threshold, 0, R0_discrete),
                     EIR_M1 = classic * m_ * a_, # so per year
                     EIR_M2 = multiple * m_ * a_, # so per year
                     EIR_discrete = discrete * m_ * a_,
                     EIR_M3 = model_3 * m_ * a_,
                     classic_10 = v.eq_I_m_full(a = a_, b = b_, c = c_, mu = 1/mu_, m = m_, n = 10, r = r_),
                     classic_10 = ifelse(classic_10 < threshold, 0, classic_10),
                     EIR_10 = classic_10 * m_ * a_,
                     #classic_10 = ifelse(classic_10 > 0.2 & temp == 20, NA, classic_10),
                     #EIR_10 = ifelse(EIR_10 > 3  & temp == 20, NA, EIR_10),
                     R0_10 = v.R0(a = a_, b = b_, c = c_, mu = 1/mu_, m = m_, n = 10, r = r_),
                     R0_10 = ifelse(R0_10 < threshold, 0, R0_10)) %>% 
  subset(m_!=20 & mu_ <= 15) %>% 
         mutate(m_facet = ifelse(m_ == 1, "m = 1", ifelse(m_ == 20, "m = 20", "m = 40")),
                temp_facet = ifelse(temp == 20, "Temperature: 20°C", "Temperature: 30°C"))

# comparing the effects of different mean EIPs
s_prev_plot <- ggplot(data = m_params_ %>% gather(key = "model", value = "S", classic, multiple, classic_10), 
                      aes(x = mu_, y = S, col = model, linetype = m_facet)) +
  geom_line(size = 1.25, alpha = 0.7) + 
  facet_wrap(~temp_facet) + 
  theme_bw() + theme(text = element_text(size = 20)) +
  scale_colour_manual(name = "Model", values = c("black", "#CC79A7", "#56B4E9"),
                      labels = c("model M1 (temp-dependent)", "model M1 (temp-independent)","model M2C")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.4)) +
  ylab(expression(atop("Equilibrium sporozoite", paste("prevalence (",italic(z),")")))) + xlab("Mean mosquito life expectancy (days)") +
  scale_linetype_manual(values = c("twodash", "solid"), name = "Vector-host ratio") +
  scale_x_continuous(breaks = seq(0, 15, 5))

EIR_plot <- ggplot(data = m_params_ %>% gather(key = "model", value = "EIR", EIR_M1, EIR_M2, EIR_10), 
       aes(x = mu_, y = EIR, col = model, linetype = m_facet)) +
  geom_line(size = 1.25, alpha = 0.7) + 
  facet_wrap(~temp_facet) + 
  theme_bw() + theme(text = element_text(size = 20)) +
  scale_colour_manual(name = "Model", values = c("#CC79A7", "black", "#56B4E9"),
                      labels = c("model M1 (temp-independent)", "model M1 (temp-dependent)", "model M2C")) +
  ylab(expression(atop("Entomological inoculation rate", paste("(EIR) / bites ", day^{-1})))) + xlab("Mean mosquito life expectancy (days)") +
  scale_linetype_manual(values = c("twodash", "solid"), name = "Vector-host ratio") +
  scale_x_continuous(breaks = seq(0, 15, 5))

png(file = "results/figures/temp_changes_mean_EIP.png", height = 700, width = 1000)
plot_grid(
  plot_grid(
    s_prev_plot + theme(legend.position = "none"),
    EIR_plot + theme(legend.position = "none"), labels = c("A", "B"), nrow = 2, label_size = 20
    ),
  get_legend(EIR_plot), rel_widths = c(0.75, 0.25)
)
dev.off()

theme_c <- theme(text = element_text(size = 16), 
                 plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                 plot.title = element_text(size = 17))

png(file = "results/figures/Fig_EIP_heterogeneity.png", height = 1000, width = 700)
plot_grid(
  ggplot(data = s_params_full %>% subset(mu_ <= 15), 
       aes(x = mu_, y = m_, z = diff_I_m, group = as.factor(temp))) + 
  facet_wrap(vars(temp_facet)) +
  geom_raster(aes(fill = diff_I_m)) + 
  stat_contour(data = na.omit(s_params_full) %>% subset(mu_ <= 15), col = "black", bins = 10, size = 0.75, alpha = 0.5) +
  #  geom_text_contour(data = na.omit(s_params_full), col = "white", check_overlap = "TRUE", skip = 0, aes(z = diff_I_m*100),
  #                    rotate = FALSE) + 
  ylab("Mosquito abundance per person") + xlab("Mean mosquito life expectancy (days)") +
  # scale_fill_gradient(low = "#56B4E9", high = "#CC79A7", 
  #                     name = "Difference\nin sporozoite\nprevalence",
  #                     na.value = 'grey80',
  #                     labels = scales::percent_format(accuracy = 0.01)) + 
  scale_fill_viridis_c(option = "viridis",
                      name = "",#"Difference\nin sporozoite\nprevalence",
                      na.value = 'grey80',
                      labels = scales::percent_format(accuracy = 0.01)) +
  scale_x_continuous(breaks = seq(5, 20, 5)) +
  theme_bw() + theme_c +
    ggtitle("Difference in sporozoite prevalence"),
  
  ggplot(data = s_params_relative %>% subset(mu_ <= 15), 
       aes(x = mu_, y = m_, z = diff_I_m, group = as.factor(temp))) + 
  facet_wrap(vars(temp_facet)) +
  geom_raster(aes(fill = diff_I_m)) + 
  stat_contour(data = na.omit(s_params_relative) %>% subset(mu_ <= 15), col = "black", bins = 10, size = 0.75, alpha = 0.5) +
  ylab("Mosquito abundance per person") + xlab("Mean mosquito life expectancy (days)") +
  scale_fill_viridis_c(option = "viridis",#low = "#56B4E9", high = "#009E73", #"#CC79A7"
                      name = "", #name = "Normalised\ndifference\nin sporozoite\nprevalence",
                      na.value = 'grey80',
                      labels = scales::percent,
                      limits = c(0, 1)) + 
  scale_x_continuous(breaks = seq(5, 20, 5)) +
  theme_bw() + theme_c +
    ggtitle("Normalised difference in sporozoite prevalence"),
  
  ggplot(data = s_params_full %>% subset(mu_ <= 15), 
         aes(x = mu_, y = m_, z = diff_EIR, group = as.factor(temp))) + 
    facet_wrap(vars(temp_facet)) +
    geom_raster(aes(fill = diff_EIR)) + 
    stat_contour(data = na.omit(s_params_full) %>% subset(mu_ <= 15), col = "black", bins = 10, size = 0.75, alpha = 0.5) +
    ylab("Mosquito abundance per person") + xlab("Mean mosquito life expectancy (days)") +
    scale_fill_viridis_c(option = "viridis", #low = "#56B4E9", high = "#E69F00", 
                        name = "",#name = "Difference\nin EIR",
                        na.value = 'grey80') + 
    scale_x_continuous(breaks = seq(5, 20, 5)) +
    theme_bw() + theme_c +
    ggtitle("Difference in the EIR"),

  ncol = 1, labels = c("A", "B", "C")
)
dev.off()

# initial values
get_vals <- function(df = m_params_, temp_in, m_in, value, variable){
  df_ <- subset(df, temp == temp_in & m_ == m_in)
  df_ <- df_[which(df_[,variable] == value),]
  return(df_)
}

max(na.omit(subset(s_params_full, temp == 30)$diff_I_m))*100
max(na.omit(subset(s_params_full, temp == 20)$diff_EIR))
max(na.omit(subset(s_params_full, temp == 20)$diff_R0))
max(na.omit(subset(s_params_full, temp == 30)$diff_EIR))
max(na.omit(subset(s_params_full, temp == 30)$diff_R0))

subset(s_params_full, temp == 30) %>% group_by(mu_) %>% summarise(m = mean(na.omit(diff_I_m))) %>% 
  ungroup() %>% top_n(1, m)

get_vals(temp_in = 20, m_in = 20, value = 10, variable = "mu_")$EIR_M1

get_vals(df = s_params_relative, temp_in = 20, m_in = 20, value = 6.9, variable = "mu_")$diff_I_m



m_p_30 <- subset(m_params_, temp == 30 & m_ == 20)
min(m_p_30[which(m_p_30$classic_10 > 0),]$mu_)

m_p_17 <- subset(m_params_, temp == 20 & m_ == 20)
min(m_p_17[which(m_p_17$classic > 0),]$mu_)

min(m_p_30[which(m_p_30$classic_10 > 0),]$mu_)

max(get_vals(temp_in = 30, m_in = 20, value = 0, variable = "classic")$mu_)

max(subset(s_params_full, temp == 17 & m_ == 37)[which(subset(s_params_full, temp == 17 & m_ == 37)$I_m_M2 == min(unique(subset(s_params_full, temp == 17 & m_ == 37)$I_m_M2))), "mu_"])

max(subset(s_params_full, temp == 17 & m_ == 37)[which(is.na(subset(s_params_full, temp == 17 & m_ == 37)$I_m_M2) == 1), "mu_"])

subset(m_params_, temp == 17 & m_ == 37)[which(subset(m_params_, temp == 17 & m_ == 37)$multiple > 0),]

max(na.omit(subset(s_params_full, temp == 20)$diff_EIR)) * 100

s_params_full[which(s_params_full$diff_I_m == max(na.omit(subset(s_params_full, temp == 20)$diff_I_m))),]

s_params_full[which(s_params_full$diff_I_m == max(na.omit(subset(s_params_full, temp == 30)$diff_I_m))),]

max(na.omit(subset(s_params_full, temp == 30)$diff_R0))
