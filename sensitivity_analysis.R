rm(list = ls())
library(deSolve); library(tidyverse); library(zipfR); library(GoFKernel); library(ggnewscale); library(cowplot)
library(lubridate); library(zoo); library(viridis); library(pROC)

source(file = "utils/data_functions.R"); source(file = "utils/model_functions.R")

###########
### EIP ###
###########
# getting the EIP values
fit <- readRDS(file = "data/fit_mSOS_temp_only_f2_f3.rds")
fit_ <- rstan::extract(fit)
# index 1 is the lowest temperature - 17
# index 11 is the highest temperature - 30
# single k value
temps <- c(17, 18, 19, 20, 21, 23, 25, 27, 28, 29, 30)
n_t <- length(temps)

# for all temperatures
temps_all <- seq(17, 30, 0.01)

# temperatures for scaling
mean_temp <- 23.06936
sd_temp <- 4.361642
scaled_temps_all <- (temps_all - mean_temp) / sd_temp # scaling so on same scale as parameter fits
scaled_temps <- (temps - mean_temp) / sd_temp

# getting the EIP value
params_temp <- rstan::extract(fit)

EIP_index <- get_EIP(params_temp, scaled_temps_all, 10000)
EIP_10 <- gen_quantiles(EIP_index$EIP_10, temps_all)
EIP_50 <- gen_quantiles(EIP_index$EIP_50, temps_all)
EIP_90 <- gen_quantiles(EIP_index$EIP_90, temps_all)


