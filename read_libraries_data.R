# script to load the commonly used packages and data

library(deSolve); library(DescTools);
library(tidyverse); library(zipfR); library(GoFKernel); library(ggnewscale)
library(cowplot); library(suncalc); library(chillR); library(ggpmisc)
library(lubridate); library(zoo); library(rio); library(RColorBrewer);library(foreach); library(doParallel)
library(readxl); library(rstan); library(dfoptim); library(viridis); library(pROC);

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
mean_temp <- 23.06936
sd_temp <- 4.361642
scaled_temps_all <- (temps_all - mean_temp) / sd_temp # scaling so on same scale as parameter fits
scaled_temps <- (temps - mean_temp) / sd_temp
# getting the EIP value
params_temp <- rstan::extract(fit)
