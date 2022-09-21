# Description: R script to fit the Anopheles stephensi and Anopheles gambiae EIP PDFs
# Version: 0.1
# Date: 27/07/2022

rm(list = ls())
library(tidyverse); library(rstan); library(shinystan); library(cowplot); library(zipfR); library(truncnorm);library(ggpmisc);

###########################
### reading in the data ###
###########################

# gambiae data
gambiae_s_data <- read.csv(file = "vector_fits/data/ES_new_constant_temp_spz_processed.csv") %>% mutate(vector_species = "gambiae",
                                                                                                        gametocytemia = round(gametocytemia, digits = 5),
                                                                                                        ref = "Eunho_Suh")
gambiae_o_data <- read.csv(file = "vector_fits/data/ES_new_constant_temp_oocyst_processed.csv") %>% mutate(vector_species = "gambiae",
                                                                                                           gametocytemia = round(gametocytemia, digits = 5),
                                                                                                           ref = "Eunho_Suh")
gambiae_s_data <- gambiae_s_data[,c(2:ncol(gambiae_s_data))]
gambiae_o_data <- gambiae_o_data[,c(2:ncol(gambiae_o_data))]

# stephensi data
stephensi_data <- read.csv("vector_fits/data/mSOS_all_parasite_data.csv") %>% subset(DTR == 0) %>% 
  mutate(gametocytemia = ifelse(ref == "Shapiro_L_Thomas" | ref == "Murdock_Thomas" | ref == "Shapiro_Thomas", 0.08, NA))

stephensi_o_data <- subset(stephensi_data, lifestage == "cyst")
stephensi_s_data <- subset(stephensi_data, lifestage == "gland-spz")

gambiae_s_data <- gambiae_s_data[,-which(colnames(gambiae_s_data)=="Cup"|colnames(gambiae_s_data)=="DTR"|colnames(gambiae_s_data) == "Experiment")]
stephensi_s_data <- stephensi_s_data[,c("day_post_inf", "presence", "temp", "gametocytemia", "vector_species", "ref")]
colnames(stephensi_s_data)[1] <- "DPI"

s_data <- rbind(gambiae_s_data, stephensi_s_data)

s_data_in <- s_data %>% group_by(DPI, temp, gametocytemia, vector_species) %>% summarise(positive = sum(presence),
                                                                            sample = n()) %>% 
  ungroup() %>% 
  mutate(prevalence = positive/sample)

ggplot(data = s_data_in, aes(x = DPI, y = prevalence, colour = vector_species, shape = factor(gametocytemia))) + geom_point() +
  facet_wrap(~temp) + theme_bw()

gambiae_o_data %>% 

ggplot(data = stephensi_o_data, 
       aes(x = day_post_inf, 
           y = number_parasites, 
           col = vector_species, 
           group = factor(gametocytemia))) +
  geom_point() +
  facet_wrap(~temp)

# there is no infection at 17 with zero  
s_data <- s_data[-which(s_data$temp == 17 & s_data$gametocytemia == 0.00024),]
o_data <- o_data[-which(o_data$temp == 17 & o_data$gametocytemia == 0.00024),]


# survival data
survival_data <- read_survival_data("data/mSOS_all_survival_data.csv") # mosquito (A. stephensi) survival data