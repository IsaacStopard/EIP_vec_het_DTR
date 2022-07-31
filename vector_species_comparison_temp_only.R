# Description: R script to compare the Anopheles stephensi and Anopheles gambiae EIP PDFs
# Version: 0.1
# Date: 27/07/2022

rm(list = ls())
library(tidyverse); library(rstan); library(shinystan); library(cowplot); library(zipfR); library(truncnorm);library(ggpmisc);

#################
### functions ###
#################

source(file = "utils/functions_temp_only.R")
#source(file = "utils/plotting_functions.R")

############
### fits ###
############

fit_gambiae <- readRDS("data/fit_mSOS_temp_only_f2_f3.rds")
fit_stephensi <- readRDS("data/hierarchical_mSOS_all_temperature_model_fit")

# number of iterations used to fit the gambiae data
iterations <- 5500
warmup <- 3000
chains <- 4

###########
### EIP ###
###########

# for scaling the parameters to the relevant data scaling
mean_temp_g <- 23.06936
sd_temp_g <- 4.361642

mean_temp_s <- 27.9032
sd_temp_s <- 3.471223

# temperature to simulate the EIP PDF
temps_s <- seq(21, 34, 0.1)
scaled_temps_s <- (temps_s - mean_temp_s) / sd_temp_s
temps_g <- seq(17, 30, 0.1)
scaled_temps_g <- (temps_g - mean_temp_g) / sd_temp_g

# calculating the EIP
EIP_index_g <- get_EIP(rstan::extract(fit_gambiae), scaled_temps_g, 10000)
EIP_index_s <- get_EIP_stephensi(rstan::extract(fit_stephensi), scaled_temps_s, 12000)

EIP_q_df <- rbind(cbind(gen_quantiles(EIP_index_g$EIP_10, temps_g), data.frame("EIP" = rep("EIP[10]", length(temps_g)), "vector" = rep("A. gambiae", length(temps_g)))),
                 cbind(gen_quantiles(EIP_index_g$EIP_50, temps_g), data.frame("EIP" = rep("EIP[50]", length(temps_g)), "vector" = rep("A. gambiae", length(temps_g)))),
                 cbind(gen_quantiles(EIP_index_g$EIP_90, temps_g), data.frame("EIP" = rep("EIP[90]", length(temps_g)), "vector" = rep("A. gambiae", length(temps_g)))),
                 
                 cbind(gen_quantiles(EIP_index_s$EIP_10, temps_s), data.frame("EIP" = rep("EIP[10]", length(temps_s)), "vector" = rep("A. stephensi", length(temps_s)))),
                 cbind(gen_quantiles(EIP_index_s$EIP_50, temps_s), data.frame("EIP" = rep("EIP[50]", length(temps_s)), "vector" = rep("A. stephensi", length(temps_s)))),
                 cbind(gen_quantiles(EIP_index_s$EIP_90, temps_s), data.frame("EIP" = rep("EIP[90]", length(temps_s)), "vector" = rep("A. stephensi", length(temps_s)))))

###########################
### EIP PDF differences ###
###########################

# pooled model
p_s <- rstan::extract(fit_stephensi)
p_g <- rstan::extract(fit_gambiae)

pr_g_s <- lapply(seq(21, 30, 0.2),
                run_pr,
                p_s = p_s,
                p_g = p_g,
                scaled_temps_s = scaled_temps_s,
                temps_s = temps_s,
                scaled_temps_g = scaled_temps_g,
                temps_g = temps_g)

saveRDS(pr_g_s, file = "data/prob_EIP_ga_st.rds")
pr_g_s <- readRDS(file = "data/prob_EIP_ga_st.rds")
pr_g_s_df <- bind_rows(pr_g_s) %>% mutate(temp = seq(21, 30, 0.2))

pr_g_s_plot <- ggplot(data = pr_g_s_df) +
  geom_hline(yintercept = 0.5, linetype = 2, size = 1) +
  geom_ribbon(aes(x = temp, ymin = lower, ymax = upper), fill = "grey70", alpha = 0.575) +
  geom_line(aes(x = temp, y = median), size = 1) +
  theme_bw() +
  theme(text = element_text(size = 15)) + 
  scale_x_continuous(breaks = seq(20, 30, 2)) +
  xlab("Temperature (°C)") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))  +
  ylab(expression(italic(p(EIP[g] > EIP[s]))))

#########################
### plotting the PDFs ###
#########################
t <- seq(0, 40, 0.1)
n_t <- length(t)
temps_PDF <- seq(21, 30, 3)

PDFs_df <- bind_rows(lapply(temps_PDF, function(x){
  return(rbind(calc_PDF(get_EIP_params(scaled_temps_g[which(temps_g == x)], "gambiae", p_g), t_ = t) %>% mutate(species = "A. gambiae"),
        calc_PDF(get_EIP_params(scaled_temps_s[which(temps_s == x)], "stephensi", p_s), t_ = t) %>% mutate(species = "A. stephensi")) %>% 
    mutate(temp = x))}
  )
  )
rownames(PDFs_df) <- NULL

EIP_q_df_g <- subset(EIP_q_df, temp %in% temps_PDF)

EIP_q_df_g[,"dens"] <- PDFs_df[match(interaction(round(EIP_q_df_g$median, digits = 1), EIP_q_df_g$temp, EIP_q_df_g$vector), interaction(round(PDFs_df$t, digits = 1), PDFs_df$temp, PDFs_df$species)), "median"]

EIP_q_df_g %>% spread(key = posterior, value = EIP_value)

EIP_10 <- subset(EIP_q_df, EIP == "EIP[10]")
EIP_90 <- subset(EIP_q_df, EIP == "EIP[90]")
PDFs_df[,"EIP_10"] <- round(EIP_10[match(interaction(PDFs_df$species, PDFs_df$temp), interaction(EIP_10$vector, EIP_10$temp)),"median"], digits = 1)
PDFs_df[,"EIP_90"] <- round(EIP_90[match(interaction(PDFs_df$species, PDFs_df$temp), interaction(EIP_90$vector, EIP_90$temp)),"median"], digits = 1)


PDF_plot <- ggplot(data = PDFs_df, aes(x = t, y = median)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, group = species), alpha = 0.4, fill = "grey70") + 
  geom_ribbon(data = PDFs_df %>% mutate(t = ifelse(t<EIP_10 | t>EIP_90, NA, t)) %>% na.omit(), 
              aes(ymin = 0, ymax = median, x = t, fill = species), inherit.aes=FALSE, alpha = 0.4) +
  facet_wrap(~ temp) + 
  theme_bw() + 
  theme(text = element_text(size = 15)) +
  scale_colour_manual(values = c( "#E69F00", "#56B4E9"), labels = c(parse(text="italic('A. gambiae')"), parse(text="italic('A. stephensi')")), name = "") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), labels = c(parse(text="italic('A. gambiae')"), parse(text="italic('A. stephensi')")), name = "") +
  xlab("EIP (days)") + ylab("Density") +
  geom_segment(data = subset(EIP_q_df_g, EIP == "EIP[50]"), 
               aes(x = median, xend = median, y = 0, yend = dens, col = vector), size = 0.9) +
  geom_line(aes(col = species), size = 1.1) 

png(file = "results/figures/PDF_plot.png", height = 725, width = 700)  
plot_grid(PDF_plot,
          plot_grid(pr_g_s_plot, NULL, rel_widths = c(1, 0.1875), nrow = 1),
          rel_heights = c(2, 1.5),
          nrow = 2,
          labels = c("A", "B"))
dev.off()

