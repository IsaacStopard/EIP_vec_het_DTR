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
source(file = "utils/vector_comp_functions.R")

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


calc_EIP_v <- function(EIP_index, temps, n_iter){
  
  set.seed(12345)
  x <- runif(n_iter)
  n_temp <- length(temps)
  
  shape_total_S <- EIP_index$shape_total_S
  rate_total_S <- EIP_index$rate_total_S
  mu <- EIP_index$mu
  k <- EIP_index$k
  
  out <- bind_rows(
    lapply(seq(1, n_temp),
         function(j, n_iter, temp, shape_total_S, rate_total_S, mu, k, x){
           print(temp[j])
           EIP <- lapply(seq(1, n_iter), 
                         function(i, j, temp, shape_total_S, rate_total_S, mu, k, x){
                           Inv_EIP_CDF(shape_total_S[i, j], rate_total_S[i, j], mu[i, j], k[i, j], x)
                           }, 
                         j = j,
                         temp = temp,
                         shape_total_S = shape_total_S, 
                         rate_total_S = rate_total_S, 
                         mu = mu, 
                         k = k, 
                         x = x)
    EIP <- unlist(EIP)
    out_i <- data.frame("EIP_50" = quantile(EIP, probs = c(0.5))[[1]],
             "EIP_10" = quantile(EIP, probs = c(0.1))[[1]],
             "EIP_90" = quantile(EIP, probs = c(0.9))[[1]],
             "mean" = mean(EIP),
             "temp" = temps[j])
    rm(list = c("EIP"))
    return(out_i)},
    
    n_iter = n_iter,
    temp = temps,
    shape_total_S = shape_total_S, 
    rate_total_S = rate_total_S, 
    mu = mu, 
    k = k, 
    x = x
  ))
  
 return(out)   
}

EIP_plot_g <- calc_EIP_v(EIP_index = EIP_index_g, )
EIP_plot_s <- calc_EIP_v(EIP_index = EIP_index_s, temps = temps_s, n_iter = 12000)
saveRDS(EIP_plot_g, file = "EIP_gambiae.rds")      
saveRDS(EIP_plot_s, file = "EIP_stephensi.rds")      

EIP_plot_g <- readRDS(file = "EIP_gambiae.rds")
EIP_plot_s <- readRDS(file = "EIP_stephensi.rds")

EIP_plot <- ggplot() +
  geom_ribbon(data = rbind(EIP_plot_g %>% mutate(vector_species = "gambiae"),
                    EIP_plot_s %>% mutate(vector_species = "stephensi")),
       aes(x = temp, ymin = EIP_10, ymax = EIP_90, fill = vector_species), 
       alpha = 0.325) +
  geom_line(data = rbind(EIP_plot_g %>% mutate(vector_species = "gambiae"),
                         EIP_plot_s %>% mutate(vector_species = "stephensi")),
            aes(x = temp, y = EIP_50, colour = vector_species), size = 1.25) +
  geom_line(data = EIP_plot_g, aes(x = temp, y = mean), col = "#56B4E9", linetype = 2, size = 1.25) +
  geom_line(data = EIP_plot_s, aes(x = temp, y = mean), col = "#E69F00", linetype = 2, size = 1.25) +
  theme_bw() + theme(text = element_text(size = 15)) +
  ylab("EIP") +
  xlab("Temperature (°C)") +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"), labels = c(parse(text="italic('A. gambiae')"), parse(text="italic('A. stephensi')")), name = "") +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), labels = c(parse(text="italic('A. gambiae')"), parse(text="italic('A. stephensi')")), name = "") +
  scale_y_continuous(breaks = seq(0, 80, 10)) +
  scale_x_continuous(breaks = seq(18, 34, 2)) 

############################
### visualising the fits ###
############################

# gambiae data
g_s_data <- read.csv(file = "vector_fits/data/ES_new_constant_temp_spz_processed.csv") %>% mutate(vector_species = "gambiae",
                                                                                                        gametocytemia = round(gametocytemia, digits = 5),
                                                                                                        ref = "Eunho_Suh")
g_s_data <- g_s_data[,c(2:ncol(g_s_data))]

g_s_data$gametocytemia <- round(g_s_data$gametocytemia, digits = 5)
g_s_data <- subset(g_s_data, gametocytemia!=0.00024)

g_s_data <- g_s_data %>% mutate(index_temp = match(g_s_data$temp, sort(unique(g_s_data$temp))))

g_s_data_in <- generate_prevalence_temp(g_s_data)

# posterior predictive distribution times
PPD_times_g <- readRDS(file = "data/PPD_times_g")
length_ppd_times <- length(PPD_times_g)

iterations <- 5500
warmup <- 3000
chains <- 4

g_fits <- run_prop_ppd_df(fit_gambiae, "pooled", "S_prevalence_ppd", length_ppd_times, PPD_times_g, unique_temp = sort(unique(g_s_data_in$temp)))

# stephensi data
stephensi_data <- read.csv("vector_fits/data/mSOS_all_parasite_data.csv") %>% subset(DTR == 0) %>% 
  mutate(gametocytemia = ifelse(ref == "Shapiro_L_Thomas" | ref == "Murdock_Thomas" | ref == "Shapiro_Thomas", 0.08, NA))

stephensi_s_data <- subset(stephensi_data, lifestage == "gland-spz")
stephensi_s_data <- subset(stephensi_s_data, vector_species!="gambiae")

s_fits <- prop_ppd_function(as.data.frame(fit_stephensi), 
                            n_unq_gt = length(sort(unique(stephensi_data$temp))), 
                            length_ppd_times = length(seq(0,30,0.1)),
                            PPD_times = seq(0, 30, 0.1),
                            iterations = 4500, 
                            warmup = 1500, 
                            chains = 4,
                            Stan_data_name = "sporozoite_prevalence_ppd")

s_fits <- s_fits %>% mutate(temp = sort(unique(stephensi_data$temp))[index_gt])

fits_df <- rbind(g_fits %>% mutate("vector_species" = "gambiae"),
                 s_fits %>% mutate("model" = "pooled",
                   "vector_species" = "stephensi"))

stephensi_s_data <- stephensi_s_data[,c("day_post_inf", "presence", "temp", "gametocytemia", "vector_species", "ref")]
colnames(stephensi_s_data)[1] <- "DPI"

s_data_in <- stephensi_s_data %>% group_by(DPI, temp, vector_species) %>% summarise(positive = sum(presence),
                                                                                         sample = n()) %>% 
  ungroup() %>% 
  mutate(prevalence = positive/sample,
         lower = prevalence - (1.96 * sqrt(prevalence * (1 - prevalence) / sample)),
         upper = prevalence + (1.96 * sqrt(prevalence * (1 - prevalence) / sample)))

s_data <- rbind(g_s_data_in[,-which(colnames(g_s_data_in)=="index_temp")] %>% mutate(vector_species = "gambiae"), s_data_in)

m_days <- as.data.frame(s_data %>% group_by(temp) %>% summarise(m_dpi = max(DPI)))

fits_plot_df <- bind_rows(lapply(seq(1, nrow(m_days)), function(i, fits_df, m_days){
   subset(fits_df, DPI < (m_days[i, "m_dpi"] + 5) & temp == m_days[i, "temp"])
},
fits_df = fits_df,
m_days = m_days))

fit_plot <- ggplot(data = fits_plot_df %>% 
                     mutate(temp_label = paste0(temp, "°C")), aes(x = DPI, y = median, ymin = lower, ymax = upper, group = vector_species)) +
  geom_pointrange(data = s_data %>% 
                    mutate(temp_label = paste0(temp, "°C")), 
                  aes(x = DPI, y = prevalence, ymin = lower, ymax = upper, 
                      fill = vector_species), alpha = 0.75, shape = 21, col = "grey50") +
  geom_ribbon(alpha = 0.325, aes(fill = vector_species)) +
  geom_line(aes(col = vector_species), size = 1) +
  theme_bw() + theme(text = element_text(size = 15)) +
  facet_wrap(~temp_label) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"), labels = c(parse(text="italic('A. gambiae')"), parse(text="italic('A. stephensi')")), name = "") +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), labels = c(parse(text="italic('A. gambiae')"), parse(text="italic('A. stephensi')")), name = "") +
  ylab("Sporozoite prevalence") +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 0.75, 0.25)) +
  xlab("Days post infection")

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

png(file = "results/figures/vector_plot.png", height = 770, width = 770)
plot_grid(
  fit_plot,
  plot_grid(
    EIP_plot + theme(legend.position = "none"),
    pr_g_s_plot + theme(legend.position = "none"), 
    ncol = 2, labels = c("B", "C")
  ),
  nrow = 2,
  labels = c("A", ""),
  rel_heights = c(1, 0.7)
)
  
dev.off()

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

