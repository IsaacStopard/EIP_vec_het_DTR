rm(list = ls())
library(tidyverse)

bt <- read.csv(file = "data/BF_biting_time.csv", header = FALSE)

colnames(bt) <- c("hour", "bites")

bt$hour <- round(bt$hour, digits = 0)
bt$bites <- round(bt$bites, digits = 2)

bt <- bt %>% group_by(hour) %>% mutate(location = ifelse(bites == min(bites), "indoor", "outdoor")) %>% ungroup()

fit <- lm(formula = bites~poly(hour, 2)+location, data = bt)

ggplot(data = bt, aes(x = hour, y = bites, col = location)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE, fullrange = TRUE) +
  xlim(19, (19 + 23)) +
  ylim(0, 4)
  
# https://www.pnas.org/doi/10.1073/pnas.2104282119

