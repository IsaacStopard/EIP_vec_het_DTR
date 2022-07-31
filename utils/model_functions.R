# Model functions for the sensitivity analysis of the impact of EIP (mean and variance) 
# and seasonality on sporozoite prevalence
# Author: Isaac J Stopard
# Version: 0.01 
# Last updated: 14/06/2021
# Notes: 

####################################
### classic model at equilibrium ###
####################################
### equilibrium of classic model
eq_I_m <- function(a, c, I_h, mu, n){
  return((a * c * I_h * exp(- n * mu))/ (a * c * I_h + mu)) 
}

v.eq_I_m <- Vectorize(eq_I_m)

eq_I_m_full <- function(a, b, c, mu, m, n, r){
  return((-r * mu + a^2 * b * c * m * exp(-n * mu)) / (a * b * m * (a * c + mu)))
}

v.eq_I_m_full <- Vectorize(eq_I_m_full)

eq_I_h <- function(a, b, c, mu, m, n, r){
  p <- exp(- n * mu)
  return((-r * mu + a^2 * b * c * m * p)/(a * c * (r + a * b * m * p)))
}

v.eq_I_h <- Vectorize(eq_I_h)

R0_M1 <- function(m, a, b, c, mu, n, r){
  return((m * a^2 * b * c * exp(-mu * n))/(r * mu))
}

v.R0_M1 <- Vectorize(R0_M1)

### discretized EIP PDF
discrete_eq_I_m <- function(a, c, I_h, mu, n, p){
  return(sum(exp(-n * mu) * p) * a * c * I_h /
           (a * c * I_h + mu))
}

discrete_eq_I_m_full <- function(a, b, c, mu, m, r, n, p){
  return((-r * mu + a^2 * b * c * m * sum(exp(-n * mu) * p)) / (a * b * m * (a * c + mu)))
}

single_discrete_eq_I_m <- function(a, c, I_h, mu, n, p){
  return(exp(-n * mu) * p * a * c * I_h /
           (a * c * I_h + mu))
}

single_discrete_eq_I_m_full <- function(a, b, c, mu, m, r, n, p, n_i, p_i){
  return((exp(-n_i * mu)*p_i*(-r * mu + a^2 * b * c * m * sum(exp(-n * mu) * p))) / (a * b * m * (a * c + mu) * sum(exp(-n * mu) * p))) 
}

v.single_discrete_eq_I_m_full <- Vectorize(single_discrete_eq_I_m_full)

# continuous PDF model
prop_surv <- function(t, mu, shape, rate, mu_PL, k){
  return(exp(-mu * t) * EIP_PDF(t, shape, rate, mu_PL, k))
}

continuous_eq_I_m <- function(a, c, I_h, mu, shape, rate, mu_PL, k){
  p <- integrate(prop_surv, lower = 0, upper = Inf, rel.tol =.Machine$double.eps^0.5, mu = mu, shape = shape, rate = rate, mu_PL = mu_PL, k = k)[[1]]
  return(a * c * I_h * p / (a * c * I_h + mu))
}

continuous_eq_I_m_full <- function(a, b, c, mu, m, r, shape, rate, mu_PL, k){
  p <- integrate(prop_surv, lower = 0, upper = Inf, rel.tol =.Machine$double.eps^0.5, mu = mu, shape = shape, rate = rate, mu_PL = mu_PL, k = k)[[1]]
  return((-r * mu + a^2 * b * c * m * p)/ (a * b * m * (a * c + mu)))
}

v.continuous_eq_I_m_full <- Vectorize(continuous_eq_I_m_full)

prop_surv_int <- function(mu, shape, rate, mu_PL, k){
  return(integrate(prop_surv, lower = 0, upper = Inf, mu = mu, shape = shape, rate = rate, mu_PL = mu_PL, k = k)[[1]])
}

v.prop_surv_int <- Vectorize(prop_surv_int)

continuous_eq_I_h_full <- function(a, b, c, mu, m, r, shape, rate, mu_PL, k){
  p <- integrate(prop_surv, lower = 0, upper = Inf, rel.tol =.Machine$double.eps^0.5, mu = mu, shape = shape, rate = rate, mu_PL = mu_PL, k = k)[[1]]
  return((-r * mu + a^2 * b * c * m * p)/(a * c * (r + a * b * m * p)))
}

v.continuous_eq_I_h_full <- Vectorize(continuous_eq_I_h_full)

R0 <- function(a, b, c, mu, m, r, n){
  return((m * a^2 * b * c * exp(-mu * n))/(r * mu))
}

v.R0 <- Vectorize(R0)

continuous_R0 <- function(a, b, c, mu, m, r, n, shape, rate, mu_PL, k){
  p <- integrate(prop_surv, lower = 0, upper = Inf, rel.tol =.Machine$double.eps^0.5, mu = mu, shape = shape, rate = rate, mu_PL = mu_PL, k = k)[[1]]
  return((m * a^2 * b * c * p)/(r * mu))
}
  
v.continuous_R0 <- Vectorize(continuous_R0)

discrete_R0 <- function(a, b, c, mu, m, r, n, p){
  return((m * a^2 * b * c * sum(exp(-n * mu) * p))/(r * mu))
}

v.discrete_R0 <- Vectorize(discrete_R0)
######################
### dynamic models ###
######################

# desolve models

# DDE model with a single EIP value
# no births
single_model_nb <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    # EIP delay
    EIP_lag <- t - n_1
    if(EIP_lag < 0){ # this needs to be the time of the first
      state_lag <- rep(0,length(state)) # the initial values
    } else{
      state_lag <- lagvalue(EIP_lag)
    }
    
    ############
    ### DDEs ###
    ############
    
    M <- S + E + I
    
    dS <- -(a * c * I_h * S)  - (mu * S)
    
    dE <- (a * c * I_h * S) - (a * c * I_h * state_lag[1] * exp(-mu * n_1)) - (mu * E)
    
    dI <- (a * c * I_h * state_lag[1] * exp(-mu * n_1)) - (mu * I)
    
    list(c(dS, dE, dI))
  })
}

# DDE model with single EIP value
# constant mosquito births
single_model_b <- function(t, state, params, init_vals = state_single){
  with(as.list(c(state, params)),{
    EIP_lag <- t - n_1
    if(EIP_lag < 0){
      state_lag <- rep(0, length(state_vals)) # the initial values
    } else{
      state_lag <- lagvalue(EIP_lag)
    }
    
    dE_m <- (a * c * I_h * (1 - E_m - I_m)) - 
      (a * c * I_h * (1 - state_lag[1] - state_lag[2]) * exp(-mu * n_1)) -
      (mu * E_m)
    
    dI_m <- (a * c * I_h * (1 - state_lag[1] - state_lag[2]) * exp(-mu * n_1)) - (mu * I_m)
    
    list(c(dE_m, dI_m))
  })
}

# ODE model
# Erlang distribution model
# no mosquito births
multi_comp_model_nb <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    M <- S + E_ + I
    
    dS <- - (a * c * I_h * S)  - (mu * S)
    
    dE1 <-  (a * c * I_h * S) - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dS, E_out, dI))
  })
}

# Erlang distribution EIP model with constant births
multi_comp_model_cb <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    M <- S + E_ + I
    
    dI_h <- a * b * m * I * (1 - I_h) - r * I_h
    
    dS <- M * (beta + mu) - (a * c * I_h * S)  - (mu * S)
    
    dE1 <-  (a * c * I_h * S) - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dS, E_out, dI, dI_h))
  })
}

# model 3 - Erlang distribution model with time-varying EIP and constant births
model_3 <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    rate <- shape/EIP_fun_t(t)
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    M <- S + E_ + I
    
    dS <- M * (beta + mu) - (a * c * I_h * S)  - (mu * S)
    
    dE1 <-  (a * c * I_h * S) - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    dI_h <- a * b * m * I * (1 - I_h) - r * I_h
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dS, E_out, dI, dI_h))
  })
}

multi_comp_model_b <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    M <- S + E_ + I
    
    beta <- beta_fun(t)
    
    dS <- M * (beta + mu) - (a * c * I_h * S)  - (mu * S)
    
    dE1 <-  (a * c * I_h * S) - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dS, E_out, dI))
  })
}

multi_comp_model_eb <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    M <- S + E_ + I
    
    beta <- beta_fun(t)
    
    dS <- beta - (a * c * I_h * S)  - (mu * S)
    
    dE1 <-  (a * c * I_h * S) - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dS, E_out, dI))
  })
}

# model 4
# SMFA model: no births
# needs to be on the log scale so the numbers aren't below zero
model_4 <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    if(v_EIP == TRUE){
      rate <- shape/EIP_fun[[EIP_CI]](temp_fun[[index]](t))
    }
    
    #rate <- shape/EIP_fun(temp_fun[[i]](t), a, b, c, t0)
    dU <- - mu * U
    
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    dE1 <- - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dU, E_out, dI))
  })
}

# model with varying 
model_4_ht <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    
    rate <- shape/EIP_fun_l[[index_EIP]](temp_fun_l[[index_temp]](t))
    #rate <- shape/EIP_fun(temp_fun[[i]](t), a, b, c, t0)
    dU <- - mu * U
    
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    dE1 <- - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dU, E_out, dI))
  })
}


multi_comp_model <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    M <- S + E_ + I
    
    dS <- - (a * c * I_h * S)  - (mu * S)
    
    dE1 <-  (a * c * I_h * S) - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dS, E_out, dI))
  })
}



# extending to the multiple compartment version
multiple_model <- function(t, state, params){
  with(as.list(c(state, params)),{
    
    EIP_lag <- rep(NA, length_n)
    state_lag <- matrix(NA, nrow = length_n + 1, ncol = length_n + 1)
    for(i in 1:length_n){
      EIP_lag[i] <- t - get(paste0("n",i))
      if(EIP_lag[i] < 0){
        state_lag[i+1,] <- init_vals
      } else{
        state_lag[i+1,] <- lagvalue(EIP_lag[i]) # values of the state variables at previous times
      }
    }
    if(t > 0){
      state_lag[1,] <- lagvalue(t) # the current time
    } else{
      state_lag[1,] <- init_vals
    }
    
    for(i in 1:length_n){
      assign(paste0("dE_m",i),
             (a * c * I_h * (1 - sum(state_lag[1,c(1:(length_n + 1))])) * get(paste0("p",i))) -
               (a * c * I_h * (1 - sum(state_lag[i+1, c(1:(length_n + 1))])) * get(paste0("p",i)) * exp(-mu * get(paste0("n",i)))) -
               (mu * get(paste0("E_m",i))))
    }
    
    dI_m <- (a * c * I_h * (1 - sum(state_lag[2, c(1:(length_n + 1))])) * get(paste0("p",1)) * exp(-mu * get(paste0("n",1))))
    for(i in 2:length_n){
      dI_m <- dI_m + 
        (a * c * I_h * (1 - sum(state_lag[i+1, c(1:(length_n + 1))])) * get(paste0("p",i)) * exp(-mu * get(paste0("n",i))))
    }
    dI_m <- dI_m - (mu * I_m)
    
    # creating the list to return
    out <- rep(NA, length_n + 1)
    for(i in 1:length_n){
      out[i] <- get(paste0("dE_m",i))
    }
    out[length_n + 1] <- dI_m
    list(out)
  })
}

#########################################
### IC malaria model helper functions ###
#########################################

# function that creates the list of parameters required for the model with a single mosquito species or multiple species
# currently infectiousness of treated patients differs - this is updated so it doesn't before running the model
# comp - number of exposed mosquito compartments - initially set to 10
# EIP is initially set to 10 days
model_param_list_create_comp <- function(
  # age, heterogeneity in exposure,
  eta = 1/(21*365),
  rho = 0.85,
  a0 = 2920,
  sigma2 = 1.67,
  max_age = 100*365,
  #  rate of leaving infection states
  rA = 1/195,
  rT = 0.2,
  rD = 0.2,
  rU = 1/110.299,
  rP = 1/15,
  #  human latent period and time lag from asexual parasites to
  dE  = 12,
  delayGam = 12.5,
  # human infectiousness to mosquitoes
  cD  = 0.0676909,
  cT  =  0.322 * cD,
  cU  = 0.006203,
  gamma1  = 1.82425,
  #  Immunity reducing probability of detection
  d1 = 0.160527,
  dID = 3650,
  ID0 = 1.577533,
  kD = 0.476614,
  uD = 9.44512,
  aD = 8001.99,
  fD0 = 0.007055,
  gammaD = 4.8183,
  alphaA = 0.75735,
  alphaU = 0.185624,
  # Immunity reducing probability of infection
  b0 = 0.590076,
  b1 = 0.5,
  dB = 3650,
  IB0 = 43.8787,
  kB = 2.15506,
  uB = 7.19919,
  # Immunity reducing probability of clinical disease
  phi0 = 0.791666,
  phi1 = 0.000737,
  dCA = 10950,
  IC0 = 18.02366,
  kC = 2.36949,
  uCA = 6.06349,
  PM = 0.774368,
  dCM = 67.6952,
  # entomological parameters
  #delayMos = 10,
  comp = 10,
  EIP = 10,
  tau1 = 0.69,
  tau2 = 2.31,
  mu0 = 0.132,
  Q0 = 0.92,
  chi = 0.86,
  bites_Bed = 0.89,
  bites_Indoors = 0.97,
  # intervention parameters
  num_int = 1,
  itn_cov = 0,
  irs_cov = 0,
  ITN_IRS_on = -1,
  DY = 365,
  d_ITN0 = 0.41,
  r_ITN0 = 0.56,
  r_ITN1 = 0.24,
  r_IRS0 = 0.6,
  d_IRS0 = 1,
  irs_half_life =   0.5 * DY,
  itn_half_life =   2.64 * DY,
  IRS_interval =   1 * DY,
  ITN_interval =   3 * DY,
  ...
  
){
  # set up param list
  mp_list <- list()
  
  # catach extra params and place in list
  extra_param_list <- list(...)
  if(length(extra_param_list)>0){
    if(is.list(extra_param_list[[1]])){
      extra_param_list <- extra_param_list[[1]]
    }
  }
  
  ## DEFAULT PARAMS
  
  # duration of year
  mp_list$DY <- DY
  
  # age, heterogeneity in exposure
  mp_list$eta <- eta
  mp_list$rho <- rho
  mp_list$a0 <- a0
  mp_list$sigma2 <- sigma2
  mp_list$max_age <- max_age
  
  # rate of leaving infection states
  mp_list$rA <- rA
  mp_list$rT <- rT
  mp_list$rD <- rD
  mp_list$rU <- rU
  mp_list$rP <- rP
  
  # human latent period and time lag from asexual parasites to
  # infectiousness
  mp_list$dE <- dE
  mp_list$delayGam <- delayGam
  
  # infectiousness to mosquitoes
  mp_list$cD <- cD
  mp_list$cT <- cT
  mp_list$cU <- cU
  mp_list$gamma1 <- gamma1
  
  # Immunity reducing probability of detection
  mp_list$d1 <- d1
  mp_list$dID <- dID
  mp_list$ID0 <- ID0
  mp_list$kD <- kD
  mp_list$uD <- uD
  mp_list$aD <- aD
  mp_list$fD0 <- fD0
  mp_list$gammaD <- gammaD
  
  # PCR prevalence parameters
  mp_list$alphaA <- alphaA
  mp_list$alphaU <- alphaU
  
  # anti-infection immunity
  mp_list$b0 <- b0
  mp_list$b1 <- b1
  mp_list$dB <- dB
  mp_list$IB0 <- IB0
  mp_list$kB <- kB
  mp_list$uB <- uB
  
  # clinical immunity
  mp_list$phi0 <- phi0
  mp_list$phi1 <- phi1
  mp_list$dCA <- dCA
  mp_list$IC0 <- IC0
  mp_list$kC <- kC
  mp_list$uCA <- uCA
  mp_list$PM <- PM
  mp_list$dCM <- dCM
  
  # entomological parameters
  #mp_list$delayMos <- delayMos
  mp_list$tau1 <- tau1
  mp_list$tau2 <- tau2
  mp_list$mu0 <- mu0
  mp_list$Q0 <- Q0
  mp_list$chi <- chi
  mp_list$bites_Bed <- bites_Bed
  mp_list$bites_Indoors <- bites_Indoors
  mp_list$fv0 <- 1 / (tau1 + tau2)
  mp_list$av0 <- Q0 * mp_list$fv0 # daily feeeding rate on humans
  #mp_list$Surv0 <- exp(-mu0 * delayMos) # probability of surviving incubation period
  mp_list$p10 <- exp(-mu0 * tau1)  # probability of surviving one feeding cycle
  mp_list$p2 <- exp(-mu0 * tau2)  # probability of surviving one resting cycle
  
  # Erlang distribution EIP
  mp_list$comp <- comp
  mp_list$delayMos <- comp/EIP
  
  # ITN/IRS parameters
  mp_list$itn_cov <- itn_cov
  mp_list$irs_cov <- irs_cov
  
  mp_list$num_int <- num_int
  # Catch all: Not defined the correct number of interventions
  if (itn_cov > 0 & num_int == 1){
    stop(message("Incorrect number of interventions for definied ITN coverage. Please ensure you have correctly
                 specified the number of interventions."))
  }
  if (irs_cov > 0 & num_int < 3){
    stop(message("Incorrect number of interventions for definied IRS coverage. Please ensure you have correctly
                 specified the number of interventions."))
  }
  
  # Sets start time of coverage
  mp_list$ITN_IRS_on <- ITN_IRS_on
  
  # Sets population split as coverage
  # {No intervention} {ITN only} {IRS only} {Both ITN and IRS}
  cov <- c((1 - itn_cov) * (1 - irs_cov), itn_cov * (1 - irs_cov), (1 - itn_cov) * irs_cov, itn_cov * irs_cov)
  cov <- cov[1:mp_list$num_int]
  mp_list$cov <- cov
  
  mp_list$d_ITN0 <- d_ITN0
  mp_list$r_ITN0 <- r_ITN0
  mp_list$r_ITN1 <- r_ITN1
  mp_list$r_IRS0 <- r_IRS0
  mp_list$d_IRS0 <- d_IRS0
  mp_list$irs_half_life <- irs_half_life
  mp_list$itn_half_life <- itn_half_life
  mp_list$IRS_interval <- IRS_interval
  mp_list$ITN_interval <- ITN_interval
  mp_list$irs_loss <- log(2)/mp_list$irs_half_life
  mp_list$itn_loss <- log(2)/mp_list$itn_half_life
  
  # check that none of the spare parameters in the extra
  if(sum(!is.na(match(names(extra_param_list),names(mp_list))))!=0){
    
    stop (message(cat("Extra params in ... share names with default param names. Please check:\n",
                      names(extra_param_list)[!is.na(match(names(extra_param_list),names(mp_list)))]
    )
    ))
  }
  
  return(append(mp_list,extra_param_list))
}

# function that generates the starting values - still need to run for a warmup period before equilibrium is reached
# sets the initial mosquito infection to 0.01
equilibrium_init_create_comp <- function(age_vector, 
                                         het_brackets,
                                         ft = 1,
                                         model_param_list,
                                         mv0_in,
                                         EIR_in,
                                         mf,
                                         multi_species = FALSE)
{
  
  # mpl is shorter :)
  mpl <- model_param_list
  
  ## Check Parameters
  if(!is.numeric(age_vector)) stop("age_vector provided is not numeric")
  if(!is.numeric(het_brackets)) stop("het_brackets provided is not numeric")
  if(!is.numeric(ft)) stop("ft provided is not numeric")
  
  ## population demographics
  age <- age_vector * mpl$DY
  na <- as.integer(length(age))  # number of age groups
  nh <- as.integer(het_brackets)  # number of heterogeneity groups
  h <- statmod::gauss.quad.prob(nh, dist = "normal")
  age0 <- 2
  age1 <- 10
  num_int <- mpl$num_int
  
  age_rate <- age_width <- age_mid_point <- den <- c()
  for (i in 1:(na-1)){
    age_width[i] <- age[i+1] - age[i]
    age_rate[i] <- 1/(age[i + 1] - age[i])  # vector of rates at which people leave each age group (1/age group width)
    age_mid_point[i] <- 0.5 * (age[i] + age[i + 1])  # set age group vector to the midpoint of the group
  }
  age_rate[na] = 0
  
  den <- 1/(1 + age_rate[1]/mpl$eta)
  for (i in 1:(na-1)){
    den[i+1] <- age_rate[i] * den[i]/(age_rate[i+1] + mpl$eta)  # proportion in each age_vector group
  }
  
  age59 <- which(age_vector * 12 > 59)[1] - 1  # index of age vector before age is >59 months
  age05 <- which(age_vector > 5)[1] - 1  # index of age vector before age is 5 years
  
  ## force of infection
  foi_age <- c()
  for (i in 1:na){
    foi_age[i] <- 1 - (mpl$rho * exp(-age[i]/mpl$a0))  #force of infection for each age group
  }
  fden <- foi_age * den
  omega <- sum(fden)  #normalising constant
  
  ## heterogeneity
  het_x <- h$nodes
  het_wt <- h$weights
  den_het <- outer(den, het_wt)
  rel_foi <- exp(-mpl$sigma2/2 + sqrt(mpl$sigma2) * het_x)/sum(het_wt * exp(-mpl$sigma2/2 + sqrt(mpl$sigma2) * het_x))
  
  ## EIR
  EIRY_eq <- EIR_in  # initial annual EIR
  EIRd_eq <- EIRY_eq/mpl$DY
  EIR_eq <- outer(foi_age, rel_foi) * EIRd_eq
  
  ## Immunity and FOI
  x_I <- den[1]/mpl$eta
  for (i in 2:na){
    x_I[i] <- den[i]/(den[i - 1] * age_rate[i - 1])  #temporary variables
  }
  fd <- 1 - (1 - mpl$fD0)/(1 + (age/mpl$aD)^mpl$gammaD)
  
  # maternal immunity begins at a level proportional to the clinical
  # immunity of a 20 year old, this code finds that level
  age20i <- rep(0, na)
  for (i in 2:na){
    age20i[i] <- ifelse(age[i] >= (20 * mpl$DY) & age[i - 1] < (20 * mpl$DY),
                        i, age20i[i - 1])
  }
  age20u <- as.integer(age20i[na])
  age20l <- as.integer(age20u - 1)
  age_20_factor <- (20 * mpl$DY - age[age20l] - 0.5 * age_width[age20l]) *
    2/(age_width[age20l] + age_width[age20u])
  
  # finding initial values for all immunity states
  IB_eq <- matrix(0, na, nh)
  FOI_eq <- matrix(0, na, nh)
  ID_eq <- matrix(0, na, nh)
  ICA_eq <- matrix(0, na, nh)
  ICM_init_eq <- vector(length = nh, mode = "numeric")
  ICM_eq <- matrix(0, na, nh)
  cA_eq <- matrix(0, na, nh)
  FOIvij_eq <- matrix(0, na, nh)
  p_det_eq <- matrix(0, na, nh)
  for (j in 1:nh)
  {
    for (i in 1:na){
      IB_eq[i, j] <- (ifelse(i == 1, 0, IB_eq[i - 1, j]) +
                        EIR_eq[i,j]/(EIR_eq[i, j] * mpl$uB + 1) * x_I[i])/(1 + x_I[i]/mpl$dB)
      FOI_eq[i, j] <- EIR_eq[i, j] * ifelse(IB_eq[i, j] == 0, mpl$b0,
                                            mpl$b0 * ((1 - mpl$b1)/(1 + (IB_eq[i, j]/mpl$IB0)^mpl$kB) + mpl$b1))
      ID_eq[i, j] <- (ifelse(i == 1, 0, ID_eq[i - 1, j]) +
                        FOI_eq[i, j]/(FOI_eq[i, j] * mpl$uD + 1) * x_I[i])/(1 + x_I[i]/mpl$dID)
      ICA_eq[i, j] <- (ifelse(i == 1, 0, ICA_eq[i - 1, j]) +
                         FOI_eq[i,j]/(FOI_eq[i, j] * mpl$uCA + 1) * x_I[i])/(1 + x_I[i]/mpl$dCA)
      p_det_eq[i, j] <- mpl$d1 + (1 - mpl$d1)/(1 + fd[i] * (ID_eq[i, j]/mpl$ID0)^mpl$kD)
      cA_eq[i, j] <- mpl$cU + (mpl$cD - mpl$cU) * p_det_eq[i, j]^mpl$gamma1
    }
  }
  # needs to be calculated after because it references ICA
  for (j in 1:nh){
    for (i in 1:na){
      ICM_init_eq[j] <- mpl$PM * (ICA_eq[age20l, j] + age_20_factor *
                                    (ICA_eq[age20u, j] - ICA_eq[age20l, j]))
      ICM_eq[i, j] <- ifelse(i == 1,
                             ICM_init_eq[j], ICM_eq[i - 1,j])/(1 + x_I[i]/mpl$dCM)
    }
  }
  
  IC_eq <- ICM_eq + ICA_eq
  phi_eq <- mpl$phi0 * ((1 - mpl$phi1)/(1 + (IC_eq/mpl$IC0)^mpl$kC) + mpl$phi1)
  
  
  # human states
  gamma <- mpl$eta + c(age_rate[1:(na - 1)], 0)
  delta <- c(mpl$eta, age_rate[1:(na - 1)])
  
  betaT <- matrix(rep(mpl$rT + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  betaD <- matrix(rep(mpl$rD + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  betaP <- matrix(rep(mpl$rP + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  
  aT <- FOI_eq * phi_eq * ft/betaT
  aD <- FOI_eq * phi_eq * (1 - ft)/betaD
  aP <- mpl$rT * aT/betaP
  
  Z_eq <- array(dim = c(na, nh, 4))
  Z_eq[1, , 1] <- den_het[1, ]/(1 + aT[1, ] + aD[1, ] + aP[1, ])
  Z_eq[1, , 2] <- aT[1, ] * Z_eq[1, , 1]
  Z_eq[1, , 3] <- aD[1, ] * Z_eq[1, , 1]
  Z_eq[1, , 4] <- aP[1, ] * Z_eq[1, , 1]
  
  for (j in 1:nh){
    for (i in 2:na){
      Z_eq[i, j, 1] <- (den_het[i, j] - delta[i] * (Z_eq[i - 1, j, 2]/betaT[i, j] +
                                                      Z_eq[i - 1, j, 3]/betaD[i, j] +
                                                      (mpl$rT *  Z_eq[i - 1, j, 2]/betaT[i, j]
                                                       + Z_eq[i - 1, j, 4])/betaP[i, j]))/(1 + aT[i, j] + aD[i, j] + aP[i, j])
      Z_eq[i, j, 2] <- aT[i, j] * Z_eq[i, j, 1] + delta[i] * Z_eq[i -
                                                                    1, j, 2]/betaT[i, j]
      Z_eq[i, j, 3] <- aD[i, j] * Z_eq[i, j, 1] + delta[i] * Z_eq[i -
                                                                    1, j, 3]/betaD[i, j]
      Z_eq[i, j, 4] <- aP[i, j] * Z_eq[i, j, 1] + delta[i] * (mpl$rT *
                                                                Z_eq[i - 1, j, 2]/betaT[i, j] + Z_eq[i - 1, j, 4])/betaP[i,j]
      
    }
  }
  
  Y_eq <- matrix(Z_eq[, , 1], nrow = na, ncol=nh)
  T_eq <- matrix(Z_eq[, , 2], nrow = na, ncol=nh)
  D_eq <- matrix(Z_eq[, , 3], nrow = na, ncol=nh)
  P_eq <- matrix(Z_eq[, , 4], nrow = na, ncol=nh)
  
  betaS <- apply(FOI_eq, MARGIN = 2, FUN = function(x, y){x + y}, y = gamma)
  betaA <- apply(FOI_eq * phi_eq + mpl$rA, MARGIN = 2, FUN = function(x, y){x + y}, y = gamma)
  betaU <- apply(FOI_eq + mpl$rU, MARGIN = 2, FUN = function(x, y){x + y}, y = gamma)
  
  A_eq <- matrix(ncol = nh, nrow = na)
  U_eq <- matrix(ncol = nh, nrow = na)
  S_eq <- matrix(ncol = nh, nrow = na)
  
  for (i in 1:na){
    for (j in 1:nh){
      A_eq[i, j] <- (delta[i] * ifelse(i == 1, 0, A_eq[i - 1, j]) +
                       FOI_eq[i, j] * (1 - phi_eq[i, j]) * Y_eq[i, j] +
                       mpl$rD * D_eq[i,j])/(betaA[i, j] + FOI_eq[i, j] * (1 - phi_eq[i, j]))
      U_eq[i, j] <- (mpl$rA * A_eq[i, j] + delta[i] * ifelse(i == 1,
                                                             0, U_eq[i - 1, j]))/betaU[i, j]
      S_eq[i, j] <- Y_eq[i, j] - A_eq[i, j] - U_eq[i, j]
      FOIvij_eq[i, j] <- foi_age[i] * mpl$av0 * (mpl$cT * T_eq[i, j] + mpl$cD *
                                                   D_eq[i, j] + cA_eq[i, j] * A_eq[i, j] + mpl$cU * U_eq[i, j]) * rel_foi[j]/omega
    }
  }
  
  # mosquito states - # need to update this for the
  FOIv_eq <- sum(FOIvij_eq)
  
  # initial number
  if(multi_species == TRUE){
      Ev_eq_f <- rep(0, mpl$comp)
      Iv_eq_f <- 0.01
      Sv_eq_f <- 1 - 0.01
      Ev_eq_g <- rep(0, mpl$comp)
      Iv_eq_g <- 0.01
      Sv_eq_g<- 1 - 0.01
      mv0_f <- mv0_in[1] * mf
      mv0_g <- mv0_in[2] * mf
  } else{
      Ev_eq <- rep(0, mpl$comp)
      Iv_eq <- 0.01
      Sv_eq <- 1 - 0.01
      #mv0 <- omega * EIRd_eq/(Iv_eq * mpl$av0)
      mv0 <- mv0_in[1] * mf
  }
  #Iv_eq_ <- FOIv_eq/10 * mpl$Surv0/(FOIv_eq + mpl$mu0) # this is not correct - not actually equilibrium
  #Iv_eq <- sum(Iv_eq_)
  #Sv_eq <- sum(mpl$mu0 * Iv_eq_/(FOIv_eq * mpl$Surv0)) # this is not correct - not actually equilibrium
  #Ev_eq <- rep((1 - Sv_eq - Iv_eq)/10, 10) # this is not correct - not actually equilibrium
  
  # mosquito density needed to give this EIR
  
  # add in final dimension - interventions
  num_int <- mpl$num_int
  cov <- mpl$cov
  
  mat <- matrix(0, na, nh)
  
  S_eq <- vapply(cov, FUN = function(x)
  {
    x * S_eq
  }, mat)
  T_eq <- vapply(cov, FUN = function(x)
  {
    x * T_eq
  }, mat)
  D_eq <- vapply(cov, FUN = function(x)
  {
    x * D_eq
  }, mat)
  A_eq <- vapply(cov, FUN = function(x)
  {
    x * A_eq
  }, mat)
  U_eq <- vapply(cov, FUN = function(x)
  {
    x * U_eq
  }, mat)
  P_eq <- vapply(cov, FUN = function(x)
  {
    x * P_eq
  }, mat)
  
  IB_eq = array(IB_eq, c(na, nh, num_int))
  ID_eq = array(ID_eq, c(na, nh, num_int))
  ICA_eq = array(ICA_eq, c(na, nh, num_int))
  ICM_eq = array(ICM_eq, c(na, nh, num_int))
  
  # TODO: Remove this part and put it as an edit to the equilibrium solution
  if(!is.null(mpl$ncc)){
    IB_eq = array(IB_eq, c(na, nh, num_int, mpl$ncc))
    ID_eq = array(ID_eq, c(na, nh, num_int, mpl$ncc))
    ICA_eq = array(ICA_eq, c(na, nh, num_int, mpl$ncc))
    ICM_eq = array(ICM_eq, c(na, nh, num_int, mpl$ncc))
    
    # add in final dimension - interventions
    all_rounds = mpl$MDA_grp_prop*mpl$MDA_cov
    ccov = c(all_rounds, 1-all_rounds)
    
    mat2 <- array(0, c(na,nh, num_int))
    S_eq <- vapply(ccov,FUN = function(x){x * S_eq},mat2)
    T_eq <- vapply(ccov,FUN = function(x){x * T_eq},mat2)
    D_eq <- vapply(ccov,FUN = function(x){x * D_eq},mat2)
    A_eq <- vapply(ccov,FUN = function(x){x * A_eq},mat2)
    U_eq <- vapply(ccov,FUN = function(x){x * U_eq},mat2)
    P_eq <- vapply(ccov,FUN = function(x){x * P_eq},mat2)
  }
  
  # better het bounds for equilbirum initialisation in individual model
  zetas <- rlnorm(n = 1e5,meanlog = -mpl$sigma2/2, sdlog = sqrt(mpl$sigma2))
  while(sum(zetas>100)>0){
    zetas[zetas>100] <- rlnorm(n = sum(zetas>100),meanlog = -mpl$sigma2/2, sdlog = sqrt(mpl$sigma2))
  }
  
  wt_cuts <- round(cumsum(het_wt)*1e5)
  zeros <- which(wt_cuts==0)
  wt_cuts[zeros] <- 1:length(zeros)
  larges <- which(wt_cuts==1e5)
  wt_cuts[larges] <- (1e5 - (length(larges)-1)):1e5
  wt_cuts <- c(0,wt_cuts)
  het_bounds <- sort(zetas)[wt_cuts]
  het_bounds[length(het_bounds)] <- (mpl$max_age/365)+1
  
  ## collate init
  # this is changed as well
  
  if(multi_species == TRUE){
    res <- list(init_S = S_eq, init_T = T_eq, init_D = D_eq, init_A = A_eq, init_U = U_eq,
              init_P = P_eq, init_Y = Y_eq, init_IB = IB_eq, init_ID = ID_eq, init_ICA = ICA_eq,
              init_ICM = ICM_eq, ICM_init_eq = ICM_init_eq,
              init_Iv_f = Iv_eq_f, init_Sv_f = Sv_eq_f, init_Ev_f = Ev_eq_f,
              init_Iv_g = Iv_eq_g, init_Sv_g = Sv_eq_g, init_Ev_g = Ev_eq_g,
              age_width = age_width, age_rate = age_rate, het_wt = het_wt, het_x = het_x,
              omega = omega, foi_age = foi_age, rel_foi = rel_foi,
              mv0_f = mv0_f, mv0_g = mv0_g,
              na = na, nh = nh, ni = num_int, x_I = x_I,
              FOI = FOI_eq, EIR_eq = EIR_eq, cA_eq = cA_eq,
              den = den, age59 = age59, age05 = age05, age = age_vector*mpl$DY, ft = ft, FOIv_eq = FOIv_eq,
              betaS = betaS, betaA = betaA, betaU = betaU, FOIvij_eq=FOIvij_eq,
              age_mid_point = age_mid_point, het_bounds = het_bounds, pi = pi,
              age20l = age20l, age20u = age20u, age_20_factor = age_20_factor)
  } else{
    res <- list(init_S = S_eq, init_T = T_eq, init_D = D_eq, init_A = A_eq, init_U = U_eq,
              init_P = P_eq, init_Y = Y_eq, init_IB = IB_eq, init_ID = ID_eq, init_ICA = ICA_eq,
              init_ICM = ICM_eq, ICM_init_eq = ICM_init_eq, init_Iv = Iv_eq, init_Sv = Sv_eq,
              init_Ev = Ev_eq,
              age_width = age_width, age_rate = age_rate, het_wt = het_wt, het_x = het_x,
              omega = omega, foi_age = foi_age, rel_foi = rel_foi,
              mv0 = mv0, na = na, nh = nh, ni = num_int, x_I = x_I,
              FOI = FOI_eq, EIR_eq = EIR_eq, cA_eq = cA_eq,
              den = den, age59 = age59, age05 = age05, age = age_vector*mpl$DY, ft = ft, FOIv_eq = FOIv_eq,
              betaS = betaS, betaA = betaA, betaU = betaU, FOIvij_eq=FOIvij_eq,
              age_mid_point = age_mid_point, het_bounds = het_bounds, pi = pi,
              age20l = age20l, age20u = age20u, age_20_factor = age_20_factor)
  }

  res <- append(res,mpl)
  
  return(res)
}

# runs the model with set EIP values and a zero beta value
# measure can be mean, min or max - this is the temperature value
run_model <- function(t_df, EIP, measure, n_warm_up, avg_t, EIP_x,
                      state_use, mpl){
  
  # EIP values
  t_df$EIP <- EIP[match(round(t_df[, measure][[1]], digits = 2), 
                        round(EIP$temp, digits = 2)), 
                  "median"]
  
  df_EIP <- data.frame("t" = seq(0, max(t_df$days_) - 1 + n_warm_up, 0.5)) # t is in days
  
  df_EIP$EIP <- t_df[match(df_EIP$t - n_warm_up + 1, t_df$days), "EIP"][[1]]
  
  df_EIP[1:(n_warm_up * 2), "EIP"] <- mean(t_df$EIP) # running with the mean EIP before reaching equilibrium
  
  df_EIP <- df_EIP %>% fill(EIP, .direction = "down")
  
  t <- df_EIP$t
  n_t <- length(t)
  
  state_use$delayMos_t <- t
  state_use$delayMos_rate <- mpl$comp / df_EIP$EIP
  state_use$ldm <- n_t
  
  state_use$beta_t <- t
  state_use$beta_rate <- rep(0, n_t)
  
  mod <- gen(user=state_use, use_dde=TRUE)
  mod_run <- mod$run(t, n_history = 10000) # takes approx - 2 seconds
  out <- mod$transform_variables(mod_run)
  
  out_df <- data.frame("time" = out$t,
                       "m" = out$mv,
                       "I" = out$Iv,
                       "S" = out$Iv / out$mv,
                       "EIP" = out$delayMos,
                       "mu" = out$mu)
  
  # mutate(E = sum(c_across(paste0("E",1):paste0("E",params["shape"]))),
  #         M = S + E + I,
  #         spz = I/M)
  n_out <- nrow(out_df)
  
  out_df$temp_measure <- rep(measure, n_out)
  out_df$avg_t <- rep(avg_t, n_out)
  out_df$EIP_x <- rep(EIP_x, n_out)
  rm(list = c("t_df", "df_EIP"))
  return(out_df)
}

set_c_temp <- function(temp, temp_df, measure){
  temp_df[, measure] <- rep(temp, nrow(temp_df))
  return(temp_df)
}
