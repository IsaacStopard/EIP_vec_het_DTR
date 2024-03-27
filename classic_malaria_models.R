classic_malaria_model <- function(t, state, parms, init_vals = state){
  with(as.list(c(state, parms)),{
    
    temp <- temp_fun(t) # temperature at time t
    if(v_EIP == TRUE){
      rate <- shape/EIP_fun(temp) # PDR at time t
    } else{
      rate <- m_rate
    }
    if(v_HMTP == TRUE){
      c <- delta_fun(t) * scale_HMTP # HMTP at time t
    } else{
      c <- m_HMTP * scale_HMTP
    }
    if(v_mu == TRUE){
      mu <- mu_fun(t)
      beta <- beta_fun_v(t)
    } else{
      mu <- mu_in
      beta <- beta_fun_c(t)
    }
    
    
    
    # calculating the total numbers of exposed mosquitoes
    E_ <- 0
    for(i in 1:shape){
      E_ <- get(paste0("E",i)) + E_
    }
    
    M <- S + E_ + I
    
    z <- I/M
    
    dI_h <- a * b * M * z * (pop - I_h) - r * I_h
    
    x <- I_h / pop
    
    dS <- beta * M - (mu * S) - a * c * S * x
    
    dE1 <- a * c * S * x - (rate * E1) - (mu * E1)
    
    for(i in 2:shape){
      assign(paste0("dE",i), (rate * get(paste0("E",i-1))) - (rate * get(paste0("E",i))) - (mu * get(paste0("E",i))))
    }
    
    dI <- (rate * get(paste0("E",shape))) - (mu * I)
    
    E_out <- c(dE1)
    for(i in 2:shape){
      E_out <- c(E_out, get(paste0("dE", i)))
    }
    
    list(c(dS, E_out, dI, dI_h), temp = temp, rate = rate, M = M, c = c, beta = beta, x = x, z = z, mu = mu)
  })
}
