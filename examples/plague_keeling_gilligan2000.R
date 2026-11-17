## Bubonic plague: rats (SIR with inherited resistance), flea index, free infectious fleas, human SIR (Keeling & Gilligan 2000) [STRUCTURE-FROM-RECALL]
## Bubonic plague (Yersinia pestis), Rattus rattus – Xenopsylla cheopis – humans; historical Europe.
## Parameter values: 0 of 15 checked against the source; the rest are recalled or illustrative.
## Reference: Keeling MJ, Gilligan CA (2000). Metapopulation dynamics of bubonic plague. Nature 407:903–906. https://doi.org/10.1038/35038073 ; and Keeling MJ, Gilligan CA (2000) Bubonic plague: a metapopulation model of a zoonosis. Proc R Soc B 267:2219–2230,…
## Source: https://doi.org/10.1038/35038073
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Srat, Irat, Rrat, Nf, Ff, Shum, Ihum, Rhum
## Flows (time unit: years):
##   (outside) -> Srat  : r_R * Srat
##   (outside) -> Srat  : r_R * (1 - p_inh) * Rrat
##   (outside) -> Rrat  : r_R * p_inh * Rrat
##   Srat -> (outside)  : (d_R + r_R * Nrat / K_R) * Srat
##   Rrat -> (outside)  : (d_R + r_R * Nrat / K_R) * Rrat
##   Srat -> Irat       : lam_rat * Srat
##   Irat -> Rrat       : m_R * g_R * Irat
##   Irat -> (outside)  : m_R * (1 - g_R) * Irat
##   Irat -> (outside)  : d_R * Irat
##   (outside) -> Nf    : r_F * Nf
##   Nf -> (outside)    : r_F * Nf / K_F * Nf
##   (outside) -> Ff    : (d_R + m_R * (1 - g_R)) * Irat * Nf
##   Ff -> (outside)    : d_F * Ff
##   (outside) -> Shum  : mu_hum * Nhum
##   Shum -> Ihum       : lam_hum * Shum
##   Ihum -> Rhum       : m_H * g_H * Ihum
##   Ihum -> (outside)  : m_H * (1 - g_H) * Ihum
##   Shum -> (outside)  : mu_hum * Shum
##   Ihum -> (outside)  : mu_hum * Ihum
##   Rhum -> (outside)  : mu_hum * Rhum

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 50, by = 0.05)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  r_R = 5,        # rat birth rate (1/yr)
  p_inh = 0.975,  # probability offspring of resistant rats inherit resistance
  K_R = 2500,     # rat carrying capacity (rats)
  d_R = 0.2,      # rat natural death (1/yr)
  bet_R = 4.7,    # rat infection rate per free flea (1/yr)
  a_s = 0.004,    # flea searching efficiency (1/rat)
  m_R = 20,       # rat infectious-period rate (≈18 d) (1/yr)
  g_R = 0.02,     # probability an infected rat survives and becomes resistant
  r_F = 20,       # flea growth rate on rats (1/yr)
  K_F = 6.57,     # flea carrying capacity per rat (fleas/rat)
  d_F = 10,       # death rate of free fleas (1/yr)
  bet_H = 0.01,   # value taken from the notes text [recalled]
  mu_hum = 0.02,  # value taken from the notes text [recalled]
  m_H = 26,       # value taken from the notes text [recalled]
  g_H = 0.4       # value taken from the notes text [recalled]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Srat = 2499,   # susceptible rats
  Irat = 1,      # infected rats
  Rrat = 0,      # resistant rats
  Nf = 6.57,     # mean fleas per rat (flea index; not a head count)
  Ff = 0,        # free infectious fleas (searching)
  Shum = 10000,  # susceptible humans
  Ihum = 0,      # infected humans
  Rhum = 0       # recovered humans
)

# ---- Model equations -------------------------------------------------
plague_keeling_gilligan2000_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Shum + Ihum + Rhum

    # Derived quantities
    Nrat <- Srat + Irat + Rrat
    Nhum <- Shum + Ihum + Rhum
    lam_rat <- bet_R * Ff * (1 - exp(-a_s * Nrat)) / Nrat  # rat force of infection: fleas that find rats infect them
    fleas_to_humans <- Ff * exp(-a_s * Nrat)               # free fleas that fail to find a rat and land on humans (K&G's human force of infection)
    lam_hum <- bet_H * fleas_to_humans / Nhum              # human force of infection

    # Rates of change
    dSrat <- r_R * Srat + r_R * (1 - p_inh) * Rrat - (d_R + r_R * Nrat / K_R) * Srat - lam_rat * Srat
    dIrat <- lam_rat * Srat - m_R * g_R * Irat - m_R * (1 - g_R) * Irat - d_R * Irat
    dRrat <- r_R * p_inh * Rrat - (d_R + r_R * Nrat / K_R) * Rrat + m_R * g_R * Irat
    dNf <- r_F * Nf - r_F * Nf / K_F * Nf
    dFf <- (d_R + m_R * (1 - g_R)) * Irat * Nf - d_F * Ff
    dShum <- mu_hum * Nhum - lam_hum * Shum - mu_hum * Shum
    dIhum <- lam_hum * Shum - m_H * g_H * Ihum - m_H * (1 - g_H) * Ihum - mu_hum * Ihum
    dRhum <- m_H * g_H * Ihum - mu_hum * Rhum

    list(c(dSrat, dIrat, dRrat, dNf, dFf, dShum, dIhum, dRhum))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = plague_keeling_gilligan2000_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Srat = "#333333", Irat = "#B71C1C", Rrat = "#4E8A3E", Nf = "#333333", Ff = "#5D4037", Shum = "#333333", Ihum = "#B71C1C", Rhum = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Bubonic plague: rats (SIR with inherited resistance), flea index, free infectious fleas, human SIR (Keeling & Gilligan 2000) [STRUCTURE-FROM-RECALL]")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
