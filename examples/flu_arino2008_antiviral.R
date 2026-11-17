## Pandemic influenza SLIAR model with antiviral treatment and pre-epidemic vaccination/prophylaxis (Arino et al. 2008)
## Pandemic influenza.
## Parameter values: 0 of 25 checked against the source; the rest are recalled or illustrative.
## Reference: Arino J, Brauer F, van den Driessche P, Watmough J, Wu J. "A model for influenza with vaccination and antiviral treatment." J Theor Biol 2008;253:118–130. https://doi.org/10.1016/j.jtbi.2008.02.026 (PDF:…
## Source: https://doi.org/10.1016/j.jtbi.2008.02.026
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, ST, L, LT, Ip, IpT, A, AT, R, Dead
## Flows (time unit: days):
##   S -> L       : foi * S
##   ST -> LT     : sig_S * foi * ST
##   L -> Ip      : p_sym * k_lat * L
##   L -> A       : (1 - p_sym) * k_lat * L
##   L -> LT      : phi_L * L
##   LT -> L      : theta_L * LT
##   LT -> IpT    : p_sym * tau_av * k_latT * LT
##   LT -> AT     : (1 - p_sym * tau_av) * k_latT * LT
##   Ip -> R      : f_rec * alpha_r * Ip
##   Ip -> Dead   : (1 - f_rec) * alpha_r * Ip
##   Ip -> IpT    : phi_I * Ip
##   IpT -> Ip    : theta_I * IpT
##   IpT -> R     : fT_rec * alpha_T * IpT
##   IpT -> Dead  : (1 - fT_rec) * alpha_T * IpT
##   A -> R       : eta_r * A
##   A -> AT      : phi_A * A
##   AT -> A      : theta_A * AT
##   AT -> R      : eta_T * AT

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 0.2)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  k_lat = 0.526,     # exit rate from L, i.e. 1/1.9 d (1/day)
  k_latT = 0.526,    # exit rate from LT (1/day)
  alpha_r = 0.244,   # exit rate from Ip, i.e. 1/4.1 d (1/day)
  eta_r = 0.244,     # exit rate from A (1/day)
  alpha_T = 0.323,   # exit rate from IpT (1/day)
  eta_T = 0.323,     # exit rate from AT (1/day)
  p_sym = 0.667,     # fraction of latents who become symptomatic
  tau_av = 0.4,      # factor by which antivirals reduce the fraction who develop symptoms
  f_rec = 0.98,      # fraction of Ip who recover (1 - CFR)
  fT_rec = 1.0,      # fraction of IpT who recover
  eps_L = 0,         # relative infectivity of latents; Longini assumes 0
  delta_a = 0.5,     # relative infectivity of asymptomatics
  phi_L = 1.61,      # treatment rate, = -ln(1-0.8), i.e. 80% treated within 1 day (1/day)
  phi_I = 1.61,      # treatment rate, = -ln(1-0.8), i.e. 80% treated within 1 day (1/day)
  phi_A = 1.61,      # treatment rate, = -ln(1-0.8), i.e. 80% treated within 1 day (1/day)
  theta_L = 0.0179,  # relapse rate from treated back to untreated (1/day)
  theta_I = 0.0179,  # relapse rate from treated back to untreated (1/day)
  theta_A = 0.0179,  # relapse rate from treated back to untreated (1/day)
  sig_S = 0.7,       # relative susceptibility of ST
  sig_L = 0.2,       # relative infectivity of treated classes
  sig_I = 0.2,       # relative infectivity of treated classes
  sig_A = 0.2,       # relative infectivity of treated classes
  q_w = 0.4,         # fraction of symptomatics who withdraw from mixing
  b0 = 0.581,        # transmission coefficient under standard incidence with q_w = 0.4 (1/day)
  g_vac = 0          # fraction of S0 treated or vaccinated before the epidemic; the paper explores 0–0.7
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 1,      # untreated susceptible
  ST = 1988,  # treated (prophylaxed or vaccinated) susceptible
  L = 0,      # latent
  LT = 0,     # treated latent
  Ip = 12,    # symptomatic infective
  IpT = 0,    # treated symptomatic infective
  A = 0,      # asymptomatic infective
  AT = 0,     # treated asymptomatic infective
  R = 0,      # recovered
  Dead = 0    # influenza deaths (cumulative)
)

# ---- Model equations -------------------------------------------------
flu_arino2008_antiviral_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + ST + L + LT + Ip + IpT + A + AT + R

    # Derived quantities
    N_M <- N - q_w * Ip - q_w * IpT                                                                                        # mixing population
    Lam <- eps_L * L + eps_L * sig_L * LT + (1 - q_w) * Ip + (1 - q_w) * sig_I * IpT + delta_a * A + delta_a * sig_A * AT  # total infectivity
    foi <- b0 * Lam / N_M                                                                                                  # force of infection on untreated susceptibles

    # Rates of change
    dS <- -foi * S
    dST <- -sig_S * foi * ST
    dL <- foi * S - p_sym * k_lat * L - (1 - p_sym) * k_lat * L - phi_L * L + theta_L * LT
    dLT <- sig_S * foi * ST + phi_L * L - theta_L * LT - p_sym * tau_av * k_latT * LT - (1 - p_sym * tau_av) * k_latT * LT
    dIp <- p_sym * k_lat * L - f_rec * alpha_r * Ip - (1 - f_rec) * alpha_r * Ip - phi_I * Ip + theta_I * IpT
    dIpT <- p_sym * tau_av * k_latT * LT + phi_I * Ip - theta_I * IpT - fT_rec * alpha_T * IpT - (1 - fT_rec) * alpha_T * IpT
    dA <- (1 - p_sym) * k_lat * L - eta_r * A - phi_A * A + theta_A * AT
    dAT <- (1 - p_sym * tau_av) * k_latT * LT + phi_A * A - theta_A * AT - eta_T * AT
    dR <- f_rec * alpha_r * Ip + fT_rec * alpha_T * IpT + eta_r * A + eta_T * AT
    dDead <- (1 - f_rec) * alpha_r * Ip + (1 - fT_rec) * alpha_T * IpT

    list(c(dS, dST, dL, dLT, dIp, dIpT, dA, dAT, dR, dDead))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = flu_arino2008_antiviral_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", ST = "#333333", L = "#E8C468", LT = "#E8C468", Ip = "#B71C1C", IpT = "#B71C1C", A = "#E57373", AT = "#E57373", R = "#4E8A3E", Dead = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Pandemic influenza SLIAR model with antiviral treatment and pre-epidemic vaccination/prophylaxis (Arino et al. 2008)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
