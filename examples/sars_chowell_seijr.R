## SARS SEIJR model with diagnosis/isolation and two susceptibility classes (Chowell et al. 2003/2004)
## SARS (the template also suits MERS).
## Parameter values: 0 of 11 checked against the source; the rest are recalled or illustrative.
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S1, S2, E, Ip, J, R, Dead
## Flows (time unit: days):
##   S1 -> E     : foi * S1
##   S2 -> E     : p_low * foi * S2
##   E -> Ip     : k_inc * E
##   Ip -> J     : alpha_d * Ip
##   Ip -> R     : g1 * Ip
##   Ip -> Dead  : delta_d * Ip
##   J -> R      : g2 * J
##   J -> Dead   : delta_d * J

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 150, by = 0.2)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta_s = 0.25,     # transmission rate (1/day)
  q_e = 0.1,         # relative infectiousness of E
  l_iso = 0.43,      # relative infectiousness of J
  p_low = 1 / 3,     # relative susceptibility of S2
  rho_hr = 0.77,     # initial fraction in S1
  k_inc = 0.15707,   # E to Ip (6.4-day incubation) (1/day)
  alpha_d = 0.2061,  # diagnosis and isolation (4.85 days); Toronto used 1/3 after controls (1/day)
  g1 = 0.035285,     # recovery from Ip (1/day)
  g2 = 0.0426,       # recovery from J (1/day)
  delta_d = 0.0279,  # disease death rate (1/day)
  N0 = 6.7e6         # value taken from the notes text [recalled]
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  S1 <- rho_hr * N0        # high-risk susceptible
  S2 <- (1 - rho_hr) * N0  # lower-risk susceptible (relative susceptibility p_low)
  E <- 0                   # exposed, low infectivity
  Ip <- 1                  # symptomatic infectious, undiagnosed
  J <- 0                   # diagnosed and isolated
  R <- 0                   # recovered
  Dead <- 0                # SARS deaths
  c(S1 = S1, S2 = S2, E = E, Ip = Ip, J = J, R = R, Dead = Dead)
})

# ---- Model equations -------------------------------------------------
sars_chowell_seijr_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S1 + S2 + E + Ip + J + R

    # Derived quantities
    foi <- beta_s * (Ip + q_e * E + l_iso * J) / N

    # Rates of change
    dS1 <- -foi * S1
    dS2 <- -p_low * foi * S2
    dE <- foi * S1 + p_low * foi * S2 - k_inc * E
    dIp <- k_inc * E - alpha_d * Ip - g1 * Ip - delta_d * Ip
    dJ <- alpha_d * Ip - g2 * J - delta_d * J
    dR <- g1 * Ip + g2 * J
    dDead <- delta_d * Ip + delta_d * J

    list(c(dS1, dS2, dE, dIp, dJ, dR, dDead))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sars_chowell_seijr_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S1 = "#333333", S2 = "#333333", E = "#E8C468", Ip = "#B71C1C", J = "#333333", R = "#4E8A3E", Dead = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SARS SEIJR model with diagnosis/isolation and two susceptibility classes (Chowell et al. 2003/2004)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
