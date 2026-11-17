## SIDARTHE model of COVID-19 in Italy (Giordano et al. 2020)
## COVID-19.
## Parameter values: 0 of 12 checked against the source; the rest are recalled or illustrative.
## Reference: Giordano G, Blanchini F, Bruno R, Colaneri P, Di Filippo A, Di Matteo A, Colaneri M. "Modelling the COVID-19 epidemic and implementation of population-wide interventions in Italy." Nat Med 2020;26:855–860. https://doi.org/10.1038/s41591-020-0883-7
## Source: https://doi.org/10.1038/s41591-020-0883-7
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, Iu, Dg, A, Rd, Th, H, Ex
## Flows (time unit: days):
##   S -> Iu   : foi * S
##   Iu -> Dg  : eps_v * Iu
##   Iu -> A   : zeta * Iu
##   Iu -> H   : lambda_h * Iu
##   Dg -> Rd  : eta * Dg
##   Dg -> H   : rho_h * Dg
##   A -> Rd   : theta * A
##   A -> Th   : mu_t * A
##   A -> H    : kappa_h * A
##   Rd -> Th  : nu_t * Rd
##   Rd -> H   : xi_h * Rd
##   Th -> Ex  : tau_t * Th
##   Th -> H   : sigma_h * Th

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 300, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  Npop = 60e6,       # value taken from the notes text [recalled]
  zeta = 0.125,      # value taken from the notes text [recalled]
  lambda_h = 0.034,  # value taken from the notes text [recalled]
  eta = 0.125,       # value taken from the notes text [recalled]
  rho_h = 0.034,     # value taken from the notes text [recalled]
  theta = 0.371,     # value taken from the notes text [recalled]
  mu_t = 0.012,      # value taken from the notes text [recalled]
  kappa_h = 0.017,   # value taken from the notes text [recalled]
  nu_t = 0.027,      # value taken from the notes text [recalled]
  xi_h = 0.017,      # value taken from the notes text [recalled]
  tau_t = 0.003,     # value taken from the notes text [recalled]
  sigma_h = 0.017    # value taken from the notes text [recalled]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 60e6 - 223,  # susceptible
  Iu = 200,        # infected, asymptomatic or mild, undetected (I)
  Dg = 20,         # diagnosed, asymptomatic or mild (D)
  A = 1,           # ailing: symptomatic, undetected (A)
  Rd = 2,          # recognised: symptomatic, detected (R)
  Th = 0,          # threatened: life-threatening, detected, in ICU (T)
  H = 0,           # healed
  Ex = 0           # extinct (dead)
)

# ---- Model equations -------------------------------------------------
covid_sidarthe_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + Iu + Dg + A + Rd + Th + H

    # Derived quantities
    alpha_v <- ifelse(t >= 3, 0.422, 0.570)
    gamma_v <- ifelse(t >= 3, 0.285, 0.456)
    bd_v <- ifelse(t >= 3, 0.0057, 0.011)
    eps_v <- ifelse(t >= 11, 0.143, 0.171)
    foi <- (alpha_v * Iu + bd_v * Dg + gamma_v * A + bd_v * Rd) / Npop

    # Rates of change
    dS <- -foi * S
    dIu <- foi * S - eps_v * Iu - zeta * Iu - lambda_h * Iu
    dDg <- eps_v * Iu - eta * Dg - rho_h * Dg
    dA <- zeta * Iu - theta * A - mu_t * A - kappa_h * A
    dRd <- eta * Dg + theta * A - nu_t * Rd - xi_h * Rd
    dTh <- mu_t * A + nu_t * Rd - tau_t * Th - sigma_h * Th
    dH <- lambda_h * Iu + rho_h * Dg + kappa_h * A + xi_h * Rd + sigma_h * Th
    dEx <- tau_t * Th

    list(c(dS, dIu, dDg, dA, dRd, dTh, dH, dEx))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = covid_sidarthe_model, parms = parameters, hmax = 0.5)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", Iu = "#B71C1C", Dg = "#111111", A = "#E57373", Rd = "#4E8A3E", Th = "#333333", H = "#6A2C9C", Ex = "#E8C468")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SIDARTHE model of COVID-19 in Italy (Giordano et al. 2020)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
