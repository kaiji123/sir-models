## African swine fever in wild boar with carcass-mediated environmental transmission
## African swine fever virus, Eurasian wild boar (Sus scrofa).
## The key structural feature is that dead infectious hosts move into a long-lived environmental compartment rather than out of the system: Ii > Cr is the flow that lets ASF persist at very low host density, because Cr decays on a timescale of
## Parameter values: 0 of 11 checked against the source; the rest are recalled or illustrative.
## Reference: O'Neill X, White A, Ruiz-Fons F, Gortázar C (2020). Modelling the transmission and persistence of African swine fever in wild boar in contrasting European scenarios. Scientific Reports 10:5895
## Source: https://doi.org/10.1038/s41598-020-62736-y
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, Ii, R, Cr
## Flows (time unit: days):
##   (outside) -> S   : pmax(dens, 0) * N
##   S -> E           : lambda * S
##   E -> Ii          : sigma * E
##   Ii -> Cr         : gamma * cfr * Ii
##   Ii -> R          : gamma * (1 - cfr) * Ii
##   S -> (outside)   : (mub + hunt) * S
##   Ii -> (outside)  : (mub + hunt) * Ii
##   Cr -> (outside)  : (decay_c + remove_c) * Cr

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 1095, by = 2)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta_dir = 0.3,    # direct boar-to-boar transmission rate [illustrative] (per day)
  beta_car = 0.05,   # carcass-to-boar transmission coefficient (rooting/sniffing contact) [illustrative] (per carcass per day)
  sigma = 1 / 4,     # rate of leaving the exposed class (incubation ≈4 d) [illustrative] (per day)
  gamma = 1 / 5,     # rate of leaving infectiousness (≈5 d) [illustrative] (per day)
  cfr = 0.95,        # case fatality (ASF genotype II is near-uniformly lethal) [unverified] (dimensionless)
  decay_c = 1 / 90,  # carcass decay/removal rate (seasonally variable; much slower in winter) [illustrative] (per day)
  remove_c = 0,      # active carcass-removal rate (the main control lever) [illustrative] (per day)
  bb = 0.004,        # boar birth rate (strongly seasonal in reality) [illustrative] (per day)
  mub = 0.002,       # boar natural death rate [illustrative] (per day)
  hunt = 0,          # hunting/culling rate [illustrative] (per day)
  Kb = 2000          # carrying capacity [illustrative] (boar)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 2000,  # susceptible wild boar
  E = 0,     # exposed (infected, not yet infectious)
  Ii = 2,    # infectious wild boar
  R = 0,     # recovered/immune survivors
  Cr = 0     # infectious carcasses in the environment
)

# ---- Model equations -------------------------------------------------
asf_wild_boar_carcass_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N_total <- S + E + Ii + R

    # Derived quantities
    N <- S + E + Ii + R                          # living boar population
    lambda <- beta_dir * Ii / N + beta_car * Cr  # total force of infection (frequency-dependent direct term plus density-dependent carcass term)
    dens <- bb * (1 - N / Kb)                    # density-dependent per-capita birth rate

    # Rates of change
    dS <- pmax(dens, 0) * N - lambda * S - (mub + hunt) * S
    dE <- lambda * S - sigma * E
    dIi <- sigma * E - gamma * cfr * Ii - gamma * (1 - cfr) * Ii - (mub + hunt) * Ii
    dR <- gamma * (1 - cfr) * Ii
    dCr <- gamma * cfr * Ii - (decay_c + remove_c) * Cr

    list(c(dS, dE, dIi, dR, dCr))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = asf_wild_boar_carcass_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", Ii = "#B71C1C", R = "#4E8A3E", Cr = "#8E0E0E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "African swine fever in wild boar with carcass-mediated environmental transmission")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
