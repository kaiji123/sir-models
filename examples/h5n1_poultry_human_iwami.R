## Avian influenza in poultry with human spillover and a mutant human-transmissible strain
## Highly pathogenic avian influenza H5N1, poultry → humans.
## The paper's point is that Ii is a dead-end (no human-to-human spread) while Im is not, so two outbreak regimes are possible depending on whether R0m exceeds 1.
## Parameter values: 0 of 13 checked against the source; the rest are recalled or illustrative.
## Reference: Iwami S, Takeuchi Y, Liu X (2007). Avian–human influenza epidemic model. Mathematical Biosciences 207(1):1–25
## Source: https://doi.org/10.1016/j.mbs.2006.08.001
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sa, Ia, S, Ii, Im, R, D
## Flows (time unit: days):
##   (outside) -> Sa  : ca
##   Sa -> Ia         : lambda_a * Sa
##   Sa -> (outside)  : (mu_a + cull) * Sa
##   Ia -> (outside)  : (mu_a + delta_a + cull) * Ia
##   (outside) -> S   : ch
##   S -> Ii          : lambda_ah * S
##   S -> Im          : lambda_m * S
##   Ii -> Im         : nu * Ii
##   Ii -> (outside)  : (mu_h + delta_i) * Ii
##   Im -> R          : gamma_m * Im
##   Im -> D          : delta_m * Im
##   Im -> (outside)  : mu_h * Im
##   S -> (outside)   : mu_h * S

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 365, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  ca = 1000,         # recruitment of susceptible birds (poultry restocking); named ca because c is reserved [illustrative] (birds per day)
  mu_a = 1 / 500,    # natural bird death rate [illustrative] (per day)
  beta_a = 2e-5,     # bird-to-bird transmission coefficient (mass action) [illustrative] (per bird per day)
  delta_a = 0.4,     # disease-induced bird death rate (HPAI kills within days) [illustrative] (per day)
  cull = 0,          # additional culling rate applied to birds (policy lever from the 2008 companion paper) [illustrative] (per day)
  ch = 100,          # human recruitment [illustrative] (humans per day)
  mu_h = 1 / 25550,  # natural human death rate (70-year life expectancy) [illustrative] (per day)
  beta_ah = 1e-8,    # bird-to-human transmission coefficient [illustrative] (per bird per day)
  beta_m = 2e-7,     # mutant human-to-human transmission coefficient [illustrative] (per human per day)
  nu = 1e-4,         # rate at which an avian-strain human case yields mutant-strain infection (mutation) [illustrative] (per day)
  delta_i = 0.5,     # disease death rate, avian-strain human cases (H5N1 CFR ≈60%) [illustrative] (per day)
  delta_m = 0.1,     # disease death rate, mutant-strain cases [illustrative] (per day)
  gamma_m = 1 / 7    # recovery rate, mutant-strain cases [illustrative] (per day)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sa = 1e5,  # susceptible birds
  Ia = 10,   # infected birds
  S = 1e6,   # susceptible humans
  Ii = 0,    # humans infected with the avian (non-transmissible) strain
  Im = 0,    # humans infected with the mutant human-transmissible strain
  R = 0,     # recovered humans
  D = 0      # cumulative human deaths
)

# ---- Model equations -------------------------------------------------
h5n1_poultry_human_iwami_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N_total <- S + Ii + Im + R

    # Derived quantities
    Na <- Sa + Ia                                             # bird population
    N <- S + Ii + Im + R                                      # living human population
    lambda_a <- beta_a * Ia                                   # bird force of infection (density-dependent)
    lambda_ah <- beta_ah * Ia                                 # spillover hazard to humans
    lambda_m <- beta_m * Im                                   # mutant-strain force of infection
    R0a <- beta_a * ca / ((mu_a + delta_a + cull) * mu_a)     # bird-strain reproduction number
    R0m <- beta_m * ch / ((mu_h + delta_m + gamma_m) * mu_h)  # mutant-strain reproduction number

    # Rates of change
    dSa <- ca - lambda_a * Sa - (mu_a + cull) * Sa
    dIa <- lambda_a * Sa - (mu_a + delta_a + cull) * Ia
    dS <- ch - lambda_ah * S - lambda_m * S - mu_h * S
    dIi <- lambda_ah * S - nu * Ii - (mu_h + delta_i) * Ii
    dIm <- lambda_m * S + nu * Ii - gamma_m * Im - delta_m * Im - mu_h * Im
    dR <- gamma_m * Im
    dD <- delta_m * Im

    list(c(dSa, dIa, dS, dIi, dIm, dR, dD))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = h5n1_poultry_human_iwami_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sa = "#333333", Ia = "#B71C1C", S = "#333333", Ii = "#B71C1C", Im = "#B71C1C", R = "#4E8A3E", D = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Avian influenza in poultry with human spillover and a mutant human-transmissible strain")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
