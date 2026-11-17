## Ebola West Africa 2014 with case isolation and vaccination (Legrand extension)
## Ebola virus disease, human, Liberia/Sierra Leone 2014.
## Overlaps substantially with the already-implemented Legrand 2007 model — the compartment skeleton (S,E,I,H,F,R) is the same; the additions are the V compartment with a vaccination lever and the 2014 West-African parameterisation.
## Parameter values: 0 of 13 checked against the source; the rest are recalled or illustrative.
## Reference: Rivers CM, Lofgren ET, Marathe M, Eubank S, Lewis BL (2014). Modeling the impact of interventions on an epidemic of Ebola in Sierra Leone and Liberia. PLOS Currents Outbreaks, Oct 16 ed
## Source: https://doi.org/10.1371/currents.outbreaks.fd38dd85078565450b0be3fcd78f5ccf
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, I, H, Fn, R, V
## Flows (time unit: days):
##   S -> E   : lambda * S
##   S -> V   : nu * eff * S
##   E -> I   : alpha * E
##   I -> H   : gamma_h * theta * I
##   I -> Fn  : gamma_d * cfr * (1 - theta) * I
##   I -> R   : gamma_i * (1 - cfr) * (1 - theta) * I
##   H -> Fn  : gamma_dh * cfr * H
##   H -> R   : gamma_dh * (1 - cfr) * H
##   Fn -> R  : gamma_f * Fn

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 365, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta_i = 0.16,       # community transmission rate [unverified] (per day)
  beta_h = 0.062,      # transmission rate from hospitalised cases (reflects imperfect infection control) [unverified] (per day)
  beta_f = 0.489,      # transmission rate at funerals/from unburied bodies [unverified] (per day)
  alpha = 1 / 11,      # rate of leaving latency (mean incubation ≈11 d) [unverified] (per day)
  theta = 0.248,       # proportion of cases hospitalised/isolated [unverified] (dimensionless)
  gamma_h = 1 / 4.12,  # rate community case → hospital [unverified] (per day)
  gamma_i = 1 / 10,    # rate of leaving community infectiousness [unverified] (per day)
  gamma_d = 1 / 8,     # rate infectious → death (community) [unverified] (per day)
  gamma_dh = 1 / 6,    # rate hospitalised → death [unverified] (per day)
  gamma_f = 1 / 2.01,  # rate of burial (removal of Fn) [unverified] (per day)
  cfr = 0.7,           # case fatality proportion [unverified] (dimensionless)
  nu = 0,              # per-capita vaccination rate (intervention lever; set >0 to simulate a campaign) [illustrative] (per day)
  eff = 0.9            # vaccine efficacy (take) [illustrative] (dimensionless)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 1e6,  # susceptible
  E = 0,    # latent (infected, not infectious)
  I = 10,   # infectious in the community
  H = 0,    # hospitalised / isolated infectious
  Fn = 0,   # deceased, not yet buried (funeral-transmission compartment)
  R = 0,    # removed (recovered or safely buried)
  V = 0     # vaccinated, immune
)

# ---- Model equations -------------------------------------------------
ebola_rivers_2014_intervention_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N_total <- S + E + I + H + Fn + R + V

    # Derived quantities
    N <- S + E + I + H + R + V                             # living population
    lambda <- (beta_i * I + beta_h * H + beta_f * Fn) / N  # total force of infection

    # Rates of change
    dS <- -lambda * S - nu * eff * S
    dE <- lambda * S - alpha * E
    dI <- alpha * E - gamma_h * theta * I - gamma_d * cfr * (1 - theta) * I - gamma_i * (1 - cfr) * (1 - theta) * I
    dH <- gamma_h * theta * I - gamma_dh * cfr * H - gamma_dh * (1 - cfr) * H
    dFn <- gamma_d * cfr * (1 - theta) * I + gamma_dh * cfr * H - gamma_f * Fn
    dR <- gamma_i * (1 - cfr) * (1 - theta) * I + gamma_dh * (1 - cfr) * H + gamma_f * Fn
    dV <- nu * eff * S

    list(c(dS, dE, dI, dH, dFn, dR, dV))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = ebola_rivers_2014_intervention_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", I = "#B71C1C", H = "#6A2C9C", Fn = "#5D4037", R = "#4E8A3E", V = "#2457C5")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Ebola West Africa 2014 with case isolation and vaccination (Legrand extension)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
