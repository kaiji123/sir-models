## Mpox with a rodent reservoir and residual smallpox vaccination immunity
## Mpox (monkeypox), rodent reservoir → humans.
## Two usable regimes from one structure: (a) the classical African setting — spillover-driven, R0h < 1, outbreaks that fade without the reservoir term; (b) a 2022-style outbreak — set beta_rh = 0, raise beta_hh above threshold, and read S as
## Parameter values: 0 of 15 checked against the source; the rest are recalled or illustrative.
## Reference: Bhunu CP, Mushayabasa S (2011). Modelling the transmission dynamics of pox-like infections. IAENG International Journal of Applied Mathematics 41(2):141–149. URL https://www.iaeng.org/IJAM/issues_v41/issue_2/IJAM_41_2_09.pdf . Historical human-to-human…
## Source: https://doi.org/10.1093/oxfordjournals.aje.a114747
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sr, Er, Ir, S, V, E, Ii, R, D
## Flows (time unit: days):
##   (outside) -> Sr  : br * Nr
##   Sr -> Er         : lambda_r * Sr
##   Sr -> (outside)  : mu_r * Sr
##   (outside) -> S   : bh
##   S -> E           : lambda_h * S
##   S -> V           : nu * S
##   V -> E           : lambda_h * (1 - eps_v) * V
##   V -> S           : wane_v * V
##   E -> Ii          : sigma_h * E
##   Ii -> R          : gamma_h * (1 - cfr) * Ii
##   Ii -> D          : gamma_h * cfr * Ii
##   S -> (outside)   : mu_h * S

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 1825, by = 2)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  br = 0.004,         # rodent birth rate [illustrative] (per day)
  mu_r = 0.004,       # rodent natural death rate [illustrative] (per day)
  beta_rr = 0.15,     # rodent-to-rodent transmission rate [illustrative] (per day)
  sigma_r = 1 / 7,    # rodent latency rate [illustrative] (per day)
  gamma_r = 1 / 21,   # rodent recovery/removal rate [illustrative] (per day)
  beta_rh = 2e-6,     # rodent-to-human spillover coefficient (hunting/bushmeat) [illustrative] (per rodent per day)
  beta_hh = 0.035,    # human-to-human transmission rate [illustrative] (per day)
  sigma_h = 1 / 8.5,  # human latency rate (incubation ≈8.5 d) [illustrative] (per day)
  gamma_h = 1 / 21,   # human removal rate (rash duration 2–4 weeks) [illustrative] (per day)
  cfr = 0.035,        # human case fatality (clade I historical ≈1–10%) [illustrative] (dimensionless)
  nu = 0,             # vaccination rate of susceptibles (ring/targeted campaign lever) [illustrative] (per day)
  eps_v = 0.85,       # vaccine protection against infection (smallpox vaccine cross-protection ≈85%) [unverified] (dimensionless)
  wane_v = 1 / 7300,  # waning of vaccine protection (20-year mean) [illustrative] (per day)
  mu_h = 1 / 21900,   # human natural death rate [illustrative] (per day)
  bh = 5              # human recruitment [illustrative] (humans per day)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sr = 10000,  # susceptible rodents
  Er = 0,      # latent rodents
  Ir = 100,    # infectious rodents
  S = 1e5,     # susceptible humans
  V = 2e4,     # humans with vaccine-derived (smallpox) protection
  E = 0,       # latent humans
  Ii = 0,      # infectious humans
  R = 0,       # recovered humans
  D = 0        # cumulative human mpox deaths
)

# ---- Model equations -------------------------------------------------
mpox_human_animal_reservoir_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N_total <- S + V + E + Ii + R

    # Derived quantities
    Nr <- Sr + Er + Ir                           # rodent population
    N <- S + V + E + Ii + R                      # living human population
    lambda_r <- beta_rr * Ir / Nr                # rodent force of infection
    lambda_h <- beta_rh * Ir + beta_hh * Ii / N  # human force of infection
    R0h <- beta_hh / (gamma_h + mu_h)            # human-to-human reproduction number (historically ≈0.3; a 2022-type MSM-outbreak scenario is produced by raising beta_hh so that R0h > 1)

    # Rates of change
    dSr <- br * Nr - lambda_r * Sr - mu_r * Sr
    dEr <- lambda_r * Sr
    dIr <- 0
    dS <- bh - lambda_h * S - nu * S + wane_v * V - mu_h * S
    dV <- nu * S - lambda_h * (1 - eps_v) * V - wane_v * V
    dE <- lambda_h * S + lambda_h * (1 - eps_v) * V - sigma_h * E
    dIi <- sigma_h * E - gamma_h * (1 - cfr) * Ii - gamma_h * cfr * Ii
    dR <- gamma_h * (1 - cfr) * Ii
    dD <- gamma_h * cfr * Ii

    list(c(dSr, dEr, dIr, dS, dV, dE, dIi, dR, dD))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = mpox_human_animal_reservoir_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sr = "#333333", Er = "#E8C468", Ir = "#B71C1C", S = "#333333", V = "#2457C5", E = "#E8C468", Ii = "#B71C1C", R = "#4E8A3E", D = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Mpox with a rodent reservoir and residual smallpox vaccination immunity")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
