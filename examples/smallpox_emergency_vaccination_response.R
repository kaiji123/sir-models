## Smallpox outbreak response: traced (ring) versus mass vaccination
## Smallpox, human; emergency public-health vaccination response.
## This is a public-health response-logistics model: the policy question it answers is whether contact-traced vaccination or population-wide vaccination averts more deaths given a finite daily vaccination capacity.
## Parameter values: 0 of 11 checked against the source; the rest are recalled or illustrative.
## Reference: Kaplan EH, Craft DL, Wein LM (2002). Emergency response to a smallpox attack: the case for mass vaccination. PNAS 99(16):10935–10940
## Source: https://doi.org/10.1073/pnas.162282799
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, Ii, Q, V, R, D, CV
## Flows (time unit: days):
##   S -> E           : lambda * S
##   S -> V           : rate_v * eps_pre * S
##   E -> V           : rate_v * eps_post * E
##   E -> Ii          : sigma * E
##   Ii -> Q          : theta * Ii
##   Q -> D           : gamma * cfr * Q
##   Q -> R           : gamma * (1 - cfr) * Q
##   (outside) -> CV  : pmin(want, cap)

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 180, by = 0.2)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.105,    # transmission rate (gives roughly 3 secondary cases per infectious person over the symptomatic period) [illustrative] (per day)
  sigma = 1 / 14,  # rate of leaving latency (mean latent period ≈14 d) [unverified] (per day)
  theta = 1 / 3,   # rate at which symptomatic cases are detected and isolated (≈3 d from rash to isolation) [illustrative] (per day)
  gamma = 1 / 14,  # rate of resolution (recovery or death) in isolation [illustrative] (per day)
  cfr = 0.3,       # case fatality of ordinary variola major [unverified] (dimensionless)
  eps_pre = 0.95,  # efficacy of vaccination of a susceptible [illustrative] (dimensionless)
  eps_post = 0.5,  # efficacy of post-exposure vaccination of a latent individual (protection declines sharply after the first few days) [illustrative] (dimensionless)
  tau_v = 7,       # delay before the vaccination campaign starts [illustrative] (days)
  nu_ring = 4,     # traced-vaccination intensity [illustrative] (contacts traced and vaccinated per isolated case per day)
  nu_mass = 0,     # per-capita mass-vaccination rate once the campaign is running [illustrative] (per day)
  cap = 2e5        # maximum daily vaccination throughput [illustrative] (doses per day)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 1e7,  # susceptible, not yet vaccinated
  E = 0,    # infected, latent (not infectious, not yet symptomatic)
  Ii = 10,  # infectious, in the community
  Q = 0,    # isolated/removed from circulation (hospitalised or quarantined)
  V = 0,    # vaccinated and protected
  R = 0,    # recovered, immune
  D = 0,    # cumulative deaths
  CV = 0    # cumulative vaccine doses delivered
)

# ---- Model equations -------------------------------------------------
smallpox_emergency_vaccination_response_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N_total <- S + E + Ii + Q + V + R

    # Derived quantities
    N <- S + E + Ii + Q + V + R                                  # living population
    active <- ifelse(t < tau_v, 0, 1)                            # campaign switch
    want <- active * (nu_mass * (S + E) + nu_ring * theta * Ii)  # doses demanded per day
    rate_v <- pmin(want, cap) / pmax(S + E, 1)                   # per-capita vaccination rate, capped by logistic throughput
    lambda <- beta * Ii / N                                      # force of infection

    # Rates of change
    dS <- -lambda * S - rate_v * eps_pre * S
    dE <- lambda * S - rate_v * eps_post * E - sigma * E
    dIi <- sigma * E - theta * Ii
    dQ <- theta * Ii - gamma * cfr * Q - gamma * (1 - cfr) * Q
    dV <- rate_v * eps_pre * S + rate_v * eps_post * E
    dR <- gamma * (1 - cfr) * Q
    dD <- gamma * cfr * Q
    dCV <- pmin(want, cap)

    list(c(dS, dE, dIi, dQ, dV, dR, dD, dCV))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = smallpox_emergency_vaccination_response_model, parms = parameters, hmax = 0.2)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", Ii = "#B71C1C", Q = "#00838F", V = "#2457C5", R = "#4E8A3E", D = "#111111", CV = "#8E0E0E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Smallpox outbreak response: traced (ring) versus mass vaccination")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
