## COVID-19 SEIR-H with two-level vaccine protection, waning and reinfection
## COVID-19.
## vh_red would need a separate exposed class for vaccinated people (E_v) to act properly.
## Parameter values: 0 of 13 checked against the source; the rest are recalled or illustrative.
## Reference: [ILLUSTRATIVE] generic structure of the kind used in many 2021 vaccination models. Values are chosen to be plausible, not taken from a single paper. Waning anchor: Andrews N et al. "Covid-19 vaccine effectiveness against the Omicron (B.1.1.529) variant." N…
## Source: https://doi.org/10.1056/NEJMoa2119451
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, V, Wv, E, Ip, Hc, R, Dead
## Flows (time unit: days):
##   S -> E      : foi * S
##   V -> E      : (1 - ve_V) * foi * V
##   Wv -> E     : (1 - ve_W) * foi * Wv
##   S -> V      : vac_rate * S
##   V -> Wv     : wane_V * V
##   Wv -> V     : boost_rate * Wv
##   R -> Wv     : wane_R * R
##   E -> Ip     : sig_lat * E
##   Ip -> Hc    : p_hosp * g_rec * Ip
##   Ip -> R     : (1 - p_hosp) * g_rec * Ip
##   Hc -> R     : (1 - p_die) * h_rec * Hc
##   Hc -> Dead  : p_die * h_rec * Hc

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 730, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  R0_eff = 3,         # illustrative
  sig_lat = 1 / 3,    # latent period of 3 days (1/day) [illustrative]
  g_rec = 1 / 5,      # infectious period of 5 days (1/day) [illustrative]
  p_hosp = 0.02,      # hospitalised fraction, unvaccinated [illustrative]
  h_rec = 1 / 8,      # hospital stay (1/day) [illustrative]
  p_die = 0.15,       # in-hospital fatality [illustrative]
  vac_rate = 0.005,   # per-capita vaccination rate of S (1/day) [illustrative]
  ve_V = 0.8,         # leaky protection against infection [illustrative]
  ve_W = 0.3,         # leaky protection against infection [illustrative]
  vh_red = 0.8,       # extra reduction in hospitalisation if infected after vaccination [illustrative]
  wane_V = 1 / 120,   # waning from V to Wv (1/day) [illustrative]
  wane_R = 1 / 180,   # waning of infection immunity, R to Wv (1/day) [illustrative]
  boost_rate = 0.002  # booster from Wv back to V (1/day) [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 0.7 * 1e6,  # susceptible, unvaccinated
  V = 0,          # recently vaccinated (high protection)
  Wv = 0,         # vaccinated, waned protection
  E = 100,        # exposed
  Ip = 100,       # infectious
  Hc = 0,         # hospitalised
  R = 0.3 * 1e6,  # recovered (hybrid or natural immunity)
  Dead = 0        # COVID deaths
)

# ---- Model equations -------------------------------------------------
covid_vax_waning_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + V + Wv + E + Ip + Hc + R

    # Derived quantities
    foi <- R0_eff * g_rec * Ip / N

    # Rates of change
    dS <- -foi * S - vac_rate * S
    dV <- -(1 - ve_V) * foi * V + vac_rate * S - wane_V * V + boost_rate * Wv
    dWv <- -(1 - ve_W) * foi * Wv + wane_V * V - boost_rate * Wv + wane_R * R
    dE <- foi * S + (1 - ve_V) * foi * V + (1 - ve_W) * foi * Wv - sig_lat * E
    dIp <- sig_lat * E - p_hosp * g_rec * Ip - (1 - p_hosp) * g_rec * Ip
    dHc <- p_hosp * g_rec * Ip - (1 - p_die) * h_rec * Hc - p_die * h_rec * Hc
    dR <- -wane_R * R + (1 - p_hosp) * g_rec * Ip + (1 - p_die) * h_rec * Hc
    dDead <- p_die * h_rec * Hc

    list(c(dS, dV, dWv, dE, dIp, dHc, dR, dDead))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = covid_vax_waning_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", V = "#2457C5", Wv = "#A5D6A7", E = "#E8C468", Ip = "#B71C1C", Hc = "#6A2C9C", R = "#4E8A3E", Dead = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "COVID-19 SEIR-H with two-level vaccine protection, waning and reinfection")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
