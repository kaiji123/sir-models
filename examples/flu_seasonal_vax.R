## Seasonal influenza SIRS with sinusoidal forcing and an annual vaccination campaign
## Seasonal influenza.
## vaccinating only S assumes vaccine is not wasted on R.
## Parameter values: 0 of 12 checked against the source; the rest are recalled or illustrative.
## Reference: [ILLUSTRATIVE] textbook SIRS-V structure (Keeling & Rohani 2008, Modeling Infectious Diseases in Humans and Animals, Princeton UP, ch. 5). The reproduction-number anchor is Biggerstaff M, Cauchemez S, Reed C, Gambhir M, Finelli L. "Estimates of the…
## Source: https://doi.org/10.1186/1471-2334-14-480
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, Ip, R, V
## Flows (time unit: days):
##   (outside) -> S  : mu * N
##   S -> Ip         : foi * S
##   Ip -> R         : g_rec * Ip
##   R -> S          : omega_R * R
##   S -> V          : ve * vrate * S
##   V -> S          : omega_V * V
##   S -> (outside)  : mu * S

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 3650, by = 5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  R0_base = 1.3,      # mean reproduction number
  g_rec = 0.25,       # recovery rate (4 days) (1/day) [illustrative]
  amp = 0.2,          # seasonal amplitude [illustrative]
  t_peak = 15,        # day of peak transmission (day of year) [illustrative]
  two_pi = 6.283185,  # constant
  omega_R = 1 / 730,  # waning of infection immunity, including drift (1/day) [illustrative]
  omega_V = 1 / 240,  # waning of vaccine protection (1/day) [illustrative]
  cov = 0.45,         # target coverage per campaign [illustrative]
  ve = 0.5,           # vaccine efficacy (all-or-nothing) [illustrative]
  camp_start = 274,   # campaign from Oct 1 for 60 days [illustrative]
  camp_len = 60,      # campaign from Oct 1 for 60 days (day of year, days) [illustrative]
  mu = 3.9e-5         # birth = death rate (1/day) [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 0.6 * 1e6,  # susceptible
  Ip = 100,       # infectious
  R = 0.4 * 1e6,  # infection-derived immunity
  V = 0           # vaccine-protected
)

# ---- Model equations -------------------------------------------------
flu_seasonal_vax_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + Ip + R + V

    # Derived quantities
    beta_t <- R0_base * g_rec * (1 + amp * cos(two_pi * (t - t_peak) / 365))
    tyear <- t - 365 * floor(t / 365)
    in_camp <- (tyear >= camp_start) * (tyear < camp_start + camp_len)
    vrate <- in_camp * (-log(1 - cov) / camp_len)
    foi <- beta_t * Ip / N

    # Rates of change
    dS <- mu * N - foi * S + omega_R * R - ve * vrate * S + omega_V * V - mu * S
    dIp <- foi * S - g_rec * Ip
    dR <- g_rec * Ip - omega_R * R
    dV <- ve * vrate * S - omega_V * V

    list(c(dS, dIp, dR, dV))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = flu_seasonal_vax_model, parms = parameters, hmax = 5)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", Ip = "#B71C1C", R = "#4E8A3E", V = "#2457C5")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Seasonal influenza SIRS with sinusoidal forcing and an annual vaccination campaign")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
