## Syphilis with incubating, primary, secondary, early/late latent stages and waning immunity
## Syphilis.
## Parameter values: 6 of 11 checked against the source; the rest are recalled or illustrative.
## Reference: Garnett GP, Aral SO, Hoyle DV, Cates W Jr, Anderson RM. The natural history of syphilis: implications for the transmission dynamics and control of infection. Sex Transm Dis 1997;24(4):185–200
## Source: https://doi.org/10.1097/00007435-199704000-00002
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, P1, P2, Le, Ll, R
## Flows (time unit: years):
##   (outside) -> S   : mu * N
##   S -> E           : foi * S
##   E -> P1          : sig_inc * E
##   P1 -> P2         : 1 / (f_prim * d_inf) * P1
##   P2 -> Le         : 1 / ((1 - f_prim) * d_inf) * P2
##   Le -> P2         : rel * Le
##   Le -> Ll         : el * Le
##   P1 -> R          : tau * P1
##   P2 -> R          : tau * P2
##   Le -> R          : tau * Le
##   Ll -> R          : tau_l * Ll
##   R -> S           : w_imm * R
##   S -> (outside)   : mu * S
##   E -> (outside)   : mu * E
##   P1 -> (outside)  : mu * P1
##   P2 -> (outside)  : mu * P2
##   Le -> (outside)  : mu * Le
##   Ll -> (outside)  : mu * Ll
##   R -> (outside)   : mu * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 30, by = 0.05)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  cp = 15,             # mean partner change rate (MSM, Ontario) (partners/yr) [verified]
  beta_s = 0.3,        # transmission probability (per partnership) [illustrative]
  sig_inc = 365 / 25,  # 1/incubation period (25 days) (/yr) [verified]
  d_inf = 154 / 365,   # total primary + secondary infectious period (yr) [verified]
  f_prim = 0.3,        # fraction of d_inf spent in primary stage [illustrative]
  tau = 0.4,           # treatment rate of infectious and early-latent cases (/yr) [verified]
  tau_l = 0.1,         # treatment rate in late latent (/yr) [illustrative]
  rel = 0.25,          # relapse rate early latent to secondary (/yr) [illustrative]
  el = 1,              # progression early to late latent (early latent about 1 yr) (/yr) [illustrative]
  w_imm = 1 / 5,       # waning of post-treatment immunity (5 yr) (/yr) [verified]
  mu = 5831 / 366279   # entry/exit rate (/yr) [verified]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 0.99,    # susceptible
  E = 0,       # incubating
  P1 = 0.005,  # primary (infectious)
  P2 = 0.005,  # secondary (infectious)
  Le = 0,      # early latent (non-infectious; relapse possible)
  Ll = 0,      # late latent
  R = 0        # treated, temporarily immune
)

# ---- Model equations -------------------------------------------------
syphilis_staged_garnett_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + P1 + P2 + Le + Ll + R

    # Derived quantities
    foi <- cp * beta_s * (P1 + P2) / N  # force of infection

    # Rates of change
    dS <- mu * N - foi * S + w_imm * R - mu * S
    dE <- foi * S - sig_inc * E - mu * E
    dP1 <- sig_inc * E - 1 / (f_prim * d_inf) * P1 - tau * P1 - mu * P1
    dP2 <- 1 / (f_prim * d_inf) * P1 - 1 / ((1 - f_prim) * d_inf) * P2 + rel * Le - tau * P2 - mu * P2
    dLe <- 1 / ((1 - f_prim) * d_inf) * P2 - rel * Le - el * Le - tau * Le - mu * Le
    dLl <- el * Le - tau_l * Ll - mu * Ll
    dR <- tau * P1 + tau * P2 + tau * Le + tau_l * Ll - w_imm * R - mu * R

    list(c(dS, dE, dP1, dP2, dLe, dLl, dR))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = syphilis_staged_garnett_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", P1 = "#F0A35E", P2 = "#F0A35E", Le = "#E8C468", Ll = "#E8C468", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Syphilis with incubating, primary, secondary, early/late latent stages and waning immunity")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
