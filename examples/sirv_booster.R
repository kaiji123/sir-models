## vaccination with waning and booster doses
## pertussis, COVID-19, HPV-like.
## N is constant. Compare b_rate = 0 against b_rate > 0.
## Parameter values: 0 of 11 checked against the source; the rest are recalled or illustrative.
## Reference: Illustrative structure. The rationale for waning and boosting follows Hethcote 2000 (SIAM Rev
## Source: https://doi.org/10.1137/S0036144500371907
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, V1, W, V2, I, R
## Flows (time unit: years):
##   (outside) -> S   : (1 - p_nb) * mu * N
##   (outside) -> V1  : p_nb * mu * N
##   V1 -> W          : omega1 * V1
##   V2 -> W          : omega2 * V2
##   R -> W           : omega_r * R
##   W -> V2          : b_rate * W
##   S -> I           : lambda * S
##   V1 -> I          : eps1 * lambda * V1
##   V2 -> I          : eps2 * lambda * V2
##   W -> I           : eps_w * lambda * W
##   I -> R           : gamma * I
##   S -> (outside)   : mu * S
##   V1 -> (outside)  : mu * V1
##   W -> (outside)   : mu * W
##   V2 -> (outside)  : mu * V2
##   I -> (outside)   : mu * I
##   R -> (outside)   : mu * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 50, by = 0.05)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  mu = 1 / 70,       # birth/death (1/yr) [illustrative]
  p_nb = 0.9,        # infant primary coverage [illustrative]
  omega1 = 1 / 4,    # waning of primary protection (1/yr) [illustrative]
  omega2 = 1 / 10,   # waning after booster (1/yr) [illustrative]
  omega_r = 1 / 20,  # waning of infection-derived immunity, flowing into W (1/yr) [illustrative]
  b_rate = 0.2,      # booster uptake among waned (1/yr) [illustrative]
  eps1 = 0.1,        # relative susceptibility in V1 [illustrative]
  eps2 = 0.1,        # relative susceptibility in V2 [illustrative]
  eps_w = 0.5,       # relative susceptibility of the waned [illustrative]
  beta = 300,        # transmission (1/yr) [illustrative]
  gamma = 365 / 21   # recovery (1/yr) [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9000,  # never vaccinated, never infected
  V1 = 0,    # primary-vaccinated (high protection)
  W = 0,     # waned (partially susceptible)
  V2 = 0,    # boosted
  I = 10,    # infectious
  R = 990    # recovered
)

# ---- Model equations -------------------------------------------------
sirv_booster_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + V1 + W + V2 + I + R

    # Derived quantities
    lambda <- beta * I / N

    # Rates of change
    dS <- (1 - p_nb) * mu * N - lambda * S - mu * S
    dV1 <- p_nb * mu * N - omega1 * V1 - eps1 * lambda * V1 - mu * V1
    dW <- omega1 * V1 + omega2 * V2 + omega_r * R - b_rate * W - eps_w * lambda * W - mu * W
    dV2 <- -omega2 * V2 + b_rate * W - eps2 * lambda * V2 - mu * V2
    dI <- lambda * S + eps1 * lambda * V1 + eps2 * lambda * V2 + eps_w * lambda * W - gamma * I - mu * I
    dR <- -omega_r * R + gamma * I - mu * R

    list(c(dS, dV1, dW, dV2, dI, dR))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sirv_booster_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", V1 = "#2457C5", W = "#A5D6A7", V2 = "#2457C5", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "vaccination with waning and booster doses")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
