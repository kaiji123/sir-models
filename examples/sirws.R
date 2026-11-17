## SIRWS (waning and boosting)
## Immunity wanes through a W stage; re-exposure while in W boosts immunity back to R (pertussis-like).
## Reference: Lavine, King & Bjornstad (2011), PNAS 108:7259-7264
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, R, W
## Flows (time unit: years):
##   (outside) -> S  : mu * N
##   S -> I          : lam * S
##   I -> R          : gamma * I
##   R -> W          : 2 * kappa * R
##   W -> S          : 2 * kappa * W
##   W -> R          : boost * lam * W
##   S -> (outside)  : mu * S
##   I -> (outside)  : mu * I
##   R -> (outside)  : mu * R
##   W -> (outside)  : mu * W

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 100, by = 0.02)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 260,        # transmission rate per year
  gamma = 365 / 21,  # recovery rate per year (21-day infectious period)
  kappa = 1 / 10,    # rate of waning per year (two stages, so mean 2 / (2 kappa) = 10 years)
  boost = 20,        # boosting factor relative to the force of infection
  mu = 1 / 70        # birth = death rate per year
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 5000,   # Susceptible
  I = 10,     # Infectious
  R = 70000,  # Immune
  W = 24990   # Waning immunity
)

# ---- Model equations -------------------------------------------------
sirws_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + R + W

    # Derived quantities
    lam <- beta * I / N  # force of infection

    # Flows recorded as outputs
    incidence <- lam * S

    # Rates of change
    dS <- mu * N - lam * S + 2 * kappa * W - mu * S
    dI <- lam * S - gamma * I - mu * I
    dR <- gamma * I - 2 * kappa * R + boost * lam * W - mu * R
    dW <- 2 * kappa * R - 2 * kappa * W - boost * lam * W - mu * W

    list(c(dS, dI, dR, dW), c(incidence = incidence))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sirws_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", R = "#4E8A3E", W = "#A5D6A7")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "SIRWS (waning and boosting)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("incidence")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (years)", ylab = "Value",
        main = "SIRWS (waning and boosting): outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
