## SIQR (isolation)
## Infectious individuals are detected and isolated (Q) at rate q and no longer transmit.
## Reference: Hethcote, Zhien & Shengbing (2002), Math. Biosci. 180:141-160
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, Q, R
## Flows (time unit: days):
##   S -> I  : lam * S
##   I -> R  : gamma * I
##   I -> Q  : q * I
##   Q -> R  : gamma * Q

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.4,     # transmission rate per day
  gamma = 1 / 7,  # recovery rate
  q = 0.1         # rate of detection and isolation
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9990,  # Susceptible
  I = 10,    # Infectious
  Q = 0,     # Isolated
  R = 0      # Recovered
)

# ---- Model equations -------------------------------------------------
siqr_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + Q + R

    # Derived quantities
    lam <- beta * I / N  # force of infection

    # Rates of change
    dS <- -lam * S
    dI <- lam * S - gamma * I - q * I
    dQ <- q * I - gamma * Q
    dR <- gamma * I + gamma * Q

    list(c(dS, dI, dQ, dR))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = siqr_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", Q = "#00838F", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SIQR (isolation)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
