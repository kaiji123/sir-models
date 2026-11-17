## SEIAR (asymptomatic infection)
## A proportion of infections is asymptomatic (A) and less infectious.
## Reference: Arino et al. (2008), J. R. Soc. Interface 5(Suppl 2)
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, I, A, R
## Flows (time unit: days):
##   S -> E  : lam * S
##   E -> I  : sigma * ps * E
##   E -> A  : sigma * (1 - ps) * E
##   I -> R  : gamma * I
##   A -> R  : gamma * A

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.5,     # transmission rate per day
  rA = 0.5,       # relative infectiousness of asymptomatics
  sigma = 1 / 4,  # 1/(latent period)
  ps = 0.6,       # proportion of infections that are symptomatic
  gamma = 1 / 6   # recovery rate
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 99990,  # Susceptible
  E = 10,     # Exposed
  I = 0,      # Symptomatic
  A = 0,      # Asymptomatic
  R = 0       # Recovered
)

# ---- Model equations -------------------------------------------------
seiar_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + I + A + R

    # Derived quantities
    lam <- beta * (I + rA * A) / N  # force of infection

    # Flows recorded as outputs
    cases <- sigma * ps * E

    # Rates of change
    dS <- -lam * S
    dE <- lam * S - sigma * ps * E - sigma * (1 - ps) * E
    dI <- sigma * ps * E - gamma * I
    dA <- sigma * (1 - ps) * E - gamma * A
    dR <- gamma * I + gamma * A

    list(c(dS, dE, dI, dA, dR), c(cases = cases))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = seiar_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", I = "#B71C1C", A = "#E57373", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SEIAR (asymptomatic infection)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("cases")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "SEIAR (asymptomatic infection): outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
