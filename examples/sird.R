## SIRD
## SIR with disease-induced deaths. D is excluded from the living population N.
## Reference: Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals, §2.2
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, R, D
## Flows (time unit: days):
##   S -> I  : lam * S
##   I -> R  : gamma * (1 - cfr) * I
##   I -> D  : gamma * cfr * I

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.3,      # transmission rate per day
  gamma = 1 / 10,  # rate of leaving the infectious state
  cfr = 0.02       # case fatality ratio
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9990,  # Susceptible
  I = 10,    # Infectious
  R = 0,     # Recovered
  D = 0      # Dead
)

# ---- Model equations -------------------------------------------------
sird_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + R

    # Derived quantities
    lam <- beta * I / N  # force of infection

    # Flows recorded as outputs
    deaths <- gamma * cfr * I

    # Rates of change
    dS <- -lam * S
    dI <- lam * S - gamma * (1 - cfr) * I - gamma * cfr * I
    dR <- gamma * (1 - cfr) * I
    dD <- gamma * cfr * I

    list(c(dS, dI, dR, dD), c(deaths = deaths))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sird_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", R = "#4E8A3E", D = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SIRD")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("deaths")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "SIRD: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
