## SVIR with a leaky vaccine
## Vaccinated individuals are infected at a reduced rate (1 - eff) * lam.
## Reference: Vynnycky & White (2010), An Introduction to Infectious Disease Modelling, ch. 7
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, V, I, R
## Flows (time unit: days):
##   S -> I  : lam * S
##   V -> I  : (1 - eff) * lam * V
##   I -> R  : gamma * I

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 300, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.4,     # transmission rate per day
  gamma = 1 / 7,  # recovery rate
  eff = 0.6       # vaccine efficacy against infection
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 3999,  # Susceptible
  V = 6000,  # Vaccinated
  I = 1,     # Infectious
  R = 0      # Recovered
)

# ---- Model equations -------------------------------------------------
leaky_vaccine_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + V + I + R

    # Derived quantities
    lam <- beta * I / N  # force of infection

    # Flows recorded as outputs
    breakthrough <- (1 - eff) * lam * V

    # Rates of change
    dS <- -lam * S
    dV <- -(1 - eff) * lam * V
    dI <- lam * S + (1 - eff) * lam * V - gamma * I
    dR <- gamma * I

    list(c(dS, dV, dI, dR), c(breakthrough = breakthrough))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = leaky_vaccine_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", V = "#2457C5", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SVIR with a leaky vaccine")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("breakthrough")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "SVIR with a leaky vaccine: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
