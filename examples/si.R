## SI
## Infection with no recovery (e.g. untreated chronic infection) in a closed population.
## Reference: Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals, ch. 2
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I
## Flows (time unit: days):
##   S -> I  : lam * S

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 100, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.2  # transmission rate per day
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 999,  # Susceptible
  I = 1     # Infected
)

# ---- Model equations -------------------------------------------------
si_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I

    # Derived quantities
    lam <- beta * I / N  # force of infection

    # Flows recorded as outputs
    incidence <- lam * S

    # Rates of change
    dS <- -lam * S
    dI <- lam * S

    list(c(dS, dI), c(incidence = incidence))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = si_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SI")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("incidence")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "SI: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
