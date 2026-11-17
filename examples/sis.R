## SIS
## Infection without lasting immunity (e.g. many bacterial STIs): recovered people become susceptible again.
## Reference: Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals, §2.4
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I
## Flows (time unit: days):
##   S -> I  : lam * S
##   I -> S  : gamma * I

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.3,     # transmission rate per day
  gamma = 1 / 10  # recovery rate = 1/(infectious period)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 999,  # Susceptible
  I = 1     # Infectious
)

# ---- Model equations -------------------------------------------------
sis_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I

    # Derived quantities
    lam <- beta * I / N  # force of infection

    # Rates of change
    dS <- -lam * S + gamma * I
    dI <- lam * S - gamma * I

    list(c(dS, dI))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sis_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SIS")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
