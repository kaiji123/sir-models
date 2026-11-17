## SIR
## Basic SIR epidemic in a closed population.
## Reference: Kermack & McKendrick (1927); Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals, §2.1
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, R
## Flows (time unit: days):
##   S -> I  : lam * S
##   I -> R  : gamma * I

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 160, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.3,     # transmission rate per day
  gamma = 1 / 10  # recovery rate = 1/(infectious period)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 990,  # Susceptible
  I = 10,   # Infectious
  R = 0     # Recovered
)

# ---- Model equations -------------------------------------------------
sir_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + R

    # Derived quantities
    lam <- beta * I / N  # force of infection

    # Flows recorded as outputs
    incidence <- lam * S

    # Rates of change
    dS <- -lam * S
    dI <- lam * S - gamma * I
    dR <- gamma * I

    list(c(dS, dI, dR), c(incidence = incidence))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sir_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SIR")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("incidence")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "SIR: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
