## SEIR
## Epidemic with a latent (infected, not yet infectious) period, closed population.
## Reference: Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals, §2.5
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, I, R
## Flows (time unit: days):
##   S -> E  : lam * S
##   E -> I  : sigma * E
##   I -> R  : gamma * I

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.5,     # transmission rate per day
  sigma = 1 / 5,  # rate of becoming infectious = 1/(latent period)
  gamma = 1 / 7   # recovery rate = 1/(infectious period)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9999,  # Susceptible
  E = 1,     # Exposed (latent)
  I = 0,     # Infectious
  R = 0      # Recovered
)

# ---- Model equations -------------------------------------------------
seir_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + I + R

    # Derived quantities
    lam <- beta * I / N  # force of infection

    # Flows recorded as outputs
    onsets <- sigma * E

    # Rates of change
    dS <- -lam * S
    dE <- lam * S - sigma * E
    dI <- sigma * E - gamma * I
    dR <- gamma * I

    list(c(dS, dE, dI, dR), c(onsets = onsets))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = seir_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SEIR")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("onsets")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "SEIR: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
