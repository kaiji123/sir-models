## SEIRD
## SEIR with disease-induced deaths.
## Reference: Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals, §2.5
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, I, R, D
## Flows (time unit: days):
##   S -> E  : lam * S
##   E -> I  : sigma * E
##   I -> R  : gamma * (1 - cfr) * I
##   I -> D  : gamma * cfr * I

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 250, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.4,     # transmission rate per day
  sigma = 1 / 5,  # 1/(latent period)
  gamma = 1 / 8,  # rate of leaving the infectious state
  cfr = 0.01      # case fatality ratio
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 99990,  # Susceptible
  E = 10,     # Exposed
  I = 0,      # Infectious
  R = 0,      # Recovered
  D = 0       # Dead
)

# ---- Model equations -------------------------------------------------
seird_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + I + R

    # Derived quantities
    lam <- beta * I / N  # force of infection

    # Flows recorded as outputs
    incidence <- sigma * E
    deaths <- gamma * cfr * I

    # Rates of change
    dS <- -lam * S
    dE <- lam * S - sigma * E
    dI <- sigma * E - gamma * (1 - cfr) * I - gamma * cfr * I
    dR <- gamma * (1 - cfr) * I
    dD <- gamma * cfr * I

    list(c(dS, dE, dI, dR, dD), c(incidence = incidence, deaths = deaths))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = seird_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", I = "#B71C1C", R = "#4E8A3E", D = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SEIRD")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("incidence", "deaths")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "SEIRD: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
