## SEIHRD (hospitalisation)
## A proportion of cases is hospitalised; deaths occur in hospital.
## Reference: Vynnycky & White (2010), An Introduction to Infectious Disease Modelling, ch. 4 (structure)
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, I, R, H, D
## Flows (time unit: days):
##   S -> E  : lam * S
##   E -> I  : sigma * E
##   I -> R  : gamma * (1 - h) * I
##   I -> H  : gamma * h * I
##   H -> R  : tauH * (1 - m) * H
##   H -> D  : tauH * m * H

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 300, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.45,    # transmission rate per day
  sigma = 1 / 5,  # 1/(latent period)
  gamma = 1 / 7,  # rate of leaving the infectious state
  h = 0.05,       # proportion hospitalised
  tauH = 1 / 10,  # rate of leaving hospital
  m = 0.2         # in-hospital fatality
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 99990,  # Susceptible
  E = 10,     # Exposed
  I = 0,      # Infectious
  R = 0,      # Recovered
  H = 0,      # Hospitalised
  D = 0       # Dead
)

# ---- Model equations -------------------------------------------------
seihrd_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + I + R + H

    # Derived quantities
    lam <- beta * I / N  # force of infection (hospitalised patients do not transmit)

    # Flows recorded as outputs
    admissions <- gamma * h * I
    deaths <- tauH * m * H

    # Rates of change
    dS <- -lam * S
    dE <- lam * S - sigma * E
    dI <- sigma * E - gamma * (1 - h) * I - gamma * h * I
    dR <- gamma * (1 - h) * I + tauH * (1 - m) * H
    dH <- gamma * h * I - tauH * (1 - m) * H - tauH * m * H
    dD <- tauH * m * H

    list(c(dS, dE, dI, dR, dH, dD), c(admissions = admissions, deaths = deaths))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = seihrd_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", I = "#B71C1C", R = "#4E8A3E", H = "#6A2C9C", D = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SEIHRD (hospitalisation)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("admissions", "deaths")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "SEIHRD (hospitalisation): outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
