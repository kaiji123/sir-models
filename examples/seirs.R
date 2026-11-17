## SEIRS
## SEIR with waning immunity, giving recurrent epidemics.
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
##   R -> S  : omega * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 1500, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.5,      # transmission rate per day
  sigma = 1 / 5,   # 1/(latent period)
  gamma = 1 / 7,   # 1/(infectious period)
  omega = 1 / 365  # 1/(duration of immunity)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9999,  # Susceptible
  E = 1,     # Exposed
  I = 0,     # Infectious
  R = 0      # Recovered
)

# ---- Model equations -------------------------------------------------
seirs_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + I + R

    # Derived quantities
    lam <- beta * I / N  # force of infection

    # Rates of change
    dS <- -lam * S + omega * R
    dE <- lam * S - sigma * E
    dI <- sigma * E - gamma * I
    dR <- gamma * I - omega * R

    list(c(dS, dE, dI, dR))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = seirs_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SEIRS")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
