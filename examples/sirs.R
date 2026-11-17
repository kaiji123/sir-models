## SIRS
## SIR with waning immunity: recovered individuals return to the susceptible class.
## Reference: Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals, §2.4
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, R
## Flows (time unit: days):
##   S -> I  : lam * S
##   I -> R  : gamma * I
##   R -> S  : omega * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 1000, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.3,      # transmission rate per day
  gamma = 1 / 10,  # recovery rate
  omega = 1 / 180  # rate of loss of immunity = 1/(duration of immunity)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 990,  # Susceptible
  I = 10,   # Infectious
  R = 0     # Recovered
)

# ---- Model equations -------------------------------------------------
sirs_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + R

    # Derived quantities
    lam <- beta * I / N  # force of infection

    # Rates of change
    dS <- -lam * S + omega * R
    dI <- lam * S - gamma * I
    dR <- gamma * I - omega * R

    list(c(dS, dI, dR))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sirs_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SIRS")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
