## SEIR with births and deaths
## Endemic SEIR parameterised by R0 (measles-like): transmission rate derived from R0.
## Reference: Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals, §2.5 (Program 2.6)
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, I, R
## Flows (time unit: years):
##   (outside) -> S  : mu * N
##   S -> E          : lam * S
##   E -> I          : sigma * E
##   I -> R          : gamma * I
##   S -> (outside)  : mu * S
##   E -> (outside)  : mu * E
##   I -> (outside)  : mu * I
##   R -> (outside)  : mu * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 50, by = 0.01)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  R0 = 13,          # basic reproduction number
  sigma = 365 / 8,  # rate of becoming infectious per year (8-day latent period)
  gamma = 365 / 5,  # recovery rate per year (5-day infectious period)
  mu = 1 / 70       # birth = death rate per year
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 6000,  # Susceptible
  E = 10,    # Exposed
  I = 10,    # Infectious
  R = 93980  # Recovered
)

# ---- Model equations -------------------------------------------------
seir_vital_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + I + R

    # Derived quantities
    beta <- R0 * (gamma + mu) * (sigma + mu) / sigma  # transmission rate implied by R0
    lam <- beta * I / N                               # force of infection

    # Flows recorded as outputs
    incidence <- sigma * E

    # Rates of change
    dS <- mu * N - lam * S - mu * S
    dE <- lam * S - sigma * E - mu * E
    dI <- sigma * E - gamma * I - mu * I
    dR <- gamma * I - mu * R

    list(c(dS, dE, dI, dR), c(incidence = incidence))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = seir_vital_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "SEIR with births and deaths")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("incidence")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (years)", ylab = "Value",
        main = "SEIR with births and deaths: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
