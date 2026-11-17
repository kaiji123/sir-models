## SIR with births and deaths
## Endemic SIR with demography; damped oscillations to an endemic equilibrium (measles-like parameters).
## Reference: Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals, §2.1.2 (Program 2.2)
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, R
## Flows (time unit: years):
##   (outside) -> S  : mu * N
##   S -> I          : lam * S
##   I -> R          : gamma * I
##   S -> (outside)  : mu * S
##   I -> (outside)  : mu * I
##   R -> (outside)  : mu * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 60, by = 0.02)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 520,       # transmission rate per year
  gamma = 365 / 7,  # recovery rate per year (7-day infectious period)
  mu = 1 / 70       # birth = death rate per year (70-year life expectancy)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 10000,  # Susceptible
  I = 10,     # Infectious
  R = 89990   # Recovered
)

# ---- Model equations -------------------------------------------------
sir_vital_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + R

    # Derived quantities
    lam <- beta * I / N  # force of infection

    # Flows recorded as outputs
    incidence <- lam * S

    # Rates of change
    dS <- mu * N - lam * S - mu * S
    dI <- lam * S - gamma * I - mu * I
    dR <- gamma * I - mu * R

    list(c(dS, dI, dR), c(incidence = incidence))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sir_vital_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "SIR with births and deaths")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("incidence")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (years)", ylab = "Value",
        main = "SIR with births and deaths: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
