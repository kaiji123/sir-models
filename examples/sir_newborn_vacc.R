## SIR with newborn vaccination
## A proportion p of newborns is vaccinated (enters V). Elimination requires p > 1 - 1/R0.
## Reference: Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals, §8.1; Anderson & May (1991), Infectious Diseases of Humans, ch. 5
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, R, V
## Flows (time unit: years):
##   (outside) -> S  : (1 - p) * mu * N
##   (outside) -> V  : p * mu * N
##   S -> I          : lam * S
##   I -> R          : gamma * I
##   S -> (outside)  : mu * S
##   I -> (outside)  : mu * I
##   R -> (outside)  : mu * R
##   V -> (outside)  : mu * V

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 80, by = 0.02)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  R0 = 15,          # basic reproduction number
  gamma = 365 / 7,  # recovery rate per year
  mu = 1 / 70,      # birth = death rate per year
  p = 0.8           # proportion of newborns vaccinated
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 10000,  # Susceptible
  I = 10,     # Infectious
  R = 89990,  # Recovered
  V = 0       # Vaccinated
)

# ---- Model equations -------------------------------------------------
sir_newborn_vacc_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + R + V

    # Derived quantities
    beta <- R0 * (gamma + mu)  # transmission rate implied by R0
    lam <- beta * I / N        # force of infection

    # Flows recorded as outputs
    incidence <- lam * S

    # Rates of change
    dS <- (1 - p) * mu * N - lam * S - mu * S
    dI <- lam * S - gamma * I - mu * I
    dR <- gamma * I - mu * R
    dV <- p * mu * N - mu * V

    list(c(dS, dI, dR, dV), c(incidence = incidence))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sir_newborn_vacc_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", R = "#4E8A3E", V = "#2457C5")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "SIR with newborn vaccination")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("incidence")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (years)", ylab = "Value",
        main = "SIR with newborn vaccination: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
