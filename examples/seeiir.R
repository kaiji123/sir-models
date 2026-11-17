## SEEIIR (gamma-distributed periods)
## Two latent and two infectious stages give gamma-distributed (Erlang k = 2) latent and infectious
## periods with the same means as the SEIR model ("linear chain trick").
## Reference: Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals, §3.3
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E1, E2, I1, I2, R
## Flows (time unit: days):
##   S -> E1   : lam * S
##   E1 -> E2  : 2 * sigma * E1
##   E2 -> I1  : 2 * sigma * E2
##   I1 -> I2  : 2 * gamma * I1
##   I2 -> R   : 2 * gamma * I2

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.5,     # transmission rate per day
  sigma = 1 / 5,  # 1/(mean latent period)
  gamma = 1 / 7   # 1/(mean infectious period)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9999,  # Susceptible
  E1 = 1,    # Latent, stage 1
  E2 = 0,    # Latent, stage 2
  I1 = 0,    # Infectious, stage 1
  I2 = 0,    # Infectious, stage 2
  R = 0      # Recovered
)

# ---- Model equations -------------------------------------------------
seeiir_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E1 + E2 + I1 + I2 + R

    # Derived quantities
    lam <- beta * (I1 + I2) / N  # force of infection
    prevalence <- I1 + I2        # number infectious

    # Flows recorded as outputs
    onsets <- 2 * sigma * E2

    # Rates of change
    dS <- -lam * S
    dE1 <- lam * S - 2 * sigma * E1
    dE2 <- 2 * sigma * E1 - 2 * sigma * E2
    dI1 <- 2 * sigma * E2 - 2 * gamma * I1
    dI2 <- 2 * gamma * I1 - 2 * gamma * I2
    dR <- 2 * gamma * I2

    list(c(dS, dE1, dE2, dI1, dI2, dR), c(prevalence = prevalence, onsets = onsets))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = seeiir_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E1 = "#E8C468", E2 = "#E8C468", I1 = "#B71C1C", I2 = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SEEIIR (gamma-distributed periods)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("prevalence", "onsets")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "SEEIIR (gamma-distributed periods): outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
