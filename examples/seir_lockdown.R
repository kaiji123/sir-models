## SEIR with a temporary lockdown
## Contacts fall by a proportion eff between days t_on and t_off.
## Reference: Flaxman et al. (2020), Nature 584:257-261 (concept)
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
times <- seq(0, 365, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  R0 = 2.5,       # basic reproduction number
  sigma = 1 / 5,  # 1/(latent period)
  gamma = 1 / 5,  # 1/(infectious period)
  eff = 0.6,      # proportional reduction in contacts during lockdown
  t_on = 50,      # day lockdown starts
  t_off = 120     # day lockdown ends
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 999990,  # Susceptible
  E = 10,      # Exposed
  I = 0,       # Infectious
  R = 0        # Recovered
)

# ---- Model equations -------------------------------------------------
seir_lockdown_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + I + R

    # Derived quantities
    contact <- 1 - eff * (t >= t_on & t < t_off)  # relative contact rate
    lam <- contact * R0 * gamma * I / N           # force of infection
    Rt <- contact * R0 * S / N                    # effective reproduction number

    # Flows recorded as outputs
    incidence <- sigma * E

    # Rates of change
    dS <- -lam * S
    dE <- lam * S - sigma * E
    dI <- sigma * E - gamma * I
    dR <- gamma * I

    list(c(dS, dE, dI, dR), c(Rt = Rt, incidence = incidence))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = seir_lockdown_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SEIR with a temporary lockdown")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("Rt", "incidence")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "SEIR with a temporary lockdown: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
