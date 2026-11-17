## SEITR influenza-type model with antiviral treatment
## pandemic/seasonal influenza, antiviral treatment.
## Closed population, so N is constant.
## Parameter values: 0 of 6 checked against the source; the rest are recalled or illustrative.
## Reference: Brauer F. (2008) "Compartmental models in epidemiology", in Brauer, van den Driessche & Wu (eds), Mathematical Epidemiology, Lecture Notes in Mathematics 1945, Springer, pp. 19–79
## Source: https://doi.org/10.1007/978-3-540-78911-6_2
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, I, Tr, R
## Flows (time unit: days):
##   S -> E   : lambda * S
##   E -> I   : kappa * E
##   I -> Tr  : tau * I
##   I -> R   : alpha * I
##   Tr -> R  : eta * Tr

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 0.2)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.5,       # transmission rate (1/day) [illustrative]
  kappa = 1 / 1.5,  # progression from latent (1/day) [illustrative]
  tau = 0.3,        # treatment uptake rate (1/day) [illustrative]
  alpha = 1 / 4,    # untreated recovery rate (1/day) [illustrative]
  eta = 1 / 2.5,    # recovery rate under treatment (1/day) [illustrative]
  delta = 0.5       # relative infectiousness of treated people [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9990,  # susceptible
  E = 0,     # latent
  I = 10,    # infectious, untreated
  Tr = 0,    # infectious, under treatment
  R = 0      # recovered
)

# ---- Model equations -------------------------------------------------
sitr_treatment_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + I + Tr + R

    # Derived quantities
    lambda <- beta * (I + delta * Tr) / N                                   # force of infection
    R0 <- beta * (1 / (alpha + tau) + delta * tau / ((alpha + tau) * eta))  # basic reproduction number (about 1.26 with these values; set beta = 0.8 for about 2.0)

    # Rates of change
    dS <- -lambda * S
    dE <- lambda * S - kappa * E
    dI <- kappa * E - tau * I - alpha * I
    dTr <- tau * I - eta * Tr
    dR <- alpha * I + eta * Tr

    list(c(dS, dE, dI, dTr, dR))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sitr_treatment_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", I = "#B71C1C", Tr = "#333333", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SEITR influenza-type model with antiviral treatment")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
