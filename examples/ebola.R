## Ebola (community, hospital, funeral)
## Transmission in the community (I), in hospital (H) and at funerals (F); structure after Legrand et al.
## Parameter values are illustrative.
## Reference: Legrand et al. (2007), Epidemiology and Infection 135:610-621
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, I, H, Fu, R
## Flows (time unit: days):
##   S -> E           : lam * S
##   E -> I           : alpha * E
##   I -> H           : gammaH * theta * I
##   I -> Fu          : gammaD * (1 - theta) * delta * I
##   I -> R           : gammaI * (1 - theta) * (1 - delta) * I
##   H -> Fu          : gammaDH * delta * H
##   H -> R           : gammaIH * (1 - delta) * H
##   Fu -> (outside)  : gammaF * Fu

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 300, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  betaI = 0.3,        # transmission rate in the community
  betaH = 0.05,       # transmission rate in hospital
  betaF = 0.4,        # transmission rate at funerals
  alpha = 1 / 7,      # 1/(incubation period)
  theta = 0.5,        # proportion of cases hospitalised
  delta = 0.7,        # case fatality ratio
  gammaH = 1 / 5,     # 1/(onset to hospitalisation)
  gammaD = 1 / 9.6,   # 1/(onset to death)
  gammaI = 1 / 10,    # 1/(onset to recovery)
  gammaDH = 1 / 4.6,  # 1/(hospitalisation to death)
  gammaIH = 1 / 5,    # 1/(hospitalisation to recovery)
  gammaF = 1 / 2      # 1/(death to burial)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 99999,  # Susceptible
  E = 1,      # Exposed
  I = 0,      # Infectious (community)
  H = 0,      # Hospitalised
  Fu = 0,     # Dead, not yet buried (funeral)
  R = 0       # Removed
)

# ---- Model equations -------------------------------------------------
ebola_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + I + H + Fu + R

    # Derived quantities
    lam <- (betaI * I + betaH * H + betaF * Fu) / N  # force of infection

    # Flows recorded as outputs
    onsets <- alpha * E
    burials <- gammaF * Fu

    # Rates of change
    dS <- -lam * S
    dE <- lam * S - alpha * E
    dI <- alpha * E - gammaH * theta * I - gammaD * (1 - theta) * delta * I - gammaI * (1 - theta) * (1 - delta) * I
    dH <- gammaH * theta * I - gammaDH * delta * H - gammaIH * (1 - delta) * H
    dFu <- gammaD * (1 - theta) * delta * I + gammaDH * delta * H - gammaF * Fu
    dR <- gammaI * (1 - theta) * (1 - delta) * I + gammaIH * (1 - delta) * H

    list(c(dS, dE, dI, dH, dFu, dR), c(onsets = onsets, burials = burials))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = ebola_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", I = "#B71C1C", H = "#6A2C9C", Fu = "#5D4037", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Ebola (community, hospital, funeral)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("onsets", "burials")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "Ebola (community, hospital, funeral): outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
