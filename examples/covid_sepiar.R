## COVID-19 style SEPIR
## Latent (E), pre-symptomatic infectious (P), then symptomatic (Is) or asymptomatic (Ia) infection.
## Reference: Davies et al. (2020), Nature Medicine 26:1205-1211 (structure)
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, P, Is, Ia, R
## Flows (time unit: days):
##   S -> E   : lam * S
##   E -> P   : sigma * E
##   P -> Is  : delta * ps * P
##   P -> Ia  : delta * (1 - ps) * P
##   Is -> R  : gamma * Is
##   Ia -> R  : gamma * Ia

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 300, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  R0 = 2.5,       # basic reproduction number
  sigma = 1 / 3,  # 1/(latent period)
  delta = 1 / 2,  # 1/(pre-symptomatic period)
  gamma = 1 / 5,  # 1/(infectious period after the pre-symptomatic stage)
  ps = 0.6,       # proportion symptomatic
  rP = 1,         # relative infectiousness, pre-symptomatic
  rA = 0.5        # relative infectiousness, asymptomatic
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 999990,  # Susceptible
  E = 10,      # Latent
  P = 0,       # Pre-symptomatic
  Is = 0,      # Symptomatic
  Ia = 0,      # Asymptomatic
  R = 0        # Recovered
)

# ---- Model equations -------------------------------------------------
covid_sepiar_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + P + Is + Ia + R

    # Derived quantities
    beta <- R0 / (rP / delta + (ps + (1 - ps) * rA) / gamma)  # transmission rate implied by R0
    lam <- beta * (rP * P + Is + rA * Ia) / N                 # force of infection

    # Flows recorded as outputs
    onsets <- delta * ps * P

    # Rates of change
    dS <- -lam * S
    dE <- lam * S - sigma * E
    dP <- sigma * E - delta * ps * P - delta * (1 - ps) * P
    dIs <- delta * ps * P - gamma * Is
    dIa <- delta * (1 - ps) * P - gamma * Ia
    dR <- gamma * Is + gamma * Ia

    list(c(dS, dE, dP, dIs, dIa, dR), c(onsets = onsets))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = covid_sepiar_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", P = "#F0A35E", Is = "#B71C1C", Ia = "#E57373", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "COVID-19 style SEPIR")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("onsets")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "COVID-19 style SEPIR: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
