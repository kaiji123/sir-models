## SEACHR with re-exposure
## Asymptomatic (A) and clinical (C) infection, hospitalisation (H), death (D),
## waning immunity and re-exposure of recovereds (ER). Parameter values are illustrative;
## the force of infection is an assumption (the diagram does not define it).
## Reference: Custom diagram
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, A, C, R, ER, H, D
## Flows (time unit: days):
##   S -> E   : lam * S
##   E -> A   : nu * (1 - gamma) * E
##   E -> C   : nu * gamma * E
##   A -> R   : tauA * A
##   C -> R   : tauC * (1 - theta) * C
##   C -> H   : tauC * theta * C
##   H -> R   : tauH * (1 - omega) * H
##   H -> D   : tauH * omega * H
##   R -> ER  : sigma * lam * R
##   ER -> A  : nu * (1 - gamma / 2) * ER
##   ER -> C  : nu * gamma / 2 * ER
##   R -> S   : alpha * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 365, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.35,     # transmission rate per day
  rA = 0.5,        # relative infectiousness of asymptomatics
  nu = 1 / 9,      # 1/(latent period)
  gamma = 0.7,     # proportion of infections that are clinical
  tauA = 1 / 7,    # recovery rate of asymptomatics
  tauC = 1 / 5,    # rate of leaving the clinical stage
  theta = 0.3,     # proportion of clinical cases hospitalised
  tauH = 1 / 10,   # rate of leaving hospital
  omega = 0.4,     # case fatality among hospitalised
  sigma = 0.2,     # relative susceptibility of recovereds
  alpha = 1 / 365  # rate of waning immunity
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  E <- 10          # Exposed
  S <- 100000 - E  # Susceptible
  A <- 0           # Asymptomatic
  C <- 0           # Clinical
  R <- 0           # Recovered
  ER <- 0          # Exposed recovered
  H <- 0           # Hospitalised
  D <- 0           # Dead
  c(S = S, E = E, A = A, C = C, R = R, ER = ER, H = H, D = D)
})

# ---- Model equations -------------------------------------------------
seachr_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + A + C + R + ER + H

    # Derived quantities
    lam <- beta * (rA * A + C) / N  # force of infection (assumed)

    # Flows recorded as outputs
    admissions <- tauC * theta * C
    deaths <- tauH * omega * H

    # Rates of change
    dS <- -lam * S + alpha * R
    dE <- lam * S - nu * (1 - gamma) * E - nu * gamma * E
    dA <- nu * (1 - gamma) * E - tauA * A + nu * (1 - gamma / 2) * ER
    dC <- nu * gamma * E - tauC * (1 - theta) * C - tauC * theta * C + nu * gamma / 2 * ER
    dR <- tauA * A + tauC * (1 - theta) * C + tauH * (1 - omega) * H - sigma * lam * R - alpha * R
    dER <- sigma * lam * R - nu * (1 - gamma / 2) * ER - nu * gamma / 2 * ER
    dH <- tauC * theta * C - tauH * (1 - omega) * H - tauH * omega * H
    dD <- tauH * omega * H

    list(c(dS, dE, dA, dC, dR, dER, dH, dD), c(admissions = admissions, deaths = deaths))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = seachr_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", A = "#C62828", C = "#8E0E0E", R = "#4E8A3E", ER = "#D9A441", H = "#6A2C9C", D = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SEACHR with re-exposure")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("admissions", "deaths")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "SEACHR with re-exposure: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
