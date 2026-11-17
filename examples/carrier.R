## SICR (chronic carriers)
## A proportion q of acute infections become long-term carriers (hepatitis B / typhoid-like).
## Reference: Anderson & May (1991), Infectious Diseases of Humans, ch. 6; Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals, §3.2
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, C, R
## Flows (time unit: years):
##   (outside) -> S  : mu * N
##   S -> I          : lam * S
##   I -> R          : gamma * (1 - q) * I
##   I -> C          : gamma * q * I
##   C -> R          : gammaC * C
##   S -> (outside)  : mu * S
##   I -> (outside)  : mu * I
##   C -> (outside)  : mu * C
##   R -> (outside)  : mu * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 150, by = 0.1)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 4,        # transmission rate per year
  eps = 0.1,       # relative infectiousness of carriers
  gamma = 6,       # rate of leaving acute infection per year (2 months)
  q = 0.1,         # proportion becoming carriers
  gammaC = 0.025,  # rate of carriers clearing infection per year
  mu = 1 / 70      # birth = death rate per year
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 99000,  # Susceptible
  I = 10,     # Acute infection
  C = 990,    # Chronic carrier
  R = 0       # Recovered
)

# ---- Model equations -------------------------------------------------
carrier_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + C + R

    # Derived quantities
    lam <- beta * (I + eps * C) / N  # force of infection

    # Flows recorded as outputs
    incidence <- lam * S

    # Rates of change
    dS <- mu * N - lam * S - mu * S
    dI <- lam * S - gamma * (1 - q) * I - gamma * q * I - mu * I
    dC <- gamma * q * I - gammaC * C - mu * C
    dR <- gamma * (1 - q) * I + gammaC * C - mu * R

    list(c(dS, dI, dC, dR), c(incidence = incidence))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = carrier_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", C = "#8E0E0E", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "SICR (chronic carriers)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("incidence")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (years)", ylab = "Value",
        main = "SICR (chronic carriers): outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
