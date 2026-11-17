## SEIQHR quarantine and isolation (Gumel et al. 2004, SARS)
## SARS 2003 (Toronto, Hong Kong, Singapore, Beijing).
## The headline result is that optimal isolation of symptomatic cases is more effective than sub-optimal isolation combined with quarantine (search summary of the abstract).
## Parameter values: 0 of 14 checked against the source; the rest are recalled or illustrative.
## Reference: Gumel A.B. et al. (2004) "Modelling strategies for controlling SARS outbreaks", Proc. R. Soc. B 271:2223–2232
## Source: https://doi.org/10.1098/rspb.2004.2800
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, Q, I, H, R, Dd
## Flows (time unit: days):
##   (outside) -> S  : Pi_b
##   S -> E          : lambda * S
##   E -> Q          : gamma1 * E
##   E -> I          : kappa1 * E
##   Q -> H          : kappa2 * Q
##   I -> H          : gamma2 * I
##   I -> R          : sigma1 * I
##   H -> R          : sigma2 * H
##   I -> Dd         : d1 * I
##   H -> Dd         : d2 * H
##   S -> (outside)  : mu * S
##   E -> (outside)  : mu * E
##   Q -> (outside)  : mu * Q
##   I -> (outside)  : mu * I
##   H -> (outside)  : mu * H
##   R -> (outside)  : mu * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 0.2)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  Pi_b = 136,           # recruitment (recalled; or set to muN0) (persons/day)
  mu = 1 / (70 * 365),  # natural death (1/day) [illustrative]
  beta = 0.2,           # transmission rate (1/day) [recalled]
  eps_E = 0.3,          # relative infectiousness of E [recalled]
  eps_Q = 0,            # relative infectiousness of Q [recalled]
  eps_H = 0.36,         # relative infectiousness of the isolated (Toronto hospital transmission) [recalled]
  gamma1 = 0.1,         # quarantine rate of exposed (1/day) [recalled]
  gamma2 = 0.5,         # isolation rate of symptomatic (1/day) [recalled]
  kappa1 = 0.1,         # progression E→I (1/day) [recalled]
  kappa2 = 0.125,       # progression Q→H (1/day) [recalled]
  sigma1 = 0.0337,      # recovery of I (1/day) [recalled]
  sigma2 = 0.0386,      # recovery of H (1/day) [recalled]
  d1 = 0.0068,          # SARS death rate in I (1/day) [recalled]
  d2 = 0.0045           # SARS death rate in H (1/day) [recalled]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 4000000,  # susceptible
  E = 10,       # exposed/asymptomatic, mildly infectious
  Q = 0,        # quarantined exposed
  I = 1,        # symptomatic infectious
  H = 0,        # isolated/hospitalised (Gumel's J)
  R = 0,        # recovered
  Dd = 0        # cumulative SARS deaths
)

# ---- Model equations -------------------------------------------------
seiqhr_sars_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + Q + I + H + R

    # Derived quantities
    lambda <- beta * (I + eps_E * E + eps_Q * Q + eps_H * H) / N

    # Rates of change
    dS <- Pi_b - lambda * S - mu * S
    dE <- lambda * S - gamma1 * E - kappa1 * E - mu * E
    dQ <- gamma1 * E - kappa2 * Q - mu * Q
    dI <- kappa1 * E - gamma2 * I - sigma1 * I - d1 * I - mu * I
    dH <- kappa2 * Q + gamma2 * I - sigma2 * H - d2 * H - mu * H
    dR <- sigma1 * I + sigma2 * H - mu * R
    dDd <- d1 * I + d2 * H

    list(c(dS, dE, dQ, dI, dH, dR, dDd))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = seiqhr_sars_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", Q = "#00838F", I = "#B71C1C", H = "#6A2C9C", R = "#4E8A3E", Dd = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SEIQHR quarantine and isolation (Gumel et al. 2004, SARS)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
