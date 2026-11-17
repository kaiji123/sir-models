## SIR with aware (behaviour-changing) susceptibles
## generic emerging infection; awareness/behaviour.
## Awareness is triggered only by prevalence, so the disease-free state has Sa = 0 and the epidemic threshold R0 = beta/gamma is unchanged.
## Parameter values: 0 of 5 checked against the source; the rest are recalled or illustrative.
## Reference: The concept and qualitative result come from Funk S., Gilad E., Watkins C., Jansen V.A.A. (2009) "The spread of awareness and its impact on epidemic outbreaks", PNAS 106(16):6872–6877
## Source: https://doi.org/10.1073/pnas.0810762106
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, Sa, I, R
## Flows (time unit: days):
##   S -> Sa  : phi_a * I / N * S
##   Sa -> S  : omega_f * Sa
##   S -> I   : lambda * S
##   Sa -> I  : sigma_a * lambda * Sa
##   I -> R   : gamma * I

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 0.2)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.4,      # transmission rate (1/day) [illustrative]
  gamma = 0.2,     # recovery (1/day) [illustrative]
  phi_a = 2.0,     # awareness generation from contact with or news of cases (1/day) [illustrative]
  omega_f = 0.05,  # awareness fading (1/day) [illustrative]
  sigma_a = 0.3    # relative susceptibility of aware susceptibles [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9990,  # unaware susceptible
  Sa = 0,    # aware susceptible (reduced contacts)
  I = 10,    # infectious
  R = 0      # recovered
)

# ---- Model equations -------------------------------------------------
sir_awareness_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + Sa + I + R

    # Derived quantities
    lambda <- beta * I / N  # force of infection

    # Rates of change
    dS <- -phi_a * I / N * S + omega_f * Sa - lambda * S
    dSa <- phi_a * I / N * S - omega_f * Sa - sigma_a * lambda * Sa
    dI <- lambda * S + sigma_a * lambda * Sa - gamma * I
    dR <- gamma * I

    list(c(dS, dSa, dI, dR))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sir_awareness_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", Sa = "#333333", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SIR with aware (behaviour-changing) susceptibles")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
