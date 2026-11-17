## SIR with saturated incidence (Capasso–Serio)
## cholera (Bari 1973 epidemic), generic.
## The incidence g(I)S with g bounded represents protective behaviour or saturation at high prevalence.
## Parameter values: 0 of 3 checked against the source; the rest are recalled or illustrative.
## Reference: Capasso V. & Serio G. (1978) "A generalization of the Kermack–McKendrick deterministic epidemic model", Math. Biosci. 42:43–61
## Source: https://doi.org/10.1016/0025-5564(78
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, R
## Flows (time unit: days):
##   S -> I  : lambda * S
##   I -> R  : gamma * I

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 365, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 5e-5,      # transmission coefficient [illustrative]
  alpha_s = 0.002,  # saturation ("psychological"/crowding) constant (1/person) [illustrative]
  gamma = 0.2       # recovery (1/day) [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9990,  # susceptible
  I = 10,    # infectious
  R = 0      # recovered
)

# ---- Model equations -------------------------------------------------
sir_saturated_incidence_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + R

    # Derived quantities
    lambda <- beta * I / (1 + alpha_s * I)  # saturating force of infection

    # Rates of change
    dS <- -lambda * S
    dI <- lambda * S - gamma * I
    dR <- gamma * I

    list(c(dS, dI, dR))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sir_saturated_incidence_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SIR with saturated incidence (Capasso–Serio)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
