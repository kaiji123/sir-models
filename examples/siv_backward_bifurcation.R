## SIS with vaccination (Kribs-Zaleta & Velasco-Hernández)
## generic SIS infection with imperfect waning vaccine.
## This is the canonical small model where an endemic equilibrium coexists with a stable disease-free state when the vaccine-reduced R_phi < 1.
## Parameter values: 0 of 6 checked against the source; the rest are recalled or illustrative.
## Reference: Kribs-Zaleta C.M. & Velasco-Hernández J.X. (2000) "A simple vaccination model with multiple endemic states", Math. Biosci. 164:183–201
## Source: https://doi.org/10.1016/S0025-5564(00
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, V
## Flows (time unit: years):
##   (outside) -> S  : mu * N
##   S -> I          : lambda * S
##   I -> S          : gamma * I
##   S -> V          : phi_v * S
##   V -> S          : theta_v * V
##   V -> I          : sigma_v * lambda * V
##   S -> (outside)  : mu * S
##   I -> (outside)  : mu * I
##   V -> (outside)  : mu * V

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 50, by = 0.05)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  mu = 0.02,      # birth/death (turnover) (1/yr) [illustrative]
  beta = 20,      # transmission (1/yr) [illustrative]
  gamma = 10,     # recovery, returning to susceptible (1/yr) [illustrative]
  phi_v = 1.0,    # vaccination rate (1/yr) [illustrative]
  theta_v = 0.1,  # vaccine waning (1/yr) [illustrative]
  sigma_v = 0.3   # relative susceptibility of the vaccinated [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9900,  # susceptible
  I = 100,   # infectious
  V = 0      # vaccinated
)

# ---- Model equations -------------------------------------------------
siv_backward_bifurcation_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + V

    # Derived quantities
    lambda <- beta * I / N

    # Rates of change
    dS <- mu * N - lambda * S + gamma * I - phi_v * S + theta_v * V - mu * S
    dI <- lambda * S - gamma * I + sigma_v * lambda * V - mu * I
    dV <- phi_v * S - theta_v * V - sigma_v * lambda * V - mu * V

    list(c(dS, dI, dV))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = siv_backward_bifurcation_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", V = "#2457C5")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "SIS with vaccination (Kribs-Zaleta & Velasco-Hernández)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
