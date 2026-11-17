## SIR with saturated (capacity-limited) treatment
## generic, limited treatment capacity.
## Mass-action (density-dependent) incidence as in the source papers.
## Parameter values: 0 of 7 checked against the source; the rest are recalled or illustrative.
## Reference: Wang W. & Ruan S. (2004) "Bifurcations in an epidemic model with constant removal rate of the infectives", J. Math. Anal. Appl. 291:775–793
## Source: https://doi.org/10.1016/j.jmaa.2003.11.043
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, R
## Flows (time unit: days):
##   (outside) -> S  : Lambda
##   S -> I          : beta * I * S
##   I -> R          : gamma * I
##   I -> R          : treat_pc * I
##   S -> (outside)  : mu * S
##   I -> (outside)  : (mu + d_dis) * I
##   R -> (outside)  : mu * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 3650, by = 5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  Lambda = 0.391,       # recruitment, equal to mu10000 (persons/day) [illustrative]
  mu = 1 / (70 * 365),  # natural death (1/day) [illustrative]
  beta = 3e-5,          # mass-action transmission coefficient [illustrative]
  gamma = 0.1,          # natural recovery (1/day) [illustrative]
  d_dis = 0.001,        # disease death (1/day) [illustrative]
  r_max = 0.1,          # maximal per-capita treatment rate at low prevalence (1/day) [illustrative]
  alpha_t = 0.01        # saturation constant (capacity effect) (1/person) [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9990,  # susceptible
  I = 10,    # infectious
  R = 0      # removed
)

# ---- Model equations -------------------------------------------------
sir_saturated_treatment_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + R

    # Derived quantities
    treat_pc <- r_max / (1 + alpha_t * I)  # per-capita treatment rate, so the total treatment is r_maxI/(1+alpha_tI)

    # Rates of change
    dS <- Lambda - beta * I * S - mu * S
    dI <- beta * I * S - gamma * I - treat_pc * I - (mu + d_dis) * I
    dR <- gamma * I + treat_pc * I - mu * R

    list(c(dS, dI, dR))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sir_saturated_treatment_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SIR with saturated (capacity-limited) treatment")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
