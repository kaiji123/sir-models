## SIR with reinfection of recovered (partial immunity; reinfection threshold)
## TB, pertussis, RSV, malaria-like partial immunity.
## Gomes et al. report that above the reinfection threshold reinfection dominates, prevalence is high and vaccination fails to protect (abstract via repository).
## Parameter values: 0 of 5 checked against the source; the rest are recalled or illustrative.
## Reference: Gomes M.G.M., White L.J., Medley G.F. (2004) "Infection, reinfection, and vaccination under suboptimal immune protection: epidemiological perspectives", J. Theor. Biol. 228(4):539–549
## Source: https://doi.org/10.1016/j.jtbi.2004.02.015
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, R, V
## Flows (time unit: years):
##   (outside) -> S  : (1 - p_nb) * mu * N
##   (outside) -> V  : p_nb * mu * N
##   S -> I          : lambda * S
##   R -> I          : sigma_r * lambda * R
##   V -> I          : sigma_r * lambda * V
##   I -> R          : gamma * I
##   S -> (outside)  : mu * S
##   I -> (outside)  : mu * I
##   R -> (outside)  : mu * R
##   V -> (outside)  : mu * V

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 0.2)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  mu = 1 / 70,     # birth/death (1/yr) [illustrative]
  gamma = 52,      # recovery (1 week) (1/yr) [illustrative]
  beta = 260,      # transmission (R0 ≈ 5) (1/yr) [illustrative]
  sigma_r = 0.25,  # relative susceptibility to reinfection (reinfection threshold at R0 = 1/sigma_r = 4) [illustrative]
  p_nb = 0         # newborn vaccine coverage (set to 0.9 to show vaccine failure above the reinfection threshold) [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9990,  # susceptible
  I = 10,    # infectious
  R = 0,     # recovered, partially immune
  V = 0      # vaccinated at birth (same partial protection as R)
)

# ---- Model equations -------------------------------------------------
sir_reinfection_partial_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + R + V

    # Derived quantities
    lambda <- beta * I / N

    # Rates of change
    dS <- (1 - p_nb) * mu * N - lambda * S - mu * S
    dI <- lambda * S + sigma_r * lambda * R + sigma_r * lambda * V - gamma * I - mu * I
    dR <- -sigma_r * lambda * R + gamma * I - mu * R
    dV <- p_nb * mu * N - sigma_r * lambda * V - mu * V

    list(c(dS, dI, dR, dV))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sir_reinfection_partial_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", R = "#4E8A3E", V = "#2457C5")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "SIR with reinfection of recovered (partial immunity; reinfection threshold)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
