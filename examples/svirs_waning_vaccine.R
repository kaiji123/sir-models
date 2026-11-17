## SVIRS with imperfect (leaky) and waning vaccine and waning natural immunity
## e.g. influenza, pertussis-like.
## N is constant. This differs from the existing "SVIR leaky vaccine" model in having vaccine waning, newborn vaccination and waning natural immunity.
## Parameter values: 0 of 8 checked against the source; the rest are recalled or illustrative.
## Reference: The structure follows Arino J., McCluskey C.C., van den Driessche P. (2003) "Global results for an epidemic model with vaccination that exhibits backward bifurcation", SIAM J. Appl. Math. 64:260–276
## Source: https://doi.org/10.1137/S0036139902413829
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, V, I, R
## Flows (time unit: years):
##   (outside) -> S  : (1 - p_nb) * mu * N
##   (outside) -> V  : p_nb * mu * N
##   S -> V          : nu_v * S
##   V -> S          : theta_v * V
##   S -> I          : lambda * S
##   V -> I          : sigma_v * lambda * V
##   I -> R          : gamma * I
##   R -> S          : omega_r * R
##   S -> (outside)  : mu * S
##   V -> (outside)  : mu * V
##   I -> (outside)  : mu * I
##   R -> (outside)  : mu * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 100, by = 0.1)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  mu = 1 / 70,      # birth = death rate (1/yr) [illustrative]
  p_nb = 0.5,       # fraction of newborns vaccinated [illustrative]
  nu_v = 0.1,       # vaccination rate of susceptibles (1/yr) [illustrative]
  theta_v = 0.2,    # vaccine waning (mean 5 y) (1/yr) [illustrative]
  sigma_v = 0.2,    # relative susceptibility of vaccinees (vaccine efficacy 80%) [illustrative]
  beta = 100,       # transmission rate (1/yr) [illustrative]
  gamma = 365 / 7,  # recovery (7-day infectious period) (1/yr) [illustrative]
  omega_r = 0.1     # natural immunity waning (1/yr) [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9000,  # susceptible
  V = 0,     # vaccinated (partially protected)
  I = 10,    # infectious
  R = 990    # recovered
)

# ---- Model equations -------------------------------------------------
svirs_waning_vaccine_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + V + I + R

    # Derived quantities
    lambda <- beta * I / N  # force of infection

    # Rates of change
    dS <- (1 - p_nb) * mu * N - nu_v * S + theta_v * V - lambda * S + omega_r * R - mu * S
    dV <- p_nb * mu * N + nu_v * S - theta_v * V - sigma_v * lambda * V - mu * V
    dI <- lambda * S + sigma_v * lambda * V - gamma * I - mu * I
    dR <- gamma * I - omega_r * R - mu * R

    list(c(dS, dV, dI, dR))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = svirs_waning_vaccine_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", V = "#2457C5", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "SVIRS with imperfect (leaky) and waning vaccine and waning natural immunity")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
