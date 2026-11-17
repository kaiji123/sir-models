## Cholera SIWR (water reservoir)
## Direct transmission plus indirect transmission through a contaminated water compartment W
## (W is a normalised pathogen concentration, so it is excluded from N).
## Reference: Tien & Earn (2010), Bull. Math. Biol. 72:1506-1533
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, R, W
## Flows (time unit: days):
##   S -> I          : lam * S
##   I -> R          : gamma * I
##   (outside) -> W  : xi * I / N
##   W -> (outside)  : xi * W

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  betaI = 0.1,    # direct transmission rate per day
  betaW = 0.4,    # water-borne transmission rate per day
  gamma = 1 / 5,  # recovery rate
  xi = 1 / 30     # shedding / decay rate of pathogen in water
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9990,  # Susceptible
  I = 10,    # Infectious
  R = 0,     # Recovered
  W = 0      # Water (pathogen)
)

# ---- Model equations -------------------------------------------------
cholera_siwr_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + R

    # Derived quantities
    lam <- betaI * I / N + betaW * W  # force of infection

    # Flows recorded as outputs
    incidence <- lam * S

    # Rates of change
    dS <- -lam * S
    dI <- lam * S - gamma * I
    dR <- gamma * I
    dW <- xi * I / N - xi * W

    list(c(dS, dI, dR, dW), c(incidence = incidence))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = cholera_siwr_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", R = "#4E8A3E", W = "#4FC3F7")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Cholera SIWR (water reservoir)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("incidence")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "Cholera SIWR (water reservoir): outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
