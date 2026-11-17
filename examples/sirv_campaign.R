## SIR with a vaccination campaign
## Mass vaccination of susceptibles at rate v starting on day t_vac (assumes a perfect vaccine).
## Reference: Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals, §8.1
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, R, V
## Flows (time unit: days):
##   S -> I  : lam * S
##   I -> R  : gamma * I
##   S -> V  : vacc * S

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.3,      # transmission rate per day
  gamma = 1 / 10,  # recovery rate
  v = 0.02,        # vaccination rate per day once the campaign starts
  t_vac = 30       # day the campaign starts
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9990,  # Susceptible
  I = 10,    # Infectious
  R = 0,     # Recovered
  V = 0      # Vaccinated
)

# ---- Model equations -------------------------------------------------
sirv_campaign_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + R + V

    # Derived quantities
    lam <- beta * I / N       # force of infection
    vacc <- v * (t >= t_vac)  # vaccination rate at time t

    # Flows recorded as outputs
    incidence <- lam * S

    # Rates of change
    dS <- -lam * S - vacc * S
    dI <- lam * S - gamma * I
    dR <- gamma * I
    dV <- vacc * S

    list(c(dS, dI, dR, dV), c(incidence = incidence))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sirv_campaign_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", R = "#4E8A3E", V = "#2457C5")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "SIR with a vaccination campaign")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("incidence")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "SIR with a vaccination campaign: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
