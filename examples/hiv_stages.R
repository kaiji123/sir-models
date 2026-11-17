## HIV with disease stages
## Acute (I1), chronic (I2) and AIDS (A) stages with stage-specific infectiousness; AIDS deaths leave the population.
## Reference: Anderson & May (1991), Infectious Diseases of Humans, ch. 11
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I1, I2, A, D
## Flows (time unit: years):
##   (outside) -> S   : mu * N
##   S -> I1          : lam * S
##   I1 -> I2         : k1 * I1
##   I2 -> A          : k2 * I2
##   A -> D           : alpha * A
##   S -> (outside)   : mu * S
##   I1 -> (outside)  : mu * I1
##   I2 -> (outside)  : mu * I2
##   A -> (outside)   : mu * A

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 50, by = 0.05)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  cp = 2,       # partner change rate per year
  b1 = 0.1,     # transmission probability per partnership, acute
  b2 = 0.02,    # transmission probability per partnership, chronic
  b3 = 0.05,    # transmission probability per partnership, AIDS
  k1 = 4,       # rate of leaving acute stage per year (3 months)
  k2 = 1 / 9,   # rate of progression to AIDS per year
  alpha = 0.5,  # AIDS death rate per year
  mu = 1 / 35   # rate of entering/leaving the sexually active population
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 99900,  # Susceptible
  I1 = 100,   # Acute HIV
  I2 = 0,     # Chronic HIV
  A = 0,      # AIDS
  D = 0       # AIDS deaths
)

# ---- Model equations -------------------------------------------------
hiv_stages_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I1 + I2 + A

    # Derived quantities
    lam <- cp * (b1 * I1 + b2 * I2 + b3 * A) / N  # force of infection

    # Flows recorded as outputs
    incidence <- lam * S
    aids_deaths <- alpha * A

    # Rates of change
    dS <- mu * N - lam * S - mu * S
    dI1 <- lam * S - k1 * I1 - mu * I1
    dI2 <- k1 * I1 - k2 * I2 - mu * I2
    dA <- k2 * I2 - alpha * A - mu * A
    dD <- alpha * A

    list(c(dS, dI1, dI2, dA, dD), c(incidence = incidence, aids_deaths = aids_deaths))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = hiv_stages_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I1 = "#B71C1C", I2 = "#B71C1C", A = "#E57373", D = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "HIV with disease stages")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("incidence", "aids_deaths")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (years)", ylab = "Value",
        main = "HIV with disease stages: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
