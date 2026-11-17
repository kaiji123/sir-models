## Rabies SEI (density-dependent)
## Wildlife rabies with density-dependent transmission (lam = beta * I, not divided by N)
## and logistic population regulation through density-dependent mortality.
## Reference: Anderson et al. (1981), Nature 289:765-771
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, I
## Flows (time unit: years):
##   (outside) -> S  : b * N
##   S -> E          : lam * S
##   E -> I          : sigma * E
##   S -> (outside)  : d * S
##   E -> (outside)  : d * E
##   I -> (outside)  : d * I
##   I -> (outside)  : alpha * I

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 30, by = 0.01)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  b = 1,        # per-capita birth rate per year
  mu = 0.5,     # baseline death rate per year
  K = 1000,     # carrying capacity
  beta = 0.08,  # density-dependent transmission coefficient per year
  sigma = 13,   # 1/(incubation period) per year (~4 weeks)
  alpha = 73    # rabies death rate per year (~5 days)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 900,  # Susceptible foxes
  E = 0,    # Incubating
  I = 5     # Rabid
)

# ---- Model equations -------------------------------------------------
rabies_density_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + I

    # Derived quantities
    d <- mu + (b - mu) * N / K  # density-dependent death rate
    lam <- beta * I             # force of infection (density dependent)

    # Flows recorded as outputs
    cases <- sigma * E
    rabies_deaths <- alpha * I

    # Rates of change
    dS <- b * N - lam * S - d * S
    dE <- lam * S - sigma * E - d * E
    dI <- sigma * E - d * I - alpha * I

    list(c(dS, dE, dI), c(cases = cases, rabies_deaths = rabies_deaths))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = rabies_density_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", I = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Rabies SEI (density-dependent)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("cases", "rabies_deaths")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (years)", ylab = "Value",
        main = "Rabies SEI (density-dependent): outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
