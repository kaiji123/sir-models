## Hepatitis E with latrines
## Hepatitis E outbreak with a latrine intervention (MMID practical, day 4).
## Reference: MMID practical: modelling interventions in R
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, I, R
## Flows (time unit: weeks):
##   (outside) -> S  : mui * P
##   S -> (outside)  : muo * S
##   S -> E          : lam * S
##   E -> (outside)  : muo * E
##   E -> I          : gamma * E
##   I -> (outside)  : muo * I
##   I -> R          : nui * I
##   R -> (outside)  : muo * R
##   R -> S          : omega * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(38, 73, by = 1 / 7)  # weeks

# ---- Parameters ------------------------------------------------------
parameters <- c(
  mui = 1 / (50 * 52),    # birth
  muo = 1 / (50 * 52),    # death
  R0 = 6.6,               # basic reproduction number
  omega = 1 / (10 * 52),  # rate of loss of immunity
  gamma = 1 / 2,          # latent -> infectious = 1/(latent period)
  nui = 1 / 4,            # recovery = 1/(duration of infection)
  report = 1 / 7,         # proportion of infections reported
  amp = 0,                # relative amplitude of seasonal forcing
  phi = 0,                # week of peak in seasonal forcing
  week_start = 38,        # first week of the simulation
  week_interv = 15,       # weeks after first case when intervention starts
  lat_cov_i = 0.24,       # intervention coverage of latrines
  lat_eff = 0.7           # efficacy of latrines in reducing transmission
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  E <- 1                 # Exposed
  I <- 0                 # Infectious
  R <- 0                 # Immune
  S <- 3700 - E - I - R  # Susceptible (non-immune)
  c(S = S, E = E, I = I, R = R)
})

# ---- Model equations -------------------------------------------------
hepE_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    P <- S + E + I + R

    # Derived quantities
    seas <- 1 + amp * cos(2 * pi * (t - phi) / 52)            # seasonal forcing
    beta <- R0 * (muo + nui) * (gamma + muo) / gamma          # transmission rate from R0
    lat_cov <- (t >= (week_start + week_interv)) * lat_cov_i  # latrine coverage (0 before intervention)
    latrine <- 1 - lat_eff * lat_cov                          # relative transmission with latrines
    lam <- latrine * beta * seas * I / P                      # force of infection
    inc <- report * gamma * E                                 # weekly reported incidence

    # Rates of change
    dS <- mui * P - muo * S - lam * S + omega * R
    dE <- lam * S - muo * E - gamma * E
    dI <- gamma * E - muo * I - nui * I
    dR <- nui * I - muo * R - omega * R

    list(c(dS, dE, dI, dR), c(inc = inc))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = hepE_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (weeks)", ylab = "Number of individuals",
        main = "Hepatitis E with latrines")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("inc")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (weeks)", ylab = "Value",
        main = "Hepatitis E with latrines: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
