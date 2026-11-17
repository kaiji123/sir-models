## Tuberculosis (fast and slow progression)
## New infections progress rapidly to disease with probability p, otherwise enter long-term latency (L)
## with slow reactivation; treated/self-cured cases can relapse.
## Reference: Vynnycky & White (2010), An Introduction to Infectious Disease Modelling, ch. 3; Blower et al. (1995), Nature Medicine 1:815-821
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, L, I, R
## Flows (time unit: years):
##   (outside) -> S  : (mu * N + muTB * I)
##   S -> I          : p * lam * S
##   S -> L          : (1 - p) * lam * S
##   L -> I          : v * L
##   I -> R          : r * I
##   R -> I          : relapse * R
##   S -> (outside)  : mu * S
##   L -> (outside)  : mu * L
##   I -> (outside)  : mu * I
##   I -> (outside)  : muTB * I
##   R -> (outside)  : mu * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 0.1)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 10,        # effective contacts per case per year
  p = 0.1,          # proportion progressing rapidly to disease
  v = 0.001,        # reactivation rate per year
  r = 0.5,          # rate of cure / self-cure per year
  relapse = 0.002,  # relapse rate per year
  mu = 1 / 70,      # background death rate per year
  muTB = 0.2        # TB mortality rate per year
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 99990,  # Susceptible
  L = 0,      # Latent infection
  I = 10,     # Active TB
  R = 0       # Recovered
)

# ---- Model equations -------------------------------------------------
tb_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + L + I + R

    # Derived quantities
    lam <- beta * I / N  # force of infection

    # Flows recorded as outputs
    fast_cases <- p * lam * S

    # Rates of change
    dS <- (mu * N + muTB * I) - p * lam * S - (1 - p) * lam * S - mu * S
    dL <- (1 - p) * lam * S - v * L - mu * L
    dI <- p * lam * S + v * L - r * I + relapse * R - mu * I - muTB * I
    dR <- r * I - relapse * R - mu * R

    list(c(dS, dL, dI, dR), c(fast_cases = fast_cases))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = tb_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", L = "#E8C468", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Tuberculosis (fast and slow progression)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("fast_cases")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (years)", ylab = "Value",
        main = "Tuberculosis (fast and slow progression): outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
