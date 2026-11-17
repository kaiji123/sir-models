## Ross–Macdonald (malaria)
## SIS humans and SI mosquitoes. Two separate populations, so totals are defined as variables.
## Reference: Smith et al. (2012), PLoS Pathogens 8:e1002588; Anderson & May (1991), Infectious Diseases of Humans, ch. 14
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sh, Ih, Sv, Iv
## Flows (time unit: days):
##   Sh -> Ih         : lam_h * Sh
##   Ih -> Sh         : r * Ih
##   (outside) -> Sv  : g * Nv
##   Sv -> Iv         : lam_v * Sv
##   Sv -> (outside)  : g * Sv
##   Iv -> (outside)  : g * Iv

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 730, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  a = 0.3,      # human bites per mosquito per day
  b = 0.5,      # probability an infectious bite infects a human
  cv = 0.5,     # probability a mosquito is infected biting an infected human
  r = 1 / 100,  # human recovery rate per day
  g = 1 / 10    # mosquito death rate per day
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sh = 9990,   # Susceptible humans
  Ih = 10,     # Infected humans
  Sv = 50000,  # Susceptible mosquitoes
  Iv = 0       # Infectious mosquitoes
)

# ---- Model equations -------------------------------------------------
ross_macdonald_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Derived quantities
    Nh <- Sh + Ih              # human population
    Nv <- Sv + Iv              # mosquito population
    lam_h <- a * b * Iv / Nh   # force of infection on humans
    lam_v <- a * cv * Ih / Nh  # force of infection on mosquitoes

    # Flows recorded as outputs
    human_incidence <- lam_h * Sh

    # Rates of change
    dSh <- -lam_h * Sh + r * Ih
    dIh <- lam_h * Sh - r * Ih
    dSv <- g * Nv - lam_v * Sv - g * Sv
    dIv <- lam_v * Sv - g * Iv

    list(c(dSh, dIh, dSv, dIv), c(human_incidence = human_incidence))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = ross_macdonald_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sh = "#333333", Ih = "#B71C1C", Sv = "#333333", Iv = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Ross–Macdonald (malaria)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("human_incidence")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "Ross–Macdonald (malaria): outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
