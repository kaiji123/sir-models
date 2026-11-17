## Dengue host SEIR – vector SEI
## Human SEIR coupled to mosquito SEI with an extrinsic incubation period.
## Reference: Manore et al. (2014), J. Theor. Biol. 356:174-191 (structure)
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sh, Eh, Ih, Rh, Sv, Ev, Iv
## Flows (time unit: days):
##   Sh -> Eh         : lam_h * Sh
##   Eh -> Ih         : sigmaH * Eh
##   Ih -> Rh         : gammaH * Ih
##   (outside) -> Sv  : muV * Nv
##   Sv -> Ev         : lam_v * Sv
##   Ev -> Iv         : sigmaV * Ev
##   Sv -> (outside)  : muV * Sv
##   Ev -> (outside)  : muV * Ev
##   Iv -> (outside)  : muV * Iv

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 365, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  a = 0.5,           # bites per mosquito per day
  bh = 0.4,          # mosquito-to-human transmission probability
  bv = 0.4,          # human-to-mosquito transmission probability
  sigmaH = 1 / 5.5,  # 1/(intrinsic incubation period)
  gammaH = 1 / 5,    # human recovery rate
  sigmaV = 1 / 10,   # 1/(extrinsic incubation period)
  muV = 1 / 14       # mosquito death rate
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sh = 99990,   # Susceptible humans
  Eh = 0,       # Exposed humans
  Ih = 10,      # Infectious humans
  Rh = 0,       # Recovered humans
  Sv = 200000,  # Susceptible mosquitoes
  Ev = 0,       # Exposed mosquitoes
  Iv = 0        # Infectious mosquitoes
)

# ---- Model equations -------------------------------------------------
dengue_vector_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Derived quantities
    Nh <- Sh + Eh + Ih + Rh    # human population
    Nv <- Sv + Ev + Iv         # mosquito population
    lam_h <- a * bh * Iv / Nh  # force of infection on humans
    lam_v <- a * bv * Ih / Nh  # force of infection on mosquitoes

    # Flows recorded as outputs
    human_cases <- sigmaH * Eh

    # Rates of change
    dSh <- -lam_h * Sh
    dEh <- lam_h * Sh - sigmaH * Eh
    dIh <- sigmaH * Eh - gammaH * Ih
    dRh <- gammaH * Ih
    dSv <- muV * Nv - lam_v * Sv - muV * Sv
    dEv <- lam_v * Sv - sigmaV * Ev - muV * Ev
    dIv <- sigmaV * Ev - muV * Iv

    list(c(dSh, dEh, dIh, dRh, dSv, dEv, dIv), c(human_cases = human_cases))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = dengue_vector_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sh = "#333333", Eh = "#E8C468", Ih = "#B71C1C", Rh = "#4E8A3E", Sv = "#333333", Ev = "#E8C468", Iv = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Dengue host SEIR – vector SEI")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("human_cases")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "Dengue host SEIR – vector SEI: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
