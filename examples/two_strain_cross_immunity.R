## two-strain SIR with partial cross-immunity (reduced susceptibility)
## influenza strains, dengue serotypes (without ADE), RSV.
## There are 8 compartments and no simultaneous coinfection.
## Parameter values: 0 of 8 checked against the source; the rest are recalled or illustrative.
## Reference: Keeling M.J. & Rohani P. (2008), Chapter 4 ("Multi-pathogen / multi-host"), Section 4.1 and companion Program 4.1 (from memory; the companion site was unreachable). Original framework: Castillo-Chavez C., Hethcote H.W., Andreasen V., Levin S.A., Liu W.M.…
## Source: https://doi.org/10.1007/BF00275810
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I1, I2, R1, R2, J1, J2, R12
## Flows (time unit: days):
##   (outside) -> S    : mu * N
##   S -> I1           : lambda1 * S
##   S -> I2           : lambda2 * S
##   I1 -> R1          : gamma1 * I1
##   I2 -> R2          : gamma2 * I2
##   R1 -> J2          : sigma12 * lambda2 * R1
##   R2 -> J1          : sigma21 * lambda1 * R2
##   J1 -> R12         : gamma1 * J1
##   J2 -> R12         : gamma2 * J2
##   S -> (outside)    : mu * S
##   I1 -> (outside)   : mu * I1
##   I2 -> (outside)   : mu * I2
##   R1 -> (outside)   : mu * R1
##   R2 -> (outside)   : mu * R2
##   J1 -> (outside)   : mu * J1
##   J2 -> (outside)   : mu * J2
##   R12 -> (outside)  : mu * R12

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 50, by = 0.05)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  mu = 1 / (50 * 365),  # birth/death (1/day) [illustrative]
  beta1 = 0.5,          # strain 1 transmission (1/day) [illustrative]
  beta2 = 0.45,         # strain 2 transmission (1/day) [illustrative]
  gamma1 = 1 / 7,       # recovery, strain 1 (1/day) [illustrative]
  gamma2 = 1 / 7,       # recovery, strain 2 (1/day) [illustrative]
  sigma12 = 0.5,        # susceptibility to strain 2 after strain 1 (1 = no cross-immunity, 0 = full) [illustrative]
  sigma21 = 0.5,        # susceptibility to strain 1 after strain 2 [illustrative]
  N0 = 100000           # value taken from the notes text [recalled]
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  S <- 0.9 * N0    # naive to both
  I1 <- 10         # first infection, strain 1
  I2 <- 10         # first infection, strain 2
  R1 <- 0.05 * N0  # recovered from 1, susceptible to 2
  R2 <- 0.05 * N0  # recovered from 2, susceptible to 1
  J1 <- 0          # infected with strain 1 after recovering from 2
  J2 <- 0          # infected with strain 2 after recovering from 1
  R12 <- 0         # immune to both
  c(S = S, I1 = I1, I2 = I2, R1 = R1, R2 = R2, J1 = J1, J2 = J2, R12 = R12)
})

# ---- Model equations -------------------------------------------------
two_strain_cross_immunity_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I1 + I2 + R1 + R2 + J1 + J2 + R12

    # Derived quantities
    lambda1 <- beta1 * (I1 + J1) / N
    lambda2 <- beta2 * (I2 + J2) / N

    # Rates of change
    dS <- mu * N - lambda1 * S - lambda2 * S - mu * S
    dI1 <- lambda1 * S - gamma1 * I1 - mu * I1
    dI2 <- lambda2 * S - gamma2 * I2 - mu * I2
    dR1 <- gamma1 * I1 - sigma12 * lambda2 * R1 - mu * R1
    dR2 <- gamma2 * I2 - sigma21 * lambda1 * R2 - mu * R2
    dJ1 <- sigma21 * lambda1 * R2 - gamma1 * J1 - mu * J1
    dJ2 <- sigma12 * lambda2 * R1 - gamma2 * J2 - mu * J2
    dR12 <- gamma1 * J1 + gamma2 * J2 - mu * R12

    list(c(dS, dI1, dI2, dR1, dR2, dJ1, dJ2, dR12))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = two_strain_cross_immunity_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I1 = "#B71C1C", I2 = "#B71C1C", R1 = "#4E8A3E", R2 = "#4E8A3E", J1 = "#333333", J2 = "#333333", R12 = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "two-strain SIR with partial cross-immunity (reduced susceptibility)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
