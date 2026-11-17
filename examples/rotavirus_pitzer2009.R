## Rotavirus with maternal immunity and repeated infections of decreasing severity
## Rotavirus gastroenteritis.
## Parameter values: 0 of 15 checked against the source; the rest are recalled or illustrative.
## Reference: Pitzer VE, Viboud C, Simonsen L, Steiner C, Panozzo CA, Alonso WJ, Miller MA, Glass RI, Glasser JW, Parashar UD, Grenfell BT. Demographic variability, vaccination, and the spatiotemporal dynamics of rotavirus epidemics. Science 2009;325(5938):290–294
## Source: https://doi.org/10.1126/science.1172330
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: M, S0, I1, R1, S1, I2, R2, S2, I3, R3, S3
## Flows (time unit: weeks):
##   (outside) -> M  : mu_h * N
##   M -> S0         : omega_m * M
##   S0 -> I1        : lam * S0
##   I1 -> R1        : gam_1 * I1
##   R1 -> S1        : omega_r * R1
##   S1 -> I2        : sus_1 * lam * S1
##   I2 -> R2        : gam_2 * I2
##   R2 -> S2        : omega_r * R2
##   S2 -> I3        : sus_2 * lam * S2
##   I3 -> R3        : gam_2 * I3
##   R3 -> S3        : omega_r * R3
##   S3 -> I3        : sus_2 * lam * S3
##   M -> (outside)  : mu_h * M

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 30, by = 0.05)  # weeks

# ---- Parameters ------------------------------------------------------
parameters <- c(
  mu_h = 0.0003,     # birth/death. The paper used state-specific US birth rates (about 14–15 per 1000 per year) [I] (week^-1)
  omega_m = 1 / 13,  # loss of maternal immunity (13 weeks) [V via PMC9746532; M that Pitzer 2009 used it] (week^-1)
  gam_1 = 1,         # recovery from the first infection (1 week) [M] (week^-1)
  gam_2 = 2,         # recovery from later infections (half as long) [M] (week^-1)
  rho_inf = 0.5,     # relative infectiousness of 2nd+ infections [M]
  sus_1 = 0.62,      # relative susceptibility after 1 infection [V, Velázquez via PMC9746532]
  sus_2 = 0.37,      # relative susceptibility after 2+ infections [V, Velázquez via PMC9746532]. Pitzer may have used 0.35 [M].
  omega_r = 1 / 52,  # waning of temporary immunity (1 year) [M] (week^-1)
  R0_rv = 15,        # R0; fitted values were roughly 10–25 [M]
  amp_s = 0.15,      # seasonal amplitude [M/I]
  phi_s = 0,         # seasonal phase [I] (week)
  ps1 = 0.47,        # proportion symptomatic for the 1st, 2nd and 3rd+ infections. Velázquez 1996-type values; the 4th infection is about 0.2 [M]. Used for outputs only.
  ps2 = 0.25,        # proportion symptomatic for the 1st, 2nd and 3rd+ infections. Velázquez 1996-type values; the 4th infection is about 0.2 [M]. Used for outputs only.
  ps3 = 0.32,        # proportion symptomatic for the 1st, 2nd and 3rd+ infections. Velázquez 1996-type values; the 4th infection is about 0.2 [M]. Used for outputs only.
  N0 = 100000        # total initial population [illustrative default]
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  M <- 0.01 * N0   # maternal immunity
  S0 <- 0.05 * N0  # naive susceptible
  I1 <- 10         # first infection
  R1 <- 0          # temporary immunity after the 1st infection
  S1 <- 0.1 * N0   # susceptible after the 1st infection
  I2 <- 0          # second infection
  R2 <- 0          # temporary immunity after the 2nd infection
  S2 <- 0.2 * N0   # susceptible after the 2nd infection
  I3 <- 0          # third or later infection
  R3 <- 0          # temporary immunity after the 3rd or later infection
  S3 <- 0.64 * N0  # susceptible after 3+ infections
  c(M = M, S0 = S0, I1 = I1, R1 = R1, S1 = S1, I2 = I2, R2 = R2, S2 = S2, I3 = I3, R3 = R3, S3 = S3)
})

# ---- Model equations -------------------------------------------------
rotavirus_pitzer2009_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- M + S0 + I1 + R1 + S1 + I2 + R2 + S2 + I3 + R3 + S3

    # Derived quantities
    b0 <- R0_rv * (gam_1 + mu_h)
    bt <- b0 * (1 + amp_s * cos(6.2831853 * (t - phi_s) / 52))
    lam <- bt * (I1 + rho_inf * (I2 + I3)) / N
    inc_sev <- lam * (ps1 * S0 + ps2 * sus_1 * S1 + ps3 * sus_2 * S2)  # incidence of rotavirus gastroenteritis (RVGE)

    # Rates of change
    dM <- mu_h * N - omega_m * M - mu_h * M
    dS0 <- omega_m * M - lam * S0
    dI1 <- lam * S0 - gam_1 * I1
    dR1 <- gam_1 * I1 - omega_r * R1
    dS1 <- omega_r * R1 - sus_1 * lam * S1
    dI2 <- sus_1 * lam * S1 - gam_2 * I2
    dR2 <- gam_2 * I2 - omega_r * R2
    dS2 <- omega_r * R2 - sus_2 * lam * S2
    dI3 <- sus_2 * lam * S2 - gam_2 * I3 + sus_2 * lam * S3
    dR3 <- gam_2 * I3 - omega_r * R3
    dS3 <- omega_r * R3 - sus_2 * lam * S3

    list(c(dM, dS0, dI1, dR1, dS1, dI2, dR2, dS2, dI3, dR3, dS3))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = rotavirus_pitzer2009_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(M = "#F3C6D3", S0 = "#333333", I1 = "#B71C1C", R1 = "#4E8A3E", S1 = "#333333", I2 = "#B71C1C", R2 = "#4E8A3E", S2 = "#333333", I3 = "#B71C1C", R3 = "#4E8A3E", S3 = "#333333")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (weeks)", ylab = "Number of individuals",
        main = "Rotavirus with maternal immunity and repeated infections of decreasing severity")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
