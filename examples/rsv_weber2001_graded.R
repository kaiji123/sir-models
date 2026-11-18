## RSV with seasonal forcing and repeated infections with graded immunity (Weber, Weber & Milligan 2001)
## RSV.
## Parameter values: 0 of 11 checked against the source; the rest are recalled or illustrative.
## Reference: Weber A, Weber M, Milligan P. "Modeling epidemics caused by respiratory syncytial virus (RSV)." Math Biosci 2001;172(2):95–113. https://doi.org/10.1016/S0025-5564(01)00066-9
## Source: https://doi.org/10.1016/S0025-5564(01)00066-9
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: M, S1, E1, I1, R1, S2, E2, I2, R2
## Flows (time unit: years):
##   (outside) -> M  : mu * N
##   M -> S1         : xi_m * M
##   S1 -> E1        : foi1 * S1
##   E1 -> I1        : sig_lat * E1
##   I1 -> R1        : nu * I1
##   R1 -> S2        : gam_w * R1
##   S2 -> E2        : f2 * foi1 * S2
##   E2 -> I2        : sig_lat * E2
##   I2 -> R2        : nu * I2
##   R2 -> S2        : gam_w * R2
##   M -> (outside)  : mu * M

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 30, by = 0.05)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  b0 = 60,            # mean transmission rate (1/year) [illustrative]
  b1 = 0.3,           # seasonal amplitude [illustrative]
  phase = 0,          # seasonal phase (rad) [illustrative]
  two_pi = 6.283185,  # constant
  sig_lat = 91.25,    # 4-day latency (1/year)
  nu = 36.5,          # 10-day infectious period (1/year)
  gam_w = 1.825,      # loss of immunity after about 200 days (1/year)
  xi_m = 2,           # loss of maternal protection (about 6 months) (1/year) [illustrative]
  f2 = 0.5,           # relative susceptibility of S2 [illustrative]
  mu = 0.0133,        # birth = death rate, i.e. 1/75 y (1/year) [illustrative]
  N0 = 1e6            # value taken from the notes text [recalled]
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  M <- 0.01 * N0   # maternally protected
  S1 <- 0.05 * N0  # susceptible, never infected
  E1 <- 0          # exposed, first infection
  I1 <- 100        # infectious, first infection
  R1 <- 0.04 * N0  # temporarily immune after the first infection
  S2 <- 0.85 * N0  # susceptible, previously infected (reduced susceptibility)
  E2 <- 0          # exposed, reinfection
  I2 <- 100        # infectious, reinfection
  R2 <- 0.05 * N0  # temporarily immune after reinfection
  c(M = M, S1 = S1, E1 = E1, I1 = I1, R1 = R1, S2 = S2, E2 = E2, I2 = I2, R2 = R2)
})

# ---- Model equations -------------------------------------------------
rsv_weber2001_graded_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- M + S1 + E1 + I1 + R1 + S2 + E2 + I2 + R2

    # Derived quantities
    beta1 <- b0 * (1 + b1 * cos(two_pi * t + phase))
    foi1 <- beta1 * (I1 + I2) / N

    # Rates of change
    dM <- mu * N - xi_m * M - mu * M
    dS1 <- xi_m * M - foi1 * S1
    dE1 <- foi1 * S1 - sig_lat * E1
    dI1 <- sig_lat * E1 - nu * I1
    dR1 <- nu * I1 - gam_w * R1
    dS2 <- gam_w * R1 - f2 * foi1 * S2 + gam_w * R2
    dE2 <- f2 * foi1 * S2 - sig_lat * E2
    dI2 <- sig_lat * E2 - nu * I2
    dR2 <- nu * I2 - gam_w * R2

    list(c(dM, dS1, dE1, dI1, dR1, dS2, dE2, dI2, dR2))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = rsv_weber2001_graded_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(M = "#F3C6D3", S1 = "#333333", E1 = "#E8C468", I1 = "#B71C1C", R1 = "#4E8A3E", S2 = "#333333", E2 = "#E8C468", I2 = "#B71C1C", R2 = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "RSV with seasonal forcing and repeated infections with graded immunity (Weber, Weber & Milligan 2001)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
