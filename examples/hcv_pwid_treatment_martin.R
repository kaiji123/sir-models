## Hepatitis C among PWID with risk groups and DAA treatment as prevention
## HCV in people who inject drugs.
## Parameter values: 9 of 11 checked against the source; the rest are recalled or illustrative.
## Reference: Martin NK, Vickerman P, Grebely J, Hellard M, Hutchinson SJ, Lima VD, Foster GR, Dillon JF, Goldberg DJ, Dore GJ, Hickman M. Hepatitis C virus treatment for prevention among people who inject drugs: modeling treatment scale-up in the age of direct-acting…
## Source: https://doi.org/10.1002/hep.26431
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: X0, A0, Ch0, Tx0, Fail0, X1, A1, Ch1, Tx1, Fail1
## Flows (time unit: years):
##   (outside) -> X0     : (1 - f_high) * exit_r * Ntot
##   (outside) -> X1     : f_high * exit_r * Ntot
##   X0 -> A0            : foi0 * X0
##   X1 -> A1            : foi1 * X1
##   A0 -> X0            : delta_c * psi_a * A0
##   A0 -> Ch0           : (1 - delta_c) * psi_a * A0
##   A1 -> X1            : delta_c * psi_a * A1
##   A1 -> Ch1           : (1 - delta_c) * psi_a * A1
##   Ch0 -> Tx0          : tx_pc * Ch0
##   Ch1 -> Tx1          : tx_pc * Ch1
##   Tx0 -> X0           : alpha_svr * omega * Tx0
##   Tx0 -> Fail0        : (1 - alpha_svr) * omega * Tx0
##   Tx1 -> X1           : alpha_svr * omega * Tx1
##   Tx1 -> Fail1        : (1 - alpha_svr) * omega * Tx1
##   X1 -> X0            : kappa * X1
##   A1 -> A0            : kappa * A1
##   Ch1 -> Ch0          : kappa * Ch1
##   Tx1 -> Tx0          : kappa * Tx1
##   Fail1 -> Fail0      : kappa * Fail1
##   X0 -> (outside)     : exit_r * X0
##   A0 -> (outside)     : exit_r * A0
##   Ch0 -> (outside)    : exit_r * Ch0
##   Tx0 -> (outside)    : exit_r * Tx0
##   Fail0 -> (outside)  : exit_r * Fail0
##   X1 -> (outside)     : exit_r * X1
##   A1 -> (outside)     : exit_r * A1
##   Ch1 -> (outside)    : exit_r * Ch1
##   Tx1 -> (outside)    : exit_r * Tx1
##   Fail1 -> (outside)  : exit_r * Fail1

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 30, by = 0.05)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  pi_inf = 0.2,     # infection rate coefficient (/yr) [illustrative]
  Pi_h = 3.6,       # relative risk for high-risk PWID [verified]
  kappa = 12 / 14,  # high to low risk (duration 14 months) (/yr) [verified]
  f_high = 0.33,    # fraction of new injectors entering high risk [illustrative]
  delta_c = 0.26,   # proportion of acute infections clearing spontaneously [verified]
  psi_a = 2,        # exit from acute (6 months) (/yr) [verified]
  mu1 = 1 / 11,     # permanent cessation of injecting (11 yr career) (/yr) [verified]
  mu2 = 0.01,       # death rate (/yr) [verified]
  phi_tx = 8,       # treatments started per year (per 1000 PWID per yr) [verified]
  omega = 52 / 12,  # 1/treatment duration (12 weeks for IFN-free DAA) (/yr) [verified]
  alpha_svr = 0.9   # SVR with IFN-free DAAs [verified]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  X0 = 503,    # uninfected, low risk
  A0 = 5,      # acute HCV, low risk
  Ch0 = 150,   # chronic HCV, low risk
  Tx0 = 0,     # on treatment, low risk
  Fail0 = 12,  # failed treatment, low risk
  X1 = 230,    # uninfected, high risk
  A1 = 5,      # acute, high risk
  Ch1 = 90,    # chronic, high risk
  Tx1 = 0,     # on treatment, high risk
  Fail1 = 5    # failed treatment, high risk
)

# ---- Model equations -------------------------------------------------
hcv_pwid_treatment_martin_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- X0 + A0 + Ch0 + Tx0 + Fail0 + X1 + A1 + Ch1 + Tx1 + Fail1

    # Derived quantities
    Ntot <- X0 + A0 + Ch0 + Tx0 + Fail0 + X1 + A1 + Ch1 + Tx1 + Fail1                  # total PWID
    Chtot <- Ch0 + Ch1                                                                 # total chronic infections
    Infectious <- A0 + A1 + Ch0 + Ch1 + Fail0 + Fail1 + (1 - alpha_svr) * (Tx0 + Tx1)  # those on treatment who will achieve SVR are non-infectious [verified assumption]
    foi0 <- pi_inf * Infectious / Ntot                                                 # force of infection, low risk
    foi1 <- Pi_h * pi_inf * Infectious / Ntot                                          # force of infection, high risk
    tx_pc <- pmin(phi_tx / 1000 * Ntot / max(Chtot, 1e-6), 52)                         # per-capita treatment rate implementing a fixed number treated per year; capped
    exit_r <- mu1 + mu2                                                                # total exit rate (cessation plus death)

    # Rates of change
    dX0 <- (1 - f_high) * exit_r * Ntot - foi0 * X0 + delta_c * psi_a * A0 + alpha_svr * omega * Tx0 + kappa * X1 - exit_r * X0
    dA0 <- foi0 * X0 - delta_c * psi_a * A0 - (1 - delta_c) * psi_a * A0 + kappa * A1 - exit_r * A0
    dCh0 <- (1 - delta_c) * psi_a * A0 - tx_pc * Ch0 + kappa * Ch1 - exit_r * Ch0
    dTx0 <- tx_pc * Ch0 - alpha_svr * omega * Tx0 - (1 - alpha_svr) * omega * Tx0 + kappa * Tx1 - exit_r * Tx0
    dFail0 <- (1 - alpha_svr) * omega * Tx0 + kappa * Fail1 - exit_r * Fail0
    dX1 <- f_high * exit_r * Ntot - foi1 * X1 + delta_c * psi_a * A1 + alpha_svr * omega * Tx1 - kappa * X1 - exit_r * X1
    dA1 <- foi1 * X1 - delta_c * psi_a * A1 - (1 - delta_c) * psi_a * A1 - kappa * A1 - exit_r * A1
    dCh1 <- (1 - delta_c) * psi_a * A1 - tx_pc * Ch1 - kappa * Ch1 - exit_r * Ch1
    dTx1 <- tx_pc * Ch1 - alpha_svr * omega * Tx1 - (1 - alpha_svr) * omega * Tx1 - kappa * Tx1 - exit_r * Tx1
    dFail1 <- (1 - alpha_svr) * omega * Tx1 - kappa * Fail1 - exit_r * Fail1

    list(c(dX0, dA0, dCh0, dTx0, dFail0, dX1, dA1, dCh1, dTx1, dFail1))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = hcv_pwid_treatment_martin_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(X0 = "#333333", A0 = "#E57373", Ch0 = "#8E0E0E", Tx0 = "#333333", Fail0 = "#5D4037", X1 = "#333333", A1 = "#E57373", Ch1 = "#8E0E0E", Tx1 = "#333333", Fail1 = "#5D4037")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Hepatitis C among PWID with risk groups and DAA treatment as prevention")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
