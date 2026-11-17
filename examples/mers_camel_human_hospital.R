## MERS-CoV with camel spillover and hospital amplification
## Middle East respiratory syndrome coronavirus, dromedary camels → humans.
## The paper's central structural claim — repeated zoonotic index cases plus sub-critical human-to-human spread with hospital-based amplification — is captured by the two separate forces of infection and the two case accumulators, which let th
## Parameter values: 0 of 12 checked against the source; the rest are recalled or illustrative.
## Reference: Chowell G, Blumberg S, Simonsen L, Miller MA, Viboud C (2014). Synthesizing data and models for the spread of MERS-CoV, 2013: key role of index cases and hospital transmission. Epidemics 9:40–51
## Source: https://doi.org/10.1016/j.epidem.2014.09.011
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sc, Ic, Rc, S, E, Ii, Hh, R, D
## Flows (time unit: days):
##   (outside) -> Sc  : mu_c * Nc
##   Sc -> Ic         : lambda_c * Sc
##   Sc -> (outside)  : mu_c * Sc
##   S -> E           : (lambda_z + lambda_s) * S
##   E -> Ii          : sigma * E
##   Ii -> Hh         : eta * Ii
##   Ii -> R          : gamma_i * Ii
##   Hh -> D          : gamma_h * cfr * Hh
##   Hh -> R          : gamma_h * (1 - cfr) * Hh

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 400, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta_cc = 0.25,    # camel-to-camel transmission rate [illustrative] (per day)
  gamma_c = 1 / 14,  # camel recovery rate [illustrative] (per day)
  mu_c = 1 / 1825,   # camel birth = death rate, replenishes susceptible calves [illustrative] (per day)
  wane_c = 1 / 730,  # loss of camel immunity [illustrative] (per day)
  beta_ch = 2e-9,    # camel-to-human spillover coefficient [illustrative] (per camel per day)
  beta_com = 0.04,   # human-to-human transmission rate in the community [unverified] (per day)
  beta_hosp = 0.55,  # human-to-human transmission rate in hospital [unverified] (per day)
  sigma = 1 / 5.2,   # latency rate (mean incubation ≈5.2 d) [unverified] (per day)
  eta = 1 / 5,       # rate community case → hospitalisation (mean ≈5 d onset-to-admission) [unverified] (per day)
  gamma_i = 1 / 10,  # community removal rate [unverified] (per day)
  gamma_h = 1 / 13,  # hospital removal rate [unverified] (per day)
  cfr = 0.4          # case fatality among hospitalised cases [unverified] (dimensionless)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sc = 9000,  # susceptible camels
  Ic = 1000,  # infectious camels (actively shedding, mostly juveniles)
  Rc = 0,     # immune camels
  S = 1e6,    # susceptible humans
  E = 0,      # latent humans
  Ii = 0,     # infectious humans in the community
  Hh = 0,     # hospitalised humans (amplification setting)
  R = 0,      # recovered humans
  D = 0       # human deaths (cumulative)
)

# ---- Model equations -------------------------------------------------
mers_camel_human_hospital_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N_total <- S + E + Ii + Hh + R

    # Derived quantities
    Nc <- Sc + Ic + Rc                                                                    # camel population
    N <- S + E + Ii + Hh + R                                                              # living human population
    lambda_c <- beta_cc * Ic / Nc                                                         # camel force of infection
    lambda_z <- beta_ch * Ic                                                              # zoonotic hazard to a susceptible human
    lambda_s <- beta_com * Ii / N + beta_hosp * Hh / N                                    # human-to-human force of infection
    Rhuman <- beta_com / (eta + gamma_i) + (eta / (eta + gamma_i)) * beta_hosp / gamma_h  # human-to-human reproduction number (authors' headline result: overall R below or near 1, ≈0.4–1.5 across published estimates, with hospital
    CZ_rate <- (lambda_z) * S                                                             # rate of new cumulative zoonotic (index) cases (integrate with cumsum for the cumulative count)
    CS_rate <- (lambda_s) * S                                                             # rate of new cumulative secondary (human-acquired) cases (integrate with cumsum for the cumulative count)

    # Rates of change
    dSc <- mu_c * Nc - lambda_c * Sc - mu_c * Sc
    dIc <- lambda_c * Sc
    dRc <- 0
    dS <- -(lambda_z + lambda_s) * S
    dE <- (lambda_z + lambda_s) * S - sigma * E
    dIi <- sigma * E - eta * Ii - gamma_i * Ii
    dHh <- eta * Ii - gamma_h * cfr * Hh - gamma_h * (1 - cfr) * Hh
    dR <- gamma_i * Ii + gamma_h * (1 - cfr) * Hh
    dD <- gamma_h * cfr * Hh

    list(c(dSc, dIc, dRc, dS, dE, dIi, dHh, dR, dD), c(CZ_rate = CZ_rate, CS_rate = CS_rate))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = mers_camel_human_hospital_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sc = "#333333", Ic = "#B71C1C", Rc = "#4E8A3E", S = "#333333", E = "#E8C468", Ii = "#B71C1C", Hh = "#6A2C9C", R = "#4E8A3E", D = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "MERS-CoV with camel spillover and hospital amplification")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("CZ_rate", "CS_rate")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "MERS-CoV with camel spillover and hospital amplification: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
