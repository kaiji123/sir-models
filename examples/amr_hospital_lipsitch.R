## hospital colonisation with sensitive and resistant bacteria and patient turnover
## nosocomial antibiotic-resistant bacteria (e.g. VRE, MRSA).
## Total admissions equal total discharges, so N is constant.
## Parameter values: 0 of 8 checked against the source; the rest are recalled or illustrative.
## Reference: Lipsitch M., Bergstrom C.T., Levin B.R. (2000) "The epidemiology of antibiotic resistance in hospitals: paradoxes and prescriptions", PNAS 97(4):1938–1943
## Source: https://doi.org/10.1073/pnas.97.4.1938
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: X, Sc, Rc
## Flows (time unit: days):
##   (outside) -> X   : (1 - m_s - m_r) * mu_d * N
##   (outside) -> Sc  : m_s * mu_d * N
##   (outside) -> Rc  : m_r * mu_d * N
##   X -> Sc          : beta_s * Sc / N * X
##   X -> Rc          : beta_r * Rc / N * X
##   Sc -> X          : (tau1 + tau2 + gamma) * Sc
##   Rc -> X          : (tau2 + gamma) * Rc
##   X -> (outside)   : mu_d * X
##   Sc -> (outside)  : mu_d * Sc
##   Rc -> (outside)  : mu_d * Rc

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 365, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  mu_d = 0.1,    # discharge rate = admission rate (10-day stay) (1/day) [illustrative]
  m_s = 0.10,    # fraction admitted colonised with sensitive [illustrative]
  m_r = 0.01,    # fraction admitted colonised with resistant [illustrative]
  beta_s = 1.0,  # transmission, sensitive (1/day) [illustrative]
  beta_r = 1.0,  # transmission, resistant (1/day) [illustrative]
  tau1 = 0.1,    # treatment rate with drug 1 (clears sensitive only) (1/day) [illustrative]
  tau2 = 0.05,   # treatment rate with drug 2 (clears both) (1/day) [illustrative]
  gamma = 0.03   # spontaneous clearance (1/day) [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  X = 90,  # uncolonised patients
  Sc = 9,  # colonised with sensitive strain
  Rc = 1   # colonised with resistant strain
)

# ---- Model equations -------------------------------------------------
amr_hospital_lipsitch_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- X + Sc + Rc

    # Derived quantities
    prop_res <- Rc / (Sc + Rc)

    # Rates of change
    dX <- (1 - m_s - m_r) * mu_d * N - beta_s * Sc / N * X - beta_r * Rc / N * X + (tau1 + tau2 + gamma) * Sc + (tau2 + gamma) * Rc - mu_d * X
    dSc <- m_s * mu_d * N + beta_s * Sc / N * X - (tau1 + tau2 + gamma) * Sc - mu_d * Sc
    dRc <- m_r * mu_d * N + beta_r * Rc / N * X - (tau2 + gamma) * Rc - mu_d * Rc

    list(c(dX, dSc, dRc))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = amr_hospital_lipsitch_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(X = "#333333", Sc = "#333333", Rc = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "hospital colonisation with sensitive and resistant bacteria and patient turnover")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
