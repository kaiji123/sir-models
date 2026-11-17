## Cryptosporidiosis with environmental oocysts and exponential dose–response
## Cryptosporidium (also adaptable to Giardia cysts).
## Parameter values: 0 of 9 checked against the source; the rest are recalled or illustrative.
## Reference: Author-constructed SIRS + environment structure [I]. The dose–response is Haas CN et al., exponential model with r ≈ 0.004 per oocyst (IOWA strain): J Am Water Works Assoc 1996;88(9):131–136 [M]. Giardia exponential r ≈ 0.0199 (Rose, Haas & Regli 1991, Am…
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, R, Oc
## Flows (time unit: days):
##   (outside) -> S   : mu_h * N
##   S -> I           : lam * S
##   I -> R           : gam_c * I
##   R -> S           : omega_c * R
##   (outside) -> Oc  : sh_vol * I
##   Oc -> (outside)  : mu_o * Oc

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 3, by = 0.005)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  exp_c = 1,              # water-consumption events per day [I] (day^-1)
  vol_d = 1,              # volume ingested per event [I] (L)
  r_dose = 0.004,         # exponential dose–response parameter [M] (per oocyst)
  gam_c = 1 / 14,         # recovery and end of shedding [M/I] (day^-1)
  sh_vol = 2e-4,          # shedding divided by water volume. [I] chosen to give R0 ≈ 2. (oocysts L^-1 day^-1 per infected)
  mu_o = 0.05,            # oocyst inactivation and flushing [I] (day^-1)
  omega_c = 1 / 365,      # waning immunity [I] (day^-1)
  mu_h = 1 / (70 * 365),  # birth/death [I] (day^-1)
  rem_t = 0               # fractional log-removal by treatment, for example 0.99 [I]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9990,  # susceptible
  I = 10,    # infected/shedding
  R = 0,     # immune (partial, waning)
  Oc = 0     # oocysts per litre in the community water source
)

# ---- Model equations -------------------------------------------------
crypto_oocyst_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + R

    # Derived quantities
    dose <- vol_d * Oc * (1 - rem_t)
    lam <- exp_c * (1 - exp(-r_dose * dose))                                # daily infection probability as a hazard
    R0_cr <- exp_c * r_dose * vol_d * sh_vol * N / (mu_o * (gam_c + mu_h))  # linearised R0 (about 2)

    # Rates of change
    dS <- mu_h * N - lam * S + omega_c * R
    dI <- lam * S - gam_c * I
    dR <- gam_c * I - omega_c * R
    dOc <- sh_vol * I - mu_o * Oc

    list(c(dS, dI, dR, dOc))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = crypto_oocyst_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", R = "#4E8A3E", Oc = "#333333")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Cryptosporidiosis with environmental oocysts and exponential dose–response")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
