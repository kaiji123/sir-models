## Leptospirosis: rodent reservoir, environmental leptospires, humans as dead-end hosts
## Leptospirosis.
## Parameter values: 0 of 13 checked against the source; the rest are recalled or illustrative.
## Reference: Structure after Holt J, Davis S, Leirs H. A model of leptospirosis infection in an African rodent to determine risk to humans: seasonal fluctuations and the impact of rodent control. Acta Trop 2006;99(2–3):218–225
## Source: https://doi.org/10.1016/j.actatropica.2006.08.003
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sh, Ih, Rh, Sr, Ir, Bl
## Flows (time unit: days):
##   (outside) -> Sr  : b_r * Nr * (1 - Nr / K_r)
##   Sr -> Ir         : (bd_r * Ir / Nr + be_r * envB) * Sr
##   Sr -> (outside)  : mu_r * Sr
##   Ir -> (outside)  : mu_r * Ir
##   (outside) -> Bl  : sh_r * Ir
##   Bl -> (outside)  : del_b * Bl
##   (outside) -> Sh  : mu_h * N
##   Sh -> Ih         : be_h * envB * Sh
##   Ih -> Rh         : gam_h * Ih
##   Ih -> (outside)  : cfr_h * Ih
##   Rh -> Sh         : omega_h * Rh

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 5, by = 0.005)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  b_r = 1 / 180,            # rodent birth rate per rodent [I] (day^-1)
  K_r = 1000,               # rodent carrying capacity [I] (rodents)
  mu_r = 1 / 365,           # rodent death [I] (day^-1)
  bd_r = 0.002,             # direct rodent transmission, scaled by Ir/Nr [I] (day^-1 per infected rodent)
  be_r = 0.05,              # rodent infection from the environment [I] (day^-1)
  sh_r = 1,                 # shedding [I] (units day^-1 per infected rodent)
  del_b = 0.1,              # leptospire decay [I] (day^-1)
  K_b = 1000,               # half-saturation [I] (units)
  be_h = 1e-4,              # human exposure rate. It scales with flooding: multiply by (1 + amp_fcos(6.2831853t/365)) [I] (day^-1)
  gam_h = 1 / 14,           # human recovery [M/I] (day^-1)
  cfr_h = 0.005,            # disease death hazard (CFR about 5–10% in severe cases) [I] (day^-1)
  omega_h = 1 / (5 * 365),  # waning serovar immunity [I] (day^-1)
  mu_h = 1 / (70 * 365)     # human birth/death [I] (day^-1)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sh = 9990,  # susceptible humans
  Ih = 0,     # infected (clinical) humans
  Rh = 0,     # recovered humans
  Sr = 900,   # susceptible rodents
  Ir = 100,   # infected (chronically shedding) rodents
  Bl = 0      # leptospires in the environment (arbitrary units)
)

# ---- Model equations -------------------------------------------------
lepto_rodent_env_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sh + Ih + Rh

    # Derived quantities
    Nr <- Sr + Ir
    envB <- Bl / (K_b + Bl)

    # Rates of change
    dSh <- mu_h * N - be_h * envB * Sh + omega_h * Rh
    dIh <- be_h * envB * Sh - gam_h * Ih - cfr_h * Ih
    dRh <- gam_h * Ih - omega_h * Rh
    dSr <- b_r * Nr * (1 - Nr / K_r) - (bd_r * Ir / Nr + be_r * envB) * Sr - mu_r * Sr
    dIr <- (bd_r * Ir / Nr + be_r * envB) * Sr - mu_r * Ir
    dBl <- sh_r * Ir - del_b * Bl

    list(c(dSh, dIh, dRh, dSr, dIr, dBl))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = lepto_rodent_env_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sh = "#333333", Ih = "#B71C1C", Rh = "#4E8A3E", Sr = "#333333", Ir = "#B71C1C", Bl = "#333333")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Leptospirosis: rodent reservoir, environmental leptospires, humans as dead-end hosts")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
