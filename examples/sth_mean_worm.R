## Soil-transmitted helminths: mean worm burden with environmental reservoir of infective stages (Anderson & May)
## Ascaris / Trichuris / hookworm.
## Parameter values: 0 of 8 checked against the source; the rest are recalled or illustrative.
## Reference: Anderson RM, May RM. Infectious Diseases of Humans: Dynamics and Control. Oxford University Press, 1991, chapters 16–17 (mean-worm-burden ODE with free-living reservoir, density-dependent fecundity and negative-binomial aggregation) [M]. Implemented for…
## Source: https://doi.org/10.1186/1756-3305-7-266
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Mw, Lr
## Flows (time unit: years):
##   (outside) -> Mw  : beta_l * Lr
##   Mw -> (outside)  : (mu_w + mu_h + tau_m) * Mw
##   (outside) -> Lr  : lam_e * Mw * fdd
##   Lr -> (outside)  : mu_l * Lr

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 30, by = 0.05)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  R0_sth = 3,   # basic reproduction number (Ascaris about 1–5) [M]
  mu_w = 1,     # adult worm death (Ascaris about 1 year) [M] (yr^-1)
  mu_h = 0.02,  # host death [I] (yr^-1)
  mu_l = 6,     # infective-stage loss (about 2 months) [I] (yr^-1)
  k_agg = 0.5,  # aggregation [M]
  z_dd = 0.93,  # per-worm density-dependent fecundity factor, exp(−γ) [M/I]
  lam_e = 1,    # egg output [I; normalisation] (stages per worm per yr, normalised)
  tau_m = 0     # continuous MDA kill rate. Pulsed alternative: ifelse((t - floor(t)) < 0.02, -log(1 - 0.750.95)/0.02, 0) (75% coverage, 95% efficacy of albendazole against Ascar
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Mw = 10,  # mean worm burden per host
  Lr = 5    # infective stages in the environment, normalised per host
)

# ---- Model equations -------------------------------------------------
sth_mean_worm_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Derived quantities
    beta_l <- R0_sth * (mu_w + mu_h) * mu_l / lam_e      # uptake rate of infective stages (18.36/yr)
    fdd <- (1 + Mw * (1 - z_dd) / k_agg)^(-(k_agg + 1))  # density-dependent fecundity under negative-binomial aggregation
    prev_sth <- 1 - (1 + Mw / k_agg)^(-k_agg)            # prevalence of infection

    # Rates of change
    dMw <- beta_l * Lr - (mu_w + mu_h + tau_m) * Mw
    dLr <- lam_e * Mw * fdd - mu_l * Lr

    list(c(dMw, dLr))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sth_mean_worm_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Mw = "#F3C6D3", Lr = "#E8C468")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Soil-transmitted helminths: mean worm burden with environmental reservoir of infective stages (Anderson & May)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
