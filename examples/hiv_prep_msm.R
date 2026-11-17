## HIV among MSM with PrEP and treatment as prevention
## HIV in MSM; PrEP.
## Parameter values: 0 of 13 checked against the source; the rest are recalled or illustrative.
## Reference: Structure is a generic composite in the style of many PrEP models. The efficacy value comes from McCormack S, Dunn DT, et al. Pre-exposure prophylaxis to prevent the acquisition of HIV-1 infection (PROUD). Lancet 2016;387:53–60
## Source: https://doi.org/10.1016/S0140-6736(15
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, Sp, Ia, Ic, Dx, Av
## Flows (time unit: years):
##   (outside) -> S   : (mu * N + nu_d * (Ic + Dx))
##   S -> Sp          : kappa_p * S
##   Sp -> S          : omega_p * Sp
##   S -> Ia          : foi * S
##   Sp -> Ia         : (1 - eff_prep) * foi * Sp
##   Ia -> Ic         : sig_a * Ia
##   Ic -> Dx         : tau_d * Ic
##   Dx -> Av         : tau_a * Dx
##   Ic -> (outside)  : nu_d * Ic
##   Dx -> (outside)  : nu_d * Dx
##   S -> (outside)   : mu * S
##   Sp -> (outside)  : mu * Sp
##   Ia -> (outside)  : mu * Ia
##   Ic -> (outside)  : mu * Ic
##   Dx -> (outside)  : mu * Dx
##   Av -> (outside)  : mu * Av

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 30, by = 0.05)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  cp = 5,           # partner change rate (partners/yr) [illustrative]
  beta_c = 0.1,     # transmission probability per partnership from chronic infection (per partnership) [illustrative]
  ra = 10,          # relative infectiousness of acute infection [illustrative]
  sig_a = 4,        # exit from acute stage (about 3 months) (/yr) [illustrative]
  eff_prep = 0.86,  # PrEP efficacy [recalled]
  eps_art = 0,      # infectiousness when virally suppressed (U=U) [illustrative]
  kappa_p = 0.2,    # PrEP uptake (/yr) [illustrative]
  omega_p = 0.5,    # PrEP discontinuation (/yr) [illustrative]
  tau_d = 0.5,      # diagnosis rate (/yr) [illustrative]
  tau_a = 2,        # ART initiation after diagnosis (/yr) [illustrative]
  nu_d = 0.1,       # HIV death if untreated (/yr) [illustrative]
  mu = 0.03,        # exit (/yr) [illustrative]
  tau_p = 4         # HIV testing rate while on PrEP (3-monthly); only relevant if a PrEP-specific infection compartment is added (/yr) [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 0.88,    # susceptible not on PrEP
  Sp = 0,      # susceptible on PrEP
  Ia = 0.005,  # acute (primary) infection
  Ic = 0.06,   # chronic undiagnosed
  Dx = 0.02,   # diagnosed, not on ART
  Av = 0.035   # on ART, virally suppressed
)

# ---- Model equations -------------------------------------------------
hiv_prep_msm_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + Sp + Ia + Ic + Dx + Av

    # Derived quantities
    foi <- cp * beta_c * (ra * Ia + Ic + Dx + eps_art * Av) / N  # force of infection

    # Rates of change
    dS <- (mu * N + nu_d * (Ic + Dx)) - kappa_p * S + omega_p * Sp - foi * S - mu * S
    dSp <- kappa_p * S - omega_p * Sp - (1 - eff_prep) * foi * Sp - mu * Sp
    dIa <- foi * S + (1 - eff_prep) * foi * Sp - sig_a * Ia - mu * Ia
    dIc <- sig_a * Ia - tau_d * Ic - nu_d * Ic - mu * Ic
    dDx <- tau_d * Ic - tau_a * Dx - nu_d * Dx - mu * Dx
    dAv <- tau_a * Dx - mu * Av

    list(c(dS, dSp, dIa, dIc, dDx, dAv))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = hiv_prep_msm_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", Sp = "#333333", Ia = "#B71C1C", Ic = "#B71C1C", Dx = "#111111", Av = "#E57373")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "HIV among MSM with PrEP and treatment as prevention")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
