## Urban yellow fever: SEIR humans with severe (toxic) phase and reactive vaccination; Aedes aegypti SEI [ILLUSTRATIVE-COMPOSITE]
## Yellow fever (urban cycle), e.g. Angola/DRC 2016-type outbreak.
## Parameter values: 0 of 13 checked against the source; the rest are recalled or illustrative.
## Reference: a standard structure. Natural-history values are recalled, not re-checked this session:
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sh, Eh, Ih, Sev, Rh, Vh, Sv, Ev, Iv
## Flows (time unit: days):
##   Sh -> Eh          : lambda_h * Sh
##   Sh -> Vh          : vacc * Sh
##   Eh -> Ih          : nu_h * Eh
##   Ih -> Sev         : p_sev * gam_h * Ih
##   Ih -> Rh          : (1 - p_sev) * gam_h * Ih
##   Sev -> (outside)  : cfr_sev * w_sev * Sev
##   Sev -> Rh         : (1 - cfr_sev) * w_sev * Sev
##   (outside) -> Sv   : mu_v * Nv
##   Sv -> Ev          : lambda_v * Sv
##   Ev -> Iv          : nu_v * Ev
##   Sv -> (outside)   : mu_v * Sv
##   Ev -> (outside)   : mu_v * Ev
##   Iv -> (outside)   : mu_v * Iv

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 365, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  a_bite = 0.5,    # (1/day)
  p_vh = 0.5,      # mosquito → human
  p_hv = 0.5,      # human → mosquito
  nu_h = 0.2326,   # 1/4.3 d intrinsic incubation (Johansson 2010) (1/day)
  gam_h = 0.2222,  # 1/4.5 d viraemia (1/day)
  p_sev = 0.12,    # fraction entering the severe phase (Johansson 2014)
  cfr_sev = 0.47,  # case fatality among severe (Johansson 2014)
  w_sev = 0.1429,  # 1/7 d duration of the severe phase (1/day)
  nu_v = 0.1,      # 1/10 d EIP (Johansson 2010) (1/day)
  mu_v = 0.0714,   # (1/day)
  v_rate = 0.01,   # campaign vaccination rate (1/day)
  v_eff = 0.99,    # vaccine efficacy (WHO)
  t_vac = 60       # campaign start (day)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sh = 79990,   # susceptible
  Eh = 0,       # incubating
  Ih = 10,      # viraemic (infectious)
  Sev = 0,      # severe/toxic phase (assumed non-infectious)
  Rh = 0,       # recovered
  Vh = 20000,   # vaccinated/immune
  Sv = 200000,  # susceptible mosquitoes
  Ev = 0,       # exposed mosquitoes
  Iv = 0        # infectious mosquitoes
)

# ---- Model equations -------------------------------------------------
yellow_fever_urban_vacc_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sh + Eh + Ih + Sev + Rh + Vh

    # Derived quantities
    Nh <- Sh + Eh + Ih + Sev + Rh + Vh
    Nv <- Sv + Ev + Iv
    lambda_h <- a_bite * p_vh * Iv / Nh            # human force of infection
    lambda_v <- a_bite * p_hv * Ih / Nh            # mosquito force of infection
    vacc <- ifelse(t >= t_vac, v_rate * v_eff, 0)  # effective vaccination rate

    # Rates of change
    dSh <- -lambda_h * Sh - vacc * Sh
    dEh <- lambda_h * Sh - nu_h * Eh
    dIh <- nu_h * Eh - p_sev * gam_h * Ih - (1 - p_sev) * gam_h * Ih
    dSev <- p_sev * gam_h * Ih - cfr_sev * w_sev * Sev - (1 - cfr_sev) * w_sev * Sev
    dRh <- (1 - p_sev) * gam_h * Ih + (1 - cfr_sev) * w_sev * Sev
    dVh <- vacc * Sh
    dSv <- mu_v * Nv - lambda_v * Sv - mu_v * Sv
    dEv <- lambda_v * Sv - nu_v * Ev - mu_v * Ev
    dIv <- nu_v * Ev - mu_v * Iv

    list(c(dSh, dEh, dIh, dSev, dRh, dVh, dSv, dEv, dIv))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = yellow_fever_urban_vacc_model, parms = parameters, hmax = 0.5)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sh = "#333333", Eh = "#E8C468", Ih = "#B71C1C", Sev = "#333333", Rh = "#4E8A3E", Vh = "#2457C5", Sv = "#333333", Ev = "#E8C468", Iv = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Urban yellow fever: SEIR humans with severe (toxic) phase and reactive vaccination; Aedes aegypti SEI [ILLUSTRATIVE-COMPOSITE]")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
