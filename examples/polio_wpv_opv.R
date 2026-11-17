## Wild poliovirus with OPV vaccine virus and circulating vaccine-derived (cVDPV) lineage
## Poliomyelitis (WPV, Sabin OPV, cVDPV).
## Parameter values: 0 of 11 checked against the source; the rest are recalled or illustrative.
## Reference: Structure inspired by Fine PEM, Carneiro IAM. Transmissibility and persistence of oral polio vaccine viruses: implications for the global poliomyelitis eradication initiative. Am J Epidemiol 1999;150(10):1001–1021
## Source: https://doi.org/10.1093/oxfordjournals.aje.a009924
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, IW, IV, IVD, R
## Flows (time unit: days):
##   (outside) -> S   : mu_h * N * (1 - p_opv * take_v)
##   (outside) -> IV  : mu_h * N * p_opv * take_v
##   S -> IW          : bw * IW / N * S
##   S -> IV          : bv * IV / N * (1 - p_rev) * S
##   S -> IVD         : (bv * IV / N * p_rev + bvd * IVD / N) * S
##   IW -> R          : gam_p * IW
##   IV -> R          : gam_p * IV
##   IVD -> R         : gam_p * IVD
##   R -> S           : omega_p * R
##   S -> (outside)   : mu_h * S

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 20, by = 0.02)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  mu_h = 1 / (50 * 365),    # birth/death [I] (day^-1)
  gam_p = 1 / 28,           # recovery (about 4 weeks of excretion) [M] (day^-1)
  R0_w = 6,                 # WPV R0; reported range about 4–13 in high-transmission settings [M]
  rel_v = 0.4,              # Sabin virus transmissibility relative to WPV [I]
  rel_vd = 0.9,             # cVDPV transmissibility relative to WPV [I]
  p_rev = 0.01,             # fraction of Sabin onward transmissions producing revertant (cVDPV) virus [I]
  p_opv = 0.8,              # routine OPV coverage at birth (cohort) [I]
  take_v = 0.5,             # per-dose take, intestinal infection [I]
  omega_p = 1 / (5 * 365),  # waning of intestinal immunity [I] (day^-1)
  par_ratio = 0.005,        # paralysis-to-infection ratio, about 1/200 for WPV1 [M, widely cited]
  N0 = 100000               # total initial population [illustrative default]
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  S <- 0.1 * N0  # susceptible to intestinal infection
  IW <- 10       # wild-type infection
  IV <- 0        # OPV (Sabin-like) infection
  IVD <- 0       # cVDPV (reverted) infection
  R <- 0.9 * N0  # intestinally immune
  c(S = S, IW = IW, IV = IV, IVD = IVD, R = R)
})

# ---- Model equations -------------------------------------------------
polio_wpv_opv_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + IW + IV + IVD + R

    # Derived quantities
    bw <- R0_w * gam_p
    bv <- rel_v * R0_w * gam_p
    bvd <- rel_vd * R0_w * gam_p
    inc_par <- par_ratio * (bw * IW + bvd * IVD) / N * S  # paralytic incidence

    # Rates of change
    dS <- mu_h * N * (1 - p_opv * take_v) - bw * IW / N * S - bv * IV / N * (1 - p_rev) * S - (bv * IV / N * p_rev + bvd * IVD / N) * S + omega_p * R - mu_h * S
    dIW <- bw * IW / N * S - gam_p * IW
    dIV <- mu_h * N * p_opv * take_v + bv * IV / N * (1 - p_rev) * S - gam_p * IV
    dIVD <- (bv * IV / N * p_rev + bvd * IVD / N) * S - gam_p * IVD
    dR <- gam_p * IW + gam_p * IV + gam_p * IVD - omega_p * R

    list(c(dS, dIW, dIV, dIVD, dR))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = polio_wpv_opv_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", IW = "#B71C1C", IV = "#B71C1C", IVD = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Wild poliovirus with OPV vaccine virus and circulating vaccine-derived (cVDPV) lineage")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
