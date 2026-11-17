## HPV (single high-risk type) SIRS-type model, two sexes, with vaccination
## HPV 16/18, vaccination.
## Parameter values: 0 of 13 checked against the source; the rest are recalled or illustrative.
## Reference: Simplified from Elbasha EH, Dasbach EJ, Insinga RP. Model for assessing human papillomavirus vaccination strategies. Emerg Infect Dis 2007;13(1):28–41
## Source: https://doi.org/10.3201/eid1301.060438
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sf, If, Rf, Vf, Sm, Im, Rm, Vm
## Flows (time unit: years):
##   (outside) -> Sf  : mu * Nf * (1 - vc_f)
##   (outside) -> Vf  : mu * Nf * vc_f
##   (outside) -> Sm  : mu * Nm * (1 - vc_m)
##   (outside) -> Vm  : mu * Nm * vc_m
##   Sf -> If         : foi_f * Sf
##   Vf -> If         : (1 - ve) * foi_f * Vf
##   If -> Rf         : p_imm_f * clr_f * If
##   If -> Sf         : (1 - p_imm_f) * clr_f * If
##   Rf -> Sf         : w_nat * Rf
##   Vf -> Sf         : w_vac * Vf
##   Sm -> Im         : foi_m * Sm
##   Vm -> Im         : (1 - ve) * foi_m * Vm
##   Im -> Rm         : p_imm_m * clr_m * Im
##   Im -> Sm         : (1 - p_imm_m) * clr_m * Im
##   Rm -> Sm         : w_nat * Rm
##   Vm -> Sm         : w_vac * Vm
##   Sf -> (outside)  : mu * Sf
##   If -> (outside)  : mu * If
##   Rf -> (outside)  : mu * Rf
##   Vf -> (outside)  : mu * Vf
##   Sm -> (outside)  : mu * Sm
##   Im -> (outside)  : mu * Im
##   Rm -> (outside)  : mu * Rm
##   Vm -> (outside)  : mu * Vm

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 100, by = 0.1)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  cp = 2,         # partner change rate (partners/yr) [illustrative]
  beta_mf = 0.8,  # male-to-female transmission (Elbasha used high values) (per partnership) [illustrative]
  beta_fm = 0.7,  # female-to-male transmission (per partnership) [illustrative]
  clr_f = 1,      # clearance in women (about 1 yr) (/yr) [illustrative]
  clr_m = 1.5,    # clearance in men (/yr) [illustrative]
  p_imm_f = 0.5,  # fraction of women gaining natural immunity on clearance [illustrative]
  p_imm_m = 0.2,  # fraction of men gaining natural immunity on clearance [illustrative]
  w_nat = 0.1,    # waning of natural immunity (/yr) [illustrative]
  w_vac = 0,      # waning of vaccine protection (assumed lifelong) (/yr) [illustrative]
  vc_f = 0.7,     # vaccine coverage of girls at entry [illustrative]
  vc_m = 0,       # vaccine coverage of boys at entry [illustrative]
  ve = 0.9,       # vaccine efficacy (degree protection) [illustrative]
  mu = 1 / 30     # exit from the sexually active population (/yr) [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sf = 0.9,   # susceptible women
  If = 0.05,  # infected women
  Rf = 0.05,  # naturally immune women
  Vf = 0,     # vaccinated women
  Sm = 0.9,   # susceptible men
  Im = 0.05,  # infected men
  Rm = 0.05,  # immune men
  Vm = 0      # vaccinated men
)

# ---- Model equations -------------------------------------------------
hpv_vaccination_simplified_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sf + If + Rf + Vf + Sm + Im + Rm + Vm

    # Derived quantities
    Nf <- Sf + If + Rf + Vf          # women total
    Nm <- Sm + Im + Rm + Vm          # men total
    foi_f <- cp * beta_mf * Im / Nm  # force of infection on women
    foi_m <- cp * beta_fm * If / Nf  # force of infection on men

    # Rates of change
    dSf <- mu * Nf * (1 - vc_f) - foi_f * Sf + (1 - p_imm_f) * clr_f * If + w_nat * Rf + w_vac * Vf - mu * Sf
    dIf <- foi_f * Sf + (1 - ve) * foi_f * Vf - p_imm_f * clr_f * If - (1 - p_imm_f) * clr_f * If - mu * If
    dRf <- p_imm_f * clr_f * If - w_nat * Rf - mu * Rf
    dVf <- mu * Nf * vc_f - (1 - ve) * foi_f * Vf - w_vac * Vf - mu * Vf
    dSm <- mu * Nm * (1 - vc_m) - foi_m * Sm + (1 - p_imm_m) * clr_m * Im + w_nat * Rm + w_vac * Vm - mu * Sm
    dIm <- foi_m * Sm + (1 - ve) * foi_m * Vm - p_imm_m * clr_m * Im - (1 - p_imm_m) * clr_m * Im - mu * Im
    dRm <- p_imm_m * clr_m * Im - w_nat * Rm - mu * Rm
    dVm <- mu * Nm * vc_m - (1 - ve) * foi_m * Vm - w_vac * Vm - mu * Vm

    list(c(dSf, dIf, dRf, dVf, dSm, dIm, dRm, dVm))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = hpv_vaccination_simplified_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sf = "#333333", If = "#B71C1C", Rf = "#4E8A3E", Vf = "#2457C5", Sm = "#333333", Im = "#B71C1C", Rm = "#4E8A3E", Vm = "#2457C5")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "HPV (single high-risk type) SIRS-type model, two sexes, with vaccination")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
