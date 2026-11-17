## Malaria SEIRS humans / SEI mosquitoes with saturating biting (Chitnis, Hyman & Cushing 2008) [VERIFIED-FULL]
## Plasmodium falciparum malaria, generic high- or low-transmission settings in sub-Saharan Africa.
## Parameter values: 0 of 17 checked against the source; the rest are recalled or illustrative.
## Reference: Chitnis N, Hyman JM, Cushing JM (2008). Determining important parameters in the spread of malaria through the sensitivity analysis of a mathematical model. Bull Math Biol 70(5):1272–1296. https://doi.org/10.1007/s11538-008-9299-0 — the equations (1a–g) and…
## Source: https://doi.org/10.1007/s11538-008-9299-0
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sh, Eh, Ih, Rh, Sv, Ev, Iv
## Flows (time unit: days):
##   (outside) -> Sh  : Lambda_h
##   (outside) -> Sh  : psi_h * Nh
##   Sh -> Eh         : lambda_h * Sh
##   Eh -> Ih         : nu_h * Eh
##   Ih -> Rh         : gam_h * Ih
##   Rh -> Sh         : rho_h * Rh
##   Ih -> (outside)  : delta_h * Ih
##   Sh -> (outside)  : fh * Sh
##   Eh -> (outside)  : fh * Eh
##   Ih -> (outside)  : fh * Ih
##   Rh -> (outside)  : fh * Rh
##   (outside) -> Sv  : psi_v * Nv
##   Sv -> Ev         : lambda_v * Sv
##   Ev -> Iv         : nu_v * Ev
##   Sv -> (outside)  : fv * Sv
##   Ev -> (outside)  : fv * Ev
##   Iv -> (outside)  : fv * Iv

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 7300, by = 10)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  Lambda_h = 0.033,  # immigration of humans (humans/day)
  psi_h = 1.1e-4,    # per-capita human birth rate (1/day)
  psi_v = 0.13,      # per-capita mosquito birth rate (1/day)
  sig_v = 0.50,      # bites one mosquito would want per day (1/day)
  sig_h = 19,        # maximum bites one human can receive per day (1/day)
  bet_hv = 0.022,    # P(infection) per bite, infectious mosquito → susceptible human
  bet_vh = 0.48,     # P(infection) per bite, infectious human → mosquito
  bet_vh_r = 0.048,  # P(infection) per bite, recovered human → mosquito
  nu_h = 0.10,       # human latent progression rate (1/day)
  nu_v = 0.091,      # mosquito latent progression rate (EIP ≈ 11–12 days) (1/day)
  gam_h = 0.0035,    # human recovery rate (1/day)
  delta_h = 9.0e-5,  # disease-induced death rate (1/day)
  rho_h = 5.5e-4,    # loss of immunity (1/day)
  mu1h = 1.6e-5,     # density-independent human death and emigration (1/day)
  mu2h = 3.0e-7,     # density-dependent human death
  mu1v = 0.033,      # density-independent mosquito death (1/day)
  mu2v = 2.0e-5      # density-dependent mosquito death
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sh = 500,   # susceptible humans
  Eh = 10,    # exposed (latent) humans
  Ih = 10,    # infectious humans
  Rh = 4,     # recovered: immune, asymptomatic, slightly infectious
  Sv = 4800,  # susceptible mosquitoes
  Ev = 30,    # exposed mosquitoes
  Iv = 20     # infectious mosquitoes
)

# ---- Model equations -------------------------------------------------
malaria_chitnis2008_seirs_sei_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sh + Eh + Ih + Rh

    # Derived quantities
    Nh <- Sh + Eh + Ih + Rh
    Nv <- Sv + Ev + Iv
    fh <- mu1h + mu2h * Nh                                              # human per-capita death rate
    fv <- mu1v + mu2v * Nv                                              # mosquito per-capita death rate
    bites <- sig_v * Nv * sig_h * Nh / (sig_v * Nv + sig_h * Nh)        # total human–mosquito contacts per day
    lambda_h <- (bites / Nh) * bet_hv * Iv / Nv                         # force of infection on humans
    lambda_v <- (bites / Nv) * (bet_vh * Ih / Nh + bet_vh_r * Rh / Nh)  # force of infection on mosquitoes

    # Rates of change
    dSh <- Lambda_h + psi_h * Nh - lambda_h * Sh + rho_h * Rh - fh * Sh
    dEh <- lambda_h * Sh - nu_h * Eh - fh * Eh
    dIh <- nu_h * Eh - gam_h * Ih - delta_h * Ih - fh * Ih
    dRh <- gam_h * Ih - rho_h * Rh - fh * Rh
    dSv <- psi_v * Nv - lambda_v * Sv - fv * Sv
    dEv <- lambda_v * Sv - nu_v * Ev - fv * Ev
    dIv <- nu_v * Ev - fv * Iv

    list(c(dSh, dEh, dIh, dRh, dSv, dEv, dIv))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = malaria_chitnis2008_seirs_sei_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sh = "#333333", Eh = "#E8C468", Ih = "#B71C1C", Rh = "#4E8A3E", Sv = "#333333", Ev = "#E8C468", Iv = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Malaria SEIRS humans / SEI mosquitoes with saturating biting (Chitnis, Hyman & Cushing 2008) [VERIFIED-FULL]")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
