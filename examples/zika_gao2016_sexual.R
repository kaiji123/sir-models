## Zika with mosquito-borne and sexual transmission (Gao et al. 2016) [VERIFIED-PARAMS]
## Zika virus, Latin America 2015–16 (fitted to Brazil, Colombia, El Salvador).
## Parameter values: 0 of 15 checked against the source; the rest are recalled or illustrative.
## Reference: Gao D, Lou Y, He D, Porco TC, Kuang Y, Chowell G, Ruan S (2016). Prevention and control of Zika as a mosquito-borne and sexually transmitted disease: a mathematical modeling analysis. Sci Rep 6:28070. https://doi.org/10.1038/srep28070 ; PMC:…
## Source: https://doi.org/10.1038/srep28070
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sh, Eh, I1, I2, Ah, Rh, Sv, Ev, Iv
## Flows (time unit: days):
##   Sh -> Eh         : lambda_h * Sh
##   Eh -> I1         : theta * nu_h * Eh
##   Eh -> Ah         : (1 - theta) * nu_h * Eh
##   I1 -> I2         : gam_h1 * I1
##   I2 -> Rh         : gam_h2 * I2
##   Ah -> Rh         : gam_ha * Ah
##   (outside) -> Sv  : mu_v * Nv
##   Sv -> Ev         : lambda_v * Sv
##   Ev -> Iv         : nu_v * Ev
##   Sv -> (outside)  : mu_v * Sv
##   Ev -> (outside)  : mu_v * Ev
##   Iv -> (outside)  : mu_v * Iv

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 365, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  a_bite = 0.5,
  b_vh = 0.4,       # P(mosquito → human) per bite
  p_hv = 0.5,       # P(symptomatic human → mosquito) per bite
  eta = 0.1,        # relative human → mosquito infectivity of Eh
  bet_sex = 0.05,   # sexual transmission rate from symptomatic humans
  kappa = 0.6,      # relative sexual infectivity of Eh
  tau_c = 0.3,      # relative sexual infectivity of convalescent I2
  theta = 0.18,     # proportion symptomatic
  m_ratio = 5,
  nu_h = 0.2,       # 1/(5 d) intrinsic incubation [2–7 d] (1/day)
  nu_v = 0.1,       # 1/(10 d) extrinsic incubation [8–12 d] (1/day)
  gam_h1 = 0.2,     # 1/(5 d) acute phase [3–7 d] (1/day)
  gam_h2 = 0.05,    # 1/(20 d) convalescent phase [14–30 d] (1/day)
  gam_ha = 0.1429,  # 1/(7 d) asymptomatic infection [5–10 d] (1/day)
  mu_v = 0.0714     # 1/(14 d) mosquito life [4–35 d] (1/day)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sh = 99990,   # susceptible humans
  Eh = 0,       # exposed humans (mildly infectious to mosquitoes and sexually)
  I1 = 10,      # symptomatic acute
  I2 = 0,       # convalescent (sexually infectious only)
  Ah = 0,       # asymptomatic
  Rh = 0,       # recovered
  Sv = 500000,  # susceptible mosquitoes
  Ev = 0,       # exposed mosquitoes
  Iv = 0        # infectious mosquitoes
)

# ---- Model equations -------------------------------------------------
zika_gao2016_sexual_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sh + Eh + I1 + I2 + Ah + Rh

    # Derived quantities
    Nh <- Sh + Eh + I1 + I2 + Ah + Rh
    Nv <- Sv + Ev + Iv
    lambda_h <- a_bite * b_vh * Iv / Nh + bet_sex * (kappa * Eh + I1 + tau_c * I2) / Nh  # vector-borne plus sexual force of infection
    lambda_v <- a_bite * p_hv * (eta * Eh + I1) / Nh                                     # mosquito force of infection

    # Rates of change
    dSh <- -lambda_h * Sh
    dEh <- lambda_h * Sh - theta * nu_h * Eh - (1 - theta) * nu_h * Eh
    dI1 <- theta * nu_h * Eh - gam_h1 * I1
    dI2 <- gam_h1 * I1 - gam_h2 * I2
    dAh <- (1 - theta) * nu_h * Eh - gam_ha * Ah
    dRh <- gam_h2 * I2 + gam_ha * Ah
    dSv <- mu_v * Nv - lambda_v * Sv - mu_v * Sv
    dEv <- lambda_v * Sv - nu_v * Ev - mu_v * Ev
    dIv <- nu_v * Ev - mu_v * Iv

    list(c(dSh, dEh, dI1, dI2, dAh, dRh, dSv, dEv, dIv))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = zika_gao2016_sexual_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sh = "#333333", Eh = "#E8C468", I1 = "#B71C1C", I2 = "#B71C1C", Ah = "#E57373", Rh = "#4E8A3E", Sv = "#333333", Ev = "#E8C468", Iv = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Zika with mosquito-borne and sexual transmission (Gao et al. 2016) [VERIFIED-PARAMS]")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
