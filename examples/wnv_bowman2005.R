## West Nile virus: mosquitoes–birds–humans (Bowman et al. 2005) [VERIFIED-FULL]
## West Nile virus, North America (single season).
## Parameter values: 0 of 18 checked against the source; the rest are recalled or illustrative.
## Reference: Bowman C, Gumel AB, van den Driessche P, Wu J, Zhu H (2005). A mathematical model for assessing control strategies against West Nile virus. Bull Math Biol 67:1107–1133. https://doi.org/10.1016/j.bulm.2005.01.002 ; full text in the YorkSpace repository. The…
## Source: https://doi.org/10.1016/j.bulm.2005.01.002
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Mu, Mi, Bu, Bi, Sh, Eh, Ih, Hh, Rh
## Flows (time unit: days):
##   (outside) -> Mu  : Pi_M
##   Mu -> Mi         : lam_M * Mu
##   Mu -> (outside)  : mu_M * Mu
##   Mi -> (outside)  : mu_M * Mi
##   (outside) -> Bu  : Pi_B
##   Bu -> Bi         : lam_B * Bu
##   Bu -> (outside)  : mu_B * Bu
##   Bi -> (outside)  : mu_B * Bi
##   Bi -> (outside)  : d_B * Bi
##   (outside) -> Sh  : Pi_H
##   Sh -> Eh         : lam_H * Sh
##   Eh -> Ih         : alpha_h * Eh
##   Ih -> Hh         : delta_h * Ih
##   Hh -> Rh         : tau_h * Hh
##   Hh -> (outside)  : d_H * Hh
##   Sh -> (outside)  : mu_H * Sh
##   Eh -> (outside)  : mu_H * Eh
##   Ih -> (outside)  : mu_H * Ih
##   Hh -> (outside)  : mu_H * Hh
##   Rh -> (outside)  : mu_H * Rh

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 0.2)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  Pi_M = 20000,      # mosquito recruitment. "Variable" in the paper; this value is illustrative. (mosquitoes/day)
  mu_M = 0.05,       # mosquito death. "Variable" in the paper; illustrative 20-day life. (1/day)
  Pi_B = 1000,       # bird recruitment (birds/day)
  mu_B = 0.001,      # 1/(1000 d) bird life (1/day)
  Pi_H = 30,         # human recruitment (humans/day)
  mu_H = 3.914e-5,   # 1/(70×365 d) (1/day)
  b1 = 0.09,         # biting rate on birds (1/day)
  b2_base = 0.09,    # biting rate on humans before protection (1/day)
  comp_c = 0,        # fraction complying with personal protection (0–1)
  q_eff = 0,         # efficacy of personal protection (0–1)
  bet1 = 0.16,       # bird → mosquito transmission probability
  bet2 = 0.88,       # mosquito → bird transmission probability
  bet3 = 0.88,       # mosquito → human transmission probability
  d_B = 5e-5,        # WNV-induced bird death (1/day)
  d_H = 5e-7,        # WNV-induced human death (1/day)
  alpha_h = 0.0714,  # 1/(14 d) incubation in humans (1/day)
  delta_h = 1,       # hospitalisation rate (1/day)
  tau_h = 0.0714     # 1/(14 d) recovery (1/day)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Mu = 400000,  # uninfected female mosquitoes
  Mi = 0,       # infected (infectious) mosquitoes
  Bu = 999990,  # uninfected birds
  Bi = 10,      # infected birds
  Sh = 766500,  # susceptible humans
  Eh = 0,       # infected, asymptomatic (incubating)
  Ih = 0,       # symptomatic
  Hh = 0,       # hospitalised
  Rh = 0        # recovered (lifelong immunity)
)

# ---- Model equations -------------------------------------------------
wnv_bowman2005_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sh + Eh + Ih + Hh + Rh

    # Derived quantities
    NM <- Mu + Mi
    NB <- Bu + Bi
    NH <- Sh + Eh + Ih + Hh + Rh
    b2 <- b2_base * (1 - comp_c * q_eff)  # biting rate on humans after protection
    lam_M <- b1 * bet1 * Bi / NB          # mosquito force of infection
    lam_B <- b1 * bet2 * Mi / NB          # bird force of infection
    lam_H <- b2 * bet3 * Mi / NH          # human force of infection

    # Rates of change
    dMu <- Pi_M - lam_M * Mu - mu_M * Mu
    dMi <- lam_M * Mu - mu_M * Mi
    dBu <- Pi_B - lam_B * Bu - mu_B * Bu
    dBi <- lam_B * Bu - mu_B * Bi - d_B * Bi
    dSh <- Pi_H - lam_H * Sh - mu_H * Sh
    dEh <- lam_H * Sh - alpha_h * Eh - mu_H * Eh
    dIh <- alpha_h * Eh - delta_h * Ih - mu_H * Ih
    dHh <- delta_h * Ih - tau_h * Hh - d_H * Hh - mu_H * Hh
    dRh <- tau_h * Hh - mu_H * Rh

    list(c(dMu, dMi, dBu, dBi, dSh, dEh, dIh, dHh, dRh))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = wnv_bowman2005_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Mu = "#F3C6D3", Mi = "#F3C6D3", Bu = "#333333", Bi = "#333333", Sh = "#333333", Eh = "#E8C468", Ih = "#B71C1C", Hh = "#6A2C9C", Rh = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "West Nile virus: mosquitoes–birds–humans (Bowman et al. 2005) [VERIFIED-FULL]")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
