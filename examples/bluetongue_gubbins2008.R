## Bluetongue: cattle + sheep + Culicoides with temperature-dependent biting, mortality and EIP (after Gubbins et al. 2008) [STRUCTURE-FROM-RECALL]
## Bluetongue virus (BTV-8), northern Europe / UK.
## Parameter values: 0 of 10 checked against the source; the rest are recalled or illustrative.
## Reference: Gubbins S, Carpenter S, Baylis M, Wood JLN, Mellor PS (2008). Assessing the risk of bluetongue to UK livestock: uncertainty and sensitivity analyses of a temperature-dependent model for the basic reproduction number. J R Soc Interface 5(20):363–371.…
## Source: https://doi.org/10.1098/rsif.2007.1110
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Scow, Icow, Rcow, Ssh, Ish, Rsh, Sv, Ev1, Ev2, Ev3, Iv
## Flows (time unit: days):
##   Scow -> Icow      : lam_cow * Scow
##   Icow -> Rcow      : r_cow * Icow
##   Ssh -> Ish        : lam_sh * Ssh
##   Ish -> Rsh        : r_sh * Ish
##   Ish -> (outside)  : d_sh * Ish
##   (outside) -> Sv   : mu_v * Nv
##   Sv -> Ev1         : lam_v * Sv
##   Ev1 -> Ev2        : 3 * nu_eip * Ev1
##   Ev2 -> Ev3        : 3 * nu_eip * Ev2
##   Ev3 -> Iv         : 3 * nu_eip * Ev3
##   Sv -> (outside)   : mu_v * Sv
##   Ev1 -> (outside)  : mu_v * Ev1
##   Ev2 -> (outside)  : mu_v * Ev2
##   Ev3 -> (outside)  : mu_v * Ev3
##   Iv -> (outside)   : mu_v * Iv

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 365, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  T_mean = 15,     # illustrative (°C)
  T_amp = 7,       # illustrative (°C)
  t_shift = 110,   # sets the temperature peak in mid-July (day)
  b_vh = 0.9,      # vector → host per-bite probability
  p_hv = 0.05,     # host → vector per-bite probability (recalled range 0.01–0.15)
  phi_c = 0.75,    # biting preference weight for cattle (illustrative)
  r_cow = 0.0485,  # 1/20.6 d cattle viraemia (recalled) (1/day)
  r_sh = 0.0625,   # 1/16 d sheep viraemia (recalled) (1/day)
  d_sh = 0.006,    # BTV-induced sheep death (illustrative) (1/day)
  k_eip = 3        # EIP stages (fixed in the flows)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Scow = 199,   # susceptible cattle
  Icow = 1,     # infected cattle
  Rcow = 0,     # recovered cattle
  Ssh = 800,    # susceptible sheep
  Ish = 0,      # infected sheep
  Rsh = 0,      # recovered sheep
  Sv = 100000,  # susceptible Culicoides
  Ev1 = 0,      # EIP stage 1
  Ev2 = 0,      # EIP stage 2
  Ev3 = 0,      # EIP stage 3
  Iv = 0        # infectious midges
)

# ---- Model equations -------------------------------------------------
bluetongue_gubbins2008_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Scow + Icow + Rcow + Ssh + Ish + Rsh

    # Derived quantities
    Temp <- T_mean + T_amp * sin(2 * 3.14159265 * (t - t_shift) / 365)                # daily temperature (°C)
    a_bite <- pmax(0.0002 * Temp * (Temp - 3.7) * pmax(41.9 - Temp, 0)^(1 / 2.7), 0)  # biting rate (1/day)
    mu_v <- 0.009 * exp(0.16 * Temp)                                                  # midge mortality (1/day)
    nu_eip <- pmax(0.018 * (Temp - 13.4), 0)                                          # EIP completion rate (1/day)
    Ncow <- Scow + Icow + Rcow
    Nsh <- Ssh + Ish + Rsh
    Nv <- Sv + Ev1 + Ev2 + Ev3 + Iv
    Hw <- phi_c * Ncow + (1 - phi_c) * Nsh                                            # preference-weighted host total
    lam_cow <- a_bite * b_vh * phi_c * Iv / Hw                                        # cattle force of infection
    lam_sh <- a_bite * b_vh * (1 - phi_c) * Iv / Hw                                   # sheep force of infection
    lam_v <- a_bite * p_hv * (phi_c * Icow + (1 - phi_c) * Ish) / Hw                  # midge force of infection

    # Rates of change
    dScow <- -lam_cow * Scow
    dIcow <- lam_cow * Scow - r_cow * Icow
    dRcow <- r_cow * Icow
    dSsh <- -lam_sh * Ssh
    dIsh <- lam_sh * Ssh - r_sh * Ish - d_sh * Ish
    dRsh <- r_sh * Ish
    dSv <- mu_v * Nv - lam_v * Sv - mu_v * Sv
    dEv1 <- lam_v * Sv - 3 * nu_eip * Ev1 - mu_v * Ev1
    dEv2 <- 3 * nu_eip * Ev1 - 3 * nu_eip * Ev2 - mu_v * Ev2
    dEv3 <- 3 * nu_eip * Ev2 - 3 * nu_eip * Ev3 - mu_v * Ev3
    dIv <- 3 * nu_eip * Ev3 - mu_v * Iv

    list(c(dScow, dIcow, dRcow, dSsh, dIsh, dRsh, dSv, dEv1, dEv2, dEv3, dIv))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = bluetongue_gubbins2008_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Scow = "#333333", Icow = "#B71C1C", Rcow = "#4E8A3E", Ssh = "#333333", Ish = "#B71C1C", Rsh = "#4E8A3E", Sv = "#333333", Ev1 = "#E8C468", Ev2 = "#E8C468", Ev3 = "#E8C468", Iv = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Bluetongue: cattle + sheep + Culicoides with temperature-dependent biting, mortality and EIP (after Gubbins et al. 2008) [STRUCTURE-FROM-RECALL]")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
