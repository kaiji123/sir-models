## Rift Valley fever: Aedes (with transovarial transmission and eggs), Culex and livestock (Gaff, Hartley & Leahy 2007) [STRUCTURE-FROM-RECALL]
## Rift Valley fever, East Africa (Kenya-type floodplain).
## Parameter values: 0 of 24 checked against the source; the rest are recalled or illustrative.
## Reference: Gaff HD, Hartley DM, Leahy NP (2007). An epidemiological model of Rift Valley fever. Electron J Differ Equ 2007(115):1–12. Open access: https://ejde.math.txstate.edu/Volumes/2007/115/abstr.html . The abstract confirms two mosquito populations (Aedes, which…
## Source: https://ejde.math.txstate.edu/Volumes/2007/115/abstr.html
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Pa, Qa, Sa, Ea, Ia, Pc, Sc, Ec, Ic, Sl, El, Il, Rl
## Flows (time unit: days):
##   (outside) -> Pa  : (b_a * (Sa + Ea) + b_a * (1 - q_tov) * Ia)
##   (outside) -> Qa  : b_a * q_tov * Ia
##   Pa -> Sa         : th_a * flood * pmax(0, 1 - Nas / K_a) * Pa
##   Qa -> Ia         : th_a * flood * pmax(0, 1 - Nas / K_a) * Qa
##   Pa -> (outside)  : mu_ea * Pa
##   Qa -> (outside)  : mu_ea * Qa
##   Sa -> Ea         : lam_a * Sa
##   Ea -> Ia         : eps_a * Ea
##   Sa -> (outside)  : mu_a * Sa
##   Ea -> (outside)  : mu_a * Ea
##   Ia -> (outside)  : mu_a * Ia
##   (outside) -> Pc  : b_c * Ncs
##   Pc -> Sc         : th_c * flood * pmax(0, 1 - Ncs / K_c) * Pc
##   Pc -> (outside)  : mu_ec * Pc
##   Sc -> Ec         : lam_c * Sc
##   Ec -> Ic         : eps_c * Ec
##   Sc -> (outside)  : mu_c * Sc
##   Ec -> (outside)  : mu_c * Ec
##   Ic -> (outside)  : mu_c * Ic
##   (outside) -> Sl  : mu_l * Nl
##   Sl -> El         : lam_l * Sl
##   El -> Il         : eps_l * El
##   Il -> Rl         : gam_l * Il
##   Il -> (outside)  : d_l * Il
##   Sl -> (outside)  : mu_l * Sl
##   El -> (outside)  : mu_l * El
##   Il -> (outside)  : mu_l * Il
##   Rl -> (outside)  : mu_l * Rl

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 730, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  b_a = 0.2,      # (eggs/female/day)
  q_tov = 0.05,   # transovarial fraction of eggs from infectious Aedes
  th_a = 0.05,    # egg hatching rate (multiplied by the flood indicator) (1/day)
  mu_ea = 0.001,  # Aedes egg mortality; desiccation-resistant eggs (1/day)
  K_a = 2e5,      # Aedes adult capacity (mosquitoes)
  mu_a = 0.05,    # Aedes adult mortality (1/day)
  eps_a = 0.143,  # 1/7 d EIP in Aedes (1/day)
  a_a = 0.3,      # Aedes biting rate (1/day)
  b_c = 0.2,      # (eggs/female/day)
  th_c = 0.1,     # Culex hatching rate (1/day)
  mu_ec = 0.1,    # Culex egg mortality (1/day)
  K_c = 2e5,      # Culex adult capacity (mosquitoes)
  mu_c = 0.05,    # Culex adult mortality (1/day)
  eps_c = 0.143,  # 1/7 d EIP in Culex (1/day)
  a_c = 0.3,      # Culex biting rate (1/day)
  p_al = 0.3,     # Aedes → livestock per-bite probability
  p_cl = 0.3,     # Culex → livestock per-bite probability
  p_la = 0.3,     # livestock → Aedes per-bite probability
  p_lc = 0.3,     # livestock → Culex per-bite probability
  eps_l = 0.25,   # 1/4 d incubation in livestock (1/day)
  gam_l = 0.25,   # 1/4 d infectious period (1/day)
  d_l = 0.025,    # RVF-induced livestock death (1/day)
  mu_l = 5.5e-4,  # livestock death (5-yr life); births = mu_l·Nl (1/day)
  flood_len = 60  # wet-season window each year (days)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Pa = 1e6,   # uninfected Aedes eggs
  Qa = 0,     # infected Aedes eggs (transovarial)
  Sa = 1e5,   # susceptible Aedes
  Ea = 0,     # exposed Aedes
  Ia = 10,    # infectious Aedes
  Pc = 1e5,   # Culex eggs
  Sc = 1e5,   # susceptible Culex
  Ec = 0,     # exposed Culex
  Ic = 0,     # infectious Culex
  Sl = 1000,  # susceptible livestock
  El = 0,     # exposed livestock
  Il = 0,     # infectious livestock
  Rl = 0      # recovered livestock
)

# ---- Model equations -------------------------------------------------
rvf_gaff2007_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sl + El + Il + Rl

    # Derived quantities
    Nas <- Sa + Ea + Ia                                          # adult Aedes
    Ncs <- Sc + Ec + Ic                                          # adult Culex
    Nl <- Sl + El + Il + Rl                                      # livestock
    flood <- ifelse(t - 365 * floor(t / 365) < flood_len, 1, 0)  # annual flooding window
    lam_l <- (a_a * p_al * Ia + a_c * p_cl * Ic) / Nl            # livestock force of infection
    lam_a <- a_a * p_la * Il / Nl                                # Aedes force of infection
    lam_c <- a_c * p_lc * Il / Nl                                # Culex force of infection

    # Rates of change
    dPa <- (b_a * (Sa + Ea) + b_a * (1 - q_tov) * Ia) - th_a * flood * pmax(0, 1 - Nas / K_a) * Pa - mu_ea * Pa
    dQa <- b_a * q_tov * Ia - th_a * flood * pmax(0, 1 - Nas / K_a) * Qa - mu_ea * Qa
    dSa <- th_a * flood * pmax(0, 1 - Nas / K_a) * Pa - lam_a * Sa - mu_a * Sa
    dEa <- lam_a * Sa - eps_a * Ea - mu_a * Ea
    dIa <- th_a * flood * pmax(0, 1 - Nas / K_a) * Qa + eps_a * Ea - mu_a * Ia
    dPc <- b_c * Ncs - th_c * flood * pmax(0, 1 - Ncs / K_c) * Pc - mu_ec * Pc
    dSc <- th_c * flood * pmax(0, 1 - Ncs / K_c) * Pc - lam_c * Sc - mu_c * Sc
    dEc <- lam_c * Sc - eps_c * Ec - mu_c * Ec
    dIc <- eps_c * Ec - mu_c * Ic
    dSl <- mu_l * Nl - lam_l * Sl - mu_l * Sl
    dEl <- lam_l * Sl - eps_l * El - mu_l * El
    dIl <- eps_l * El - gam_l * Il - d_l * Il - mu_l * Il
    dRl <- gam_l * Il - mu_l * Rl

    list(c(dPa, dQa, dSa, dEa, dIa, dPc, dSc, dEc, dIc, dSl, dEl, dIl, dRl))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = rvf_gaff2007_model, parms = parameters, hmax = 1)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Pa = "#F0A35E", Qa = "#00838F", Sa = "#333333", Ea = "#E8C468", Ia = "#B71C1C", Pc = "#F0A35E", Sc = "#333333", Ec = "#E8C468", Ic = "#B71C1C", Sl = "#333333", El = "#E8C468", Il = "#B71C1C", Rl = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Rift Valley fever: Aedes (with transovarial transmission and eggs), Culex and livestock (Gaff, Hartley & Leahy 2007) [STRUCTURE-FROM-RECALL]")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
