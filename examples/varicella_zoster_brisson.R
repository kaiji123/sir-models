## Varicella with latent VZV, exogenous boosting and zoster reactivation (after Brisson et al. 2000)
## Varicella / herpes zoster.
## Parameter values: 0 of 12 checked against the source; the rest are recalled or illustrative.
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: M, S, E, Ip, Lb, Zs, Zo, Rz, V
## Flows (time unit: years):
##   (outside) -> M  : (1 - p_vac) * mu * N
##   (outside) -> V  : p_vac * mu * N
##   M -> S          : xi_m * M
##   S -> E          : foi * S
##   E -> Ip         : sig_lat * E
##   Ip -> Lb        : g_rec * Ip
##   Lb -> Zs        : w_boost * Lb
##   Zs -> Lb        : boost_f * foi * Zs
##   Zs -> Zo        : react * Zs
##   Zo -> Rz        : z_rec * Zo
##   M -> (outside)  : mu * M

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 100, by = 0.1)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  R0_v = 8,            # illustrative (varicella R0 is often given as 7–12)
  sig_lat = 365 / 14,  # 14-day latent period (1/year) [illustrative]
  g_rec = 365 / 7,     # 7-day infectious period (1/year) [illustrative]
  xi_m = 2,            # loss of maternal protection (1/year) [illustrative]
  w_boost = 1 / 20,    # loss of boosted protection (1/year)
  boost_f = 1,         # fraction of varicella exposures that boost latent people [illustrative]
  react = 0.004,       # zoster reactivation rate in Zs (1/year) [illustrative]
  z_rec = 365 / 21,    # zoster duration (1/year) [illustrative]
  rel_z = 0.2,         # infectiousness of zoster relative to varicella [illustrative]
  p_vac = 0,           # fraction of births vaccinated; vary between 0 and 0.9 [illustrative]
  mu = 0.0133,         # birth = death rate (1/year) [illustrative]
  N0 = 1e6             # value taken from the notes text [recalled]
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  M <- 0.005 * N0  # maternally protected
  S <- 0.08 * N0   # susceptible to varicella
  E <- 0           # latent (pre-infectious) varicella
  Ip <- 100        # varicella, infectious
  Lb <- 0.5 * N0   # VZV-latent, boosted (protected from zoster)
  Zs <- 0.41 * N0  # VZV-latent, susceptible to zoster
  Zo <- 0          # active zoster
  Rz <- 0          # post-zoster, immune
  V <- 0           # varicella-vaccinated, immune; vaccine-strain latency ignored
  c(M = M, S = S, E = E, Ip = Ip, Lb = Lb, Zs = Zs, Zo = Zo, Rz = Rz, V = V)
})

# ---- Model equations -------------------------------------------------
varicella_zoster_brisson_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- M + S + E + Ip + Lb + Zs + Zo + Rz + V

    # Derived quantities
    foi <- R0_v * (g_rec + mu) * (Ip + rel_z * Zo) / N

    # Rates of change
    dM <- (1 - p_vac) * mu * N - xi_m * M - mu * M
    dS <- xi_m * M - foi * S
    dE <- foi * S - sig_lat * E
    dIp <- sig_lat * E - g_rec * Ip
    dLb <- g_rec * Ip - w_boost * Lb + boost_f * foi * Zs
    dZs <- w_boost * Lb - boost_f * foi * Zs - react * Zs
    dZo <- react * Zs - z_rec * Zo
    dRz <- z_rec * Zo
    dV <- p_vac * mu * N

    list(c(dM, dS, dE, dIp, dLb, dZs, dZo, dRz, dV))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = varicella_zoster_brisson_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(M = "#F3C6D3", S = "#333333", E = "#E8C468", Ip = "#B71C1C", Lb = "#E8C468", Zs = "#333333", Zo = "#333333", Rz = "#4E8A3E", V = "#2457C5")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Varicella with latent VZV, exogenous boosting and zoster reactivation (after Brisson et al. 2000)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
