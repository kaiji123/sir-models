## Measles MSEIR with two-dose vaccination (dose-1 failures eligible for dose 2)
## Measles.
## Parameter values: 0 of 10 checked against the source; the rest are recalled or illustrative.
## Reference: [ILLUSTRATIVE] MSEIR structure from Hethcote HW, "The mathematics of infectious diseases," SIAM Review 2000;42:599–653, https://doi.org/10.1137/S0036144500371907. The values below are from my recall and were not fetched this session:
## Source: https://doi.org/10.1137/S0036144500371907
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: M, S0, Sf, V, E, Ip, R
## Flows (time unit: years):
##   (outside) -> M  : mu * N
##   M -> S0         : xi_m * M
##   S0 -> V         : e1 * nu1 * S0
##   S0 -> Sf        : (1 - e1) * nu1 * S0
##   Sf -> V         : e2 * nu2 * Sf
##   S0 -> E         : foi * S0
##   Sf -> E         : foi * Sf
##   E -> Ip         : sig_lat * E
##   Ip -> R         : g_rec * Ip
##   M -> (outside)  : mu * M

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 50, by = 0.05)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  R0_m = 15,          # Guerra 2017 range (recall)
  sig_lat = 365 / 8,  # 8-day latent period (1/year) [illustrative]
  g_rec = 365 / 5,    # 5-day infectious period (1/year) [illustrative]
  xi_m = 2,           # loss of maternal antibody at about 6 months (1/year) [illustrative]
  nu1 = 1.5,          # per-capita dose-1 rate in S0 (1/year) [illustrative]
  nu2 = 0.3,          # dose-2 rate in Sf (1/year) [illustrative]
  e1 = 0.93,          # efficacy of dose 1, and of dose 2 given to a dose-1 failure
  e2 = 0.97,          # efficacy of dose 1, and of dose 2 given to a dose-1 failure
  mu = 0.0133,        # birth = death rate (1/year) [illustrative]
  N0 = 1e6            # value taken from the notes text [recalled]
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  M <- 0.01 * N0   # maternally protected
  S0 <- 0.05 * N0  # susceptible, unvaccinated
  Sf <- 0.01 * N0  # susceptible after failed or missed dose 1, eligible for dose 2
  V <- 0.6 * N0    # vaccine-immune
  E <- 10          # exposed
  Ip <- 10         # infectious
  R <- 0.33 * N0   # recovered
  c(M = M, S0 = S0, Sf = Sf, V = V, E = E, Ip = Ip, R = R)
})

# ---- Model equations -------------------------------------------------
measles_mseir_2dose_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- M + S0 + Sf + V + E + Ip + R

    # Derived quantities
    foi <- R0_m * (g_rec + mu) * Ip / N

    # Rates of change
    dM <- mu * N - xi_m * M - mu * M
    dS0 <- xi_m * M - e1 * nu1 * S0 - (1 - e1) * nu1 * S0 - foi * S0
    dSf <- (1 - e1) * nu1 * S0 - e2 * nu2 * Sf - foi * Sf
    dV <- e1 * nu1 * S0 + e2 * nu2 * Sf
    dE <- foi * S0 + foi * Sf - sig_lat * E
    dIp <- sig_lat * E - g_rec * Ip
    dR <- g_rec * Ip

    list(c(dM, dS0, dSf, dV, dE, dIp, dR))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = measles_mseir_2dose_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(M = "#F3C6D3", S0 = "#333333", Sf = "#333333", V = "#2457C5", E = "#E8C468", Ip = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Measles MSEIR with two-dose vaccination (dose-1 failures eligible for dose 2)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
