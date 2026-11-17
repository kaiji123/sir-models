## Mumps SEIR with vaccine-derived immunity that wanes (Lewnard & Grad 2018)
## Mumps.
## the booster flow acts on all of S, so it also reaches never-vaccinated people.
## Parameter values: 0 of 8 checked against the source; the rest are recalled or illustrative.
## Reference: Lewnard JA, Grad YH. "Vaccine waning and mumps re-emergence in the United States." Sci Transl Med 2018;10(433):eaao5945. https://doi.org/10.1126/scitranslmed.aao5945
## Source: https://doi.org/10.1126/scitranslmed.aao5945
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, Ip, R, V
## Flows (time unit: years):
##   (outside) -> S  : (1 - p_vac) * mu * N
##   (outside) -> V  : p_vac * mu * N
##   S -> E          : foi * S
##   E -> Ip         : sig_lat * E
##   Ip -> R         : g_rec * Ip
##   V -> S          : w_vac * V
##   S -> V          : boost_v * S
##   S -> (outside)  : mu * S

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 100, by = 0.1)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  R0_mu = 7,           # illustrative
  sig_lat = 365 / 17,  # 17-day latent period (1/year) [illustrative]
  g_rec = 365 / 8,     # 8-day infectious period (1/year) [illustrative]
  p_vac = 0.92,        # fraction of births effectively vaccinated (coverage × efficacy) [illustrative]
  w_vac = 1 / 27,      # waning of vaccine immunity (1/year)
  boost_v = 0,         # rate of a booster (third) dose that returns S to V; set > 0 for scenarios (1/year) [illustrative]
  mu = 0.0133,         # birth = death rate (1/year) [illustrative]
  N0 = 1e6             # value taken from the notes text [recalled]
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  S <- 0.05 * N0  # susceptible
  E <- 10         # exposed
  Ip <- 10        # infectious
  R <- 0.3 * N0   # infection-derived immunity, lifelong
  V <- 0.65 * N0  # vaccine-derived immunity
  c(S = S, E = E, Ip = Ip, R = R, V = V)
})

# ---- Model equations -------------------------------------------------
mumps_waning_lewnard_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + Ip + R + V

    # Derived quantities
    foi <- R0_mu * (g_rec + mu) * Ip / N

    # Rates of change
    dS <- (1 - p_vac) * mu * N - foi * S + w_vac * V - boost_v * S - mu * S
    dE <- foi * S - sig_lat * E
    dIp <- sig_lat * E - g_rec * Ip
    dR <- g_rec * Ip
    dV <- p_vac * mu * N - w_vac * V + boost_v * S

    list(c(dS, dE, dIp, dR, dV))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = mumps_waning_lewnard_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", Ip = "#B71C1C", R = "#4E8A3E", V = "#2457C5")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Mumps SEIR with vaccine-derived immunity that wanes (Lewnard & Grad 2018)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
