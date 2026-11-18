## Pertussis with vaccination, waning, and mild reinfections in partially immune people (Hethcote-style)
## Pertussis.
## Hethcote's 1997 paper is age-structured with many age groups.
## Parameter values: 0 of 8 checked against the source; the rest are recalled or illustrative.
## Reference: Hethcote HW. "An age-structured model for pertussis transmission." Math Biosci 1997;145:89–136. https://doi.org/10.1016/S0025-5564(97)00014-X
## Source: https://doi.org/10.1016/S0025-5564(97)00014-X
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, Ip, R, V, W, Im
## Flows (time unit: years):
##   (outside) -> S  : (1 - p_vac) * mu * N
##   (outside) -> V  : p_vac * mu * N
##   S -> Ip         : foi * S
##   Ip -> R         : g_rec * Ip
##   R -> W          : w_nat * R
##   V -> W          : w_vac * V
##   W -> Im         : foi * W
##   Im -> R         : g_rec * Im
##   S -> (outside)  : mu * S

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 100, by = 0.1)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  R0_p = 15,       # illustrative; the literature for pertussis spans roughly 5–17
  g_rec = 17.4,    # 21-day infectious period (1/year) [illustrative]
  rel_m = 0.5,     # relative infectiousness of mild cases [illustrative]
  w_nat = 1 / 30,  # waning of natural immunity (1/year) [illustrative]
  w_vac = 1 / 10,  # waning of vaccine immunity (1/year) [illustrative]
  p_vac = 0.9,     # fraction of births effectively vaccinated [illustrative]
  mu = 0.0133,     # birth = death rate (1/year) [illustrative]
  N0 = 1e6         # value taken from the notes text [recalled]
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  S <- 0.05 * N0  # fully susceptible (never infected or vaccinated)
  Ip <- 50        # typical (severe) primary infection
  R <- 0.5 * N0   # fully immune after infection
  V <- 0.3 * N0   # vaccine-immune
  W <- 0.15 * N0  # waned, partially immune
  Im <- 50        # mild infection in a waned person
  c(S = S, Ip = Ip, R = R, V = V, W = W, Im = Im)
})

# ---- Model equations -------------------------------------------------
pertussis_hethcote_mild_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + Ip + R + V + W + Im

    # Derived quantities
    foi <- R0_p * (g_rec + mu) * (Ip + rel_m * Im) / N

    # Rates of change
    dS <- (1 - p_vac) * mu * N - foi * S - mu * S
    dIp <- foi * S - g_rec * Ip
    dR <- g_rec * Ip - w_nat * R + g_rec * Im
    dV <- p_vac * mu * N - w_vac * V
    dW <- w_nat * R + w_vac * V - foi * W
    dIm <- foi * W - g_rec * Im

    list(c(dS, dIp, dR, dV, dW, dIm))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = pertussis_hethcote_mild_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", Ip = "#B71C1C", R = "#4E8A3E", V = "#2457C5", W = "#A5D6A7", Im = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Pertussis with vaccination, waning, and mild reinfections in partially immune people (Hethcote-style)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
