## Haiti cholera with water reservoir, vaccination and water-sanitation intervention
## Cholera, Haiti 2010–11.
## Parameter values: 0 of 15 checked against the source; the rest are recalled or illustrative.
## Reference: Andrews JR, Basu S. Transmission dynamics and control of cholera in Haiti: an epidemic model. Lancet 2011;377(9773):1248–1255
## Source: https://doi.org/10.1016/S0140-6736(11
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, A, R, V, W
## Flows (time unit: days):
##   S -> I          : p_sym * lam * S
##   S -> A          : (1 - p_sym) * lam * S
##   I -> R          : gam_c * I
##   I -> (outside)  : cfr_rate * I
##   A -> R          : gam_c * A
##   R -> S          : omega_r * R
##   S -> V          : nu_v * eff_v * S
##   V -> S          : omega_v * V
##   (outside) -> W  : xi_w * (I + rel_a * A)
##   W -> (outside)  : del_w * W

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 300, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta_w = 1.0,        # contaminated-water exposure rate (fitted in the paper) [I] (day^-1)
  beta_h = 0.1,        # direct (short-cycle) transmission [I] (day^-1)
  kap_w = 1e6,         # half-saturation dose, as in Codeço [M] (cells/ml)
  p_sym = 0.2,         # fraction of infections that are symptomatic [I]
  rel_a = 0.1,         # relative shedding/infectiousness of asymptomatic infections [I]
  gam_c = 0.2,         # recovery rate (5 days) [M] (day^-1)
  cfr_rate = 0.002,    # symptomatic death rate (fitted in the paper) [I] (day^-1)
  xi_w = 10,           # shedding [I] (cells ml^-1 day^-1 per symptomatic person)
  del_w = 0.033,       # vibrio decay rate (30-day survival) [M, as used in Hartley/Morris/Smith 2006-style models] (day^-1)
  omega_r = 1 / 1095,  # waning of natural immunity (about 3 years) [M] (day^-1)
  nu_v = 0,            # vaccination rate; try 0.001 [I] (day^-1)
  eff_v = 0.65,        # vaccine efficacy [I]
  omega_v = 1 / 1095,  # vaccine waning [I] (day^-1)
  red_w = 0.001436,    # equals −log(0.99)/7, a 1%/week reduction in contaminated-water consumption [V for the 1%/week scenario] (day^-1)
  t_int = 120          # start of the water intervention [I] (day)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9,     # susceptible
  I = 1000,  # symptomatic infected
  A = 0,     # asymptomatic infected
  R = 0,     # recovered/immune
  V = 0,     # vaccinated (OCV)
  W = 1e3    # V. cholerae concentration in the water reservoir (cells/ml)
)

# ---- Model equations -------------------------------------------------
cholera_andrews_basu2011_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + A + R + V

    # Derived quantities
    water_mult <- exp(-red_w * max(t - t_int, 0))
    lam <- beta_w * water_mult * W / (kap_w + W) + beta_h * (I + rel_a * A) / N

    # Rates of change
    dS <- -p_sym * lam * S - (1 - p_sym) * lam * S + omega_r * R - nu_v * eff_v * S + omega_v * V
    dI <- p_sym * lam * S - gam_c * I - cfr_rate * I
    dA <- (1 - p_sym) * lam * S - gam_c * A
    dR <- gam_c * I + gam_c * A - omega_r * R
    dV <- nu_v * eff_v * S - omega_v * V
    dW <- xi_w * (I + rel_a * A) - del_w * W

    list(c(dS, dI, dA, dR, dV, dW))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = cholera_andrews_basu2011_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", A = "#E57373", R = "#4E8A3E", V = "#2457C5", W = "#A5D6A7")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Haiti cholera with water reservoir, vaccination and water-sanitation intervention")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
