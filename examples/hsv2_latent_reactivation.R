## Genital herpes (HSV-2) with primary episode, latency and recurrent reactivation
## HSV-2.
## Parameter values: 0 of 9 checked against the source; the rest are recalled or illustrative.
## Reference: Structure follows Blower SM, Porco TC, Darby G. Predicting and preventing the emergence of antiviral drug resistance in HSV-2. Nat Med 1998;4(6):673–678
## Source: https://doi.org/10.1038/nm0698-673
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, Ip, Lt, Rc
## Flows (time unit: years):
##   (outside) -> S   : mu * N
##   S -> Ip          : foi * S
##   Ip -> Lt         : 1 / d_p * Ip
##   Lt -> Rc         : react * (1 - tx) * Lt
##   Rc -> Lt         : 1 / d_r * Rc
##   S -> (outside)   : mu * S
##   Ip -> (outside)  : mu * Ip
##   Lt -> (outside)  : mu * Lt
##   Rc -> (outside)  : mu * Rc

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 50, by = 0.05)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  cp = 2,         # partner change rate (partners/yr) [illustrative]
  beta_p = 0.5,   # transmission from primary episode (per partnership) [illustrative]
  beta_r = 0.3,   # transmission from recurrence (per partnership) [illustrative]
  beta_l = 0.05,  # transmission from latency (asymptomatic shedding) (per partnership) [illustrative]
  d_p = 3 / 52,   # duration of primary episode (about 3 weeks) (yr) [illustrative]
  react = 4,      # reactivation rate (about 4 recurrences/yr in early years) (/yr) [illustrative]
  d_r = 5 / 365,  # duration of a recurrence (yr) [illustrative]
  mu = 1 / 30,    # exit from the sexually active population (/yr) [illustrative]
  tx = 0          # fraction of recurrences shortened by suppressive antivirals [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 0.8,    # susceptible
  Ip = 0.01,  # primary infection episode (high shedding)
  Lt = 0.17,  # latent (low-level asymptomatic shedding)
  Rc = 0.02   # recurrent episode (reactivated shedding)
)

# ---- Model equations -------------------------------------------------
hsv2_latent_reactivation_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + Ip + Lt + Rc

    # Derived quantities
    foi <- cp * (beta_p * Ip + beta_r * Rc + beta_l * Lt) / N  # force of infection

    # Rates of change
    dS <- mu * N - foi * S - mu * S
    dIp <- foi * S - 1 / d_p * Ip - mu * Ip
    dLt <- 1 / d_p * Ip - react * (1 - tx) * Lt + 1 / d_r * Rc - mu * Lt
    dRc <- react * (1 - tx) * Lt - 1 / d_r * Rc - mu * Rc

    list(c(dS, dIp, dLt, dRc))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = hsv2_latent_reactivation_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", Ip = "#B71C1C", Lt = "#E8C468", Rc = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Genital herpes (HSV-2) with primary episode, latency and recurrent reactivation")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
