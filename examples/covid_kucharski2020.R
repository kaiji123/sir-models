## Early Wuhan SARS-CoV-2 SEIR with Erlang stages and reporting delay (deterministic version of Kucharski et al. 2020)
## COVID-19.
## Parameter values: 0 of 7 checked against the source; the rest are recalled or illustrative.
## Reference: Kucharski AJ, Russell TW, Diamond C, Liu Y, Edmunds J, Funk S, Eggo RM; CMMID COVID-19 working group. "Early dynamics of transmission and control of COVID-19: a mathematical modelling study." Lancet Infect Dis 2020;20:553–558.…
## Source: https://doi.org/10.1016/S1473-3099(20)30144-4
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E1, E2, I1, I2, Rm, Qr, Cr
## Flows (time unit: days):
##   S -> E1          : foi * S
##   E1 -> E2         : 2 / inc_per * E1
##   E2 -> I1         : 2 / inc_per * E2
##   I1 -> I2         : 2 / inf_per * I1
##   I2 -> Rm         : 2 / inf_per * I2
##   (outside) -> Qr  : p_rep * (2 / inc_per) * E2
##   Qr -> Cr         : 1 / rep_delay * Qr

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 120, by = 0.2)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  inc_per = 5.2,    # mean incubation (days)
  inf_per = 2.9,    # mean infectious period (days)
  rep_delay = 6.1,  # mean onset-to-confirmation delay (days)
  R_pre = 2.35,     # R(t) one week before and one week after the 23 Jan 2020 travel restrictions
  R_post = 1.05,    # R(t) one week before and one week after the 23 Jan 2020 travel restrictions
  t_switch = 23,    # date of the change (day) [illustrative]
  p_rep = 1         # fraction of cases eventually confirmed [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 11e6,  # susceptible
  E1 = 20,   # incubating, Erlang-2
  E2 = 0,    # incubating, Erlang-2
  I1 = 10,   # infectious, Erlang-2
  I2 = 0,    # infectious, Erlang-2
  Rm = 0,    # removed
  Qr = 0,    # symptomatic cases awaiting confirmation
  Cr = 0     # cumulative confirmed cases
)

# ---- Model equations -------------------------------------------------
covid_kucharski2020_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E1 + E2 + I1 + I2 + Rm

    # Derived quantities
    Rt <- ifelse(t >= t_switch, R_post, R_pre)
    foi <- Rt / inf_per * (I1 + I2) / N

    # Rates of change
    dS <- -foi * S
    dE1 <- foi * S - 2 / inc_per * E1
    dE2 <- 2 / inc_per * E1 - 2 / inc_per * E2
    dI1 <- 2 / inc_per * E2 - 2 / inf_per * I1
    dI2 <- 2 / inf_per * I1 - 2 / inf_per * I2
    dRm <- 2 / inf_per * I2
    dQr <- p_rep * (2 / inc_per) * E2 - 1 / rep_delay * Qr
    dCr <- 1 / rep_delay * Qr

    list(c(dS, dE1, dE2, dI1, dI2, dRm, dQr, dCr))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = covid_kucharski2020_model, parms = parameters, hmax = 0.2)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E1 = "#E8C468", E2 = "#E8C468", I1 = "#B71C1C", I2 = "#B71C1C", Rm = "#4E8A3E", Qr = "#00838F", Cr = "#8E0E0E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Early Wuhan SARS-CoV-2 SEIR with Erlang stages and reporting delay (deterministic version of Kucharski et al. 2020)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
