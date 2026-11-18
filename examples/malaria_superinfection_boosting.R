## Ross–Macdonald with superinfection clearance and boosted acquired immunity (Macdonald/Dietz superinfection + Aron 1988 boosting) [ILLUSTRATIVE-COMPOSITE]
## P. falciparum malaria, stable endemic areas.
## Parameter values: 0 of 9 checked against the source; the rest are recalled or illustrative.
## Reference: Smith DL, Battle KE, Hay SI, et al. (2012) Ross, Macdonald, and a theory for the dynamics and control of mosquito-transmitted pathogens. PLoS Pathog 8(4):e1002588
## Source: https://doi.org/10.1371/journal.ppat.1002588
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sh, Ih, Rh, Sv, Ev, Iv
## Flows (time unit: days):
##   (outside) -> Sh  : mu_h * Nh
##   Sh -> Ih         : h_foi * Sh
##   Ih -> Rh         : r_eff * Ih
##   Rh -> Sh         : w_loss * Rh
##   Sh -> (outside)  : mu_h * Sh
##   Ih -> (outside)  : mu_h * Ih
##   Rh -> (outside)  : mu_h * Rh
##   (outside) -> Sv  : mu_v * Nv
##   Sv -> Ev         : lambda_v * Sv
##   Ev -> Iv         : nu_v * Ev
##   Sv -> (outside)  : mu_v * Sv
##   Ev -> (outside)  : mu_v * Ev
##   Iv -> (outside)  : mu_v * Iv

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 7300, by = 10)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  a_bite = 0.3,     # human-biting rate per mosquito (1/day)
  p_mh = 0.1,       # mosquito → human transmission probability per bite
  p_hm = 0.5,       # human → mosquito transmission probability per bite
  r_clear = 0.005,  # clearance rate of a single infection (about 200 days) (1/day)
  tau_imm = 1825,   # immune memory without boosting (illustrative) (days)
  mu_h = 5.5e-5,    # human birth = death rate (50-year life) (1/day)
  mu_v = 0.1,       # mosquito death rate (1/day)
  nu_v = 0.1,       # 1/EIP (1/day)
  m_ratio = 5       # sets the initial Nv; mosquito births balance deaths (mosquitoes/human)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sh = 990,   # susceptible humans
  Ih = 10,    # infected/infectious humans (possibly carrying several broods)
  Rh = 0,     # immune humans: not infectious, immunity kept up by re-exposure
  Sv = 4950,  # susceptible mosquitoes
  Ev = 0,     # mosquitoes in extrinsic incubation
  Iv = 50     # infectious mosquitoes
)

# ---- Model equations -------------------------------------------------
malaria_superinfection_boosting_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sh + Ih + Rh

    # Derived quantities
    Nh <- Sh + Ih + Rh
    Nv <- Sv + Ev + Iv
    h_foi <- pmax(a_bite * p_mh * Iv / Nh, 1e-10)  # human force of infection (EIR × b); the floor avoids 0/0
    r_eff <- h_foi / (exp(h_foi / r_clear) - 1)    # superinfection-adjusted clearance rate (tends to r_clear as h_foi → 0)
    w_loss <- h_foi / (exp(h_foi * tau_imm) - 1)   # immunity loss rate with boosting (tends to 1/tau_imm as h_foi → 0)
    lambda_v <- a_bite * p_hm * Ih / Nh            # mosquito force of infection

    # Rates of change
    dSh <- mu_h * Nh - h_foi * Sh + w_loss * Rh - mu_h * Sh
    dIh <- h_foi * Sh - r_eff * Ih - mu_h * Ih
    dRh <- r_eff * Ih - w_loss * Rh - mu_h * Rh
    dSv <- mu_v * Nv - lambda_v * Sv - mu_v * Sv
    dEv <- lambda_v * Sv - nu_v * Ev - mu_v * Ev
    dIv <- nu_v * Ev - mu_v * Iv

    list(c(dSh, dIh, dRh, dSv, dEv, dIv))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = malaria_superinfection_boosting_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sh = "#333333", Ih = "#B71C1C", Rh = "#4E8A3E", Sv = "#333333", Ev = "#E8C468", Iv = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Ross–Macdonald with superinfection clearance and boosted acquired immunity (Macdonald/Dietz superinfection + Aron 1988 boosting) [ILLUSTRATIVE-COMPOSITE]")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
