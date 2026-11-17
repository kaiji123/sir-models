## Nipah virus bat → pig → human spillover with nosocomial transmission
## Nipah virus, Pteropus bats / pigs / humans.
## 13 compartments — at the structural limit, so the pig Rp and bat Rb classes can be dropped for a leaner variant.
## Parameter values: 0 of 16 checked against the source; the rest are recalled or illustrative.
## Reference: Structure after Pulliam JRC, Epstein JH, Dushoff J, Rahman SA, Bunning M, Jamaluddin AA, Hyatt AD, Field HE, Dobson AP, Daszak P, Henipavirus Ecology Research Group (2012). Agricultural intensification, priming for persistence and the emergence of Nipah…
## Source: https://doi.org/10.1098/rsif.2011.0223
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sb, Ib, Rb, Sp, Ep, Ip, Rp, Sh, Eh, Ih, Hh, Rh, Dh
## Flows (time unit: days):
##   (outside) -> Sb  : mu_b * Nb
##   Sb -> Ib         : lambda_b * Sb
##   Ib -> Rb         : gamma_b * Ib
##   Rb -> Sb         : wane_b * Rb
##   Sb -> (outside)  : mu_b * Sb
##   (outside) -> Sp  : mu_p * Np
##   Sp -> Ep         : lambda_p * Sp
##   Ep -> Ip         : sigma_p * Ep
##   Ip -> Rp         : gamma_p * Ip
##   Sp -> (outside)  : mu_p * Sp
##   Sh -> Eh         : lambda_h * Sh
##   Eh -> Ih         : sigma_h * Eh
##   Ih -> Hh         : eta_h * Ih
##   Ih -> Rh         : gamma_h * Ih
##   Hh -> Dh         : gamma_h * cfr_h * Hh
##   Hh -> Rh         : gamma_h * (1 - cfr_h) * Hh

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 730, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta_bb = 0.02,    # bat-to-bat transmission rate [illustrative] (per day)
  gamma_b = 1 / 14,  # bat recovery rate [illustrative] (per day)
  mu_b = 1 / 2920,   # bat birth = death rate (8-year lifespan) [illustrative] (per day)
  wane_b = 1 / 365,  # loss of bat immunity (drives periodic reintroduction) [illustrative] (per day)
  beta_bp = 5e-5,    # bat-to-pig spillover coefficient (contaminated fruit/urine route) [illustrative] (per bat per day)
  beta_pp = 1.2,     # pig-to-pig transmission rate (highly transmissible in intensive piggeries) [illustrative] (per day)
  sigma_p = 1 / 7,   # pig latency rate [illustrative] (per day)
  gamma_p = 1 / 10,  # pig recovery rate [illustrative] (per day)
  mu_p = 1 / 180,    # pig turnover (birth = offtake), drives susceptible replenishment [illustrative] (per day)
  beta_ph = 1e-6,    # pig-to-human transmission coefficient (occupational contact) [illustrative] (per pig per day)
  beta_hh = 0.02,    # community human-to-human transmission rate [illustrative] (per day)
  beta_hosp = 0.15,  # nosocomial/caregiver transmission rate from Hh [illustrative] (per day)
  sigma_h = 1 / 9,   # human latency rate (incubation ≈9 d) [illustrative] (per day)
  eta_h = 1 / 4,     # rate community case → hospital [illustrative] (per day)
  gamma_h = 1 / 14,  # human removal rate [illustrative] (per day)
  cfr_h = 0.7        # human case fatality (Bangladesh outbreaks ≈70%) [unverified] (dimensionless)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sb = 9000,   # susceptible bats
  Ib = 100,    # infectious bats
  Rb = 900,    # immune bats
  Sp = 10000,  # susceptible pigs
  Ep = 0,      # latent pigs
  Ip = 0,      # infectious pigs
  Rp = 0,      # recovered/immune pigs
  Sh = 50000,  # susceptible humans
  Eh = 0,      # latent humans
  Ih = 0,      # infectious humans in the community
  Hh = 0,      # hospitalised humans (source of nosocomial transmission)
  Rh = 0,      # recovered humans
  Dh = 0       # human deaths (cumulative)
)

# ---- Model equations -------------------------------------------------
nipah_bat_pig_human_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N_total <- Sh + Eh + Ih + Hh + Rh

    # Derived quantities
    Nb <- Sb + Ib + Rb                                                # bat population
    Np <- Sp + Ep + Ip + Rp                                           # pig population
    N <- Sh + Eh + Ih + Hh + Rh                                       # living human population
    lambda_b <- beta_bb * Ib / Nb                                     # bat force of infection
    lambda_p <- beta_pp * Ip / Np + beta_bp * Ib                      # pig force of infection (within-herd frequency-dependent plus spillover)
    lambda_h <- beta_ph * Ip + beta_hh * Ih / N + beta_hosp * Hh / N  # human force of infection

    # Rates of change
    dSb <- mu_b * Nb - lambda_b * Sb + wane_b * Rb - mu_b * Sb
    dIb <- lambda_b * Sb - gamma_b * Ib
    dRb <- gamma_b * Ib - wane_b * Rb
    dSp <- mu_p * Np - lambda_p * Sp - mu_p * Sp
    dEp <- lambda_p * Sp - sigma_p * Ep
    dIp <- sigma_p * Ep - gamma_p * Ip
    dRp <- gamma_p * Ip
    dSh <- -lambda_h * Sh
    dEh <- lambda_h * Sh - sigma_h * Eh
    dIh <- sigma_h * Eh - eta_h * Ih - gamma_h * Ih
    dHh <- eta_h * Ih - gamma_h * cfr_h * Hh - gamma_h * (1 - cfr_h) * Hh
    dRh <- gamma_h * Ih + gamma_h * (1 - cfr_h) * Hh
    dDh <- gamma_h * cfr_h * Hh

    list(c(dSb, dIb, dRb, dSp, dEp, dIp, dRp, dSh, dEh, dIh, dHh, dRh, dDh))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = nipah_bat_pig_human_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sb = "#333333", Ib = "#B71C1C", Rb = "#4E8A3E", Sp = "#333333", Ep = "#E8C468", Ip = "#B71C1C", Rp = "#4E8A3E", Sh = "#333333", Eh = "#E8C468", Ih = "#B71C1C", Hh = "#6A2C9C", Rh = "#4E8A3E", Dh = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Nipah virus bat → pig → human spillover with nosocomial transmission")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
