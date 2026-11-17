## HIV universal test-and-treat (4 untreated + 4 ART stages)
## HIV, generalised heterosexual epidemic (South Africa).
## Parameter values: 0 of 10 checked against the source; the rest are recalled or illustrative.
## Reference: Granich RM, Gilks CF, Dye C, De Cock KM, Williams BG. Universal voluntary HIV testing with immediate antiretroviral therapy as a strategy for elimination of HIV transmission: a mathematical model. Lancet 2009;373(9657):48–57
## Source: https://doi.org/10.1016/S0140-6736(08
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I1, I2, I3, I4, A1, A2, A3, A4
## Flows (time unit: years):
##   (outside) -> S   : (mu * N + rho_p * I4 + rho_a * A4)
##   S -> I1          : foi * S
##   I1 -> I2         : rho_p * I1
##   I2 -> I3         : rho_p * I2
##   I3 -> I4         : rho_p * I3
##   I4 -> (outside)  : rho_p * I4
##   I1 -> A1         : tau_t * I1
##   I2 -> A2         : tau_t * I2
##   I3 -> A3         : tau_t * I3
##   I4 -> A4         : tau_t * I4
##   A1 -> A2         : rho_a * A1
##   A2 -> A3         : rho_a * A2
##   A3 -> A4         : rho_a * A3
##   A4 -> (outside)  : rho_a * A4
##   A1 -> I1         : dropout * A1
##   A2 -> I2         : dropout * A2
##   A3 -> I3         : dropout * A3
##   A4 -> I4         : dropout * A4
##   S -> (outside)   : mu * S
##   I1 -> (outside)  : mu * I1
##   I2 -> (outside)  : mu * I2
##   I3 -> (outside)  : mu * I3
##   I4 -> (outside)  : mu * I4
##   A1 -> (outside)  : mu * A1
##   A2 -> (outside)  : mu * A2
##   A3 -> (outside)  : mu * A3
##   A4 -> (outside)  : mu * A4

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 50, by = 0.05)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  lambda0 = 0.6,    # transmission coefficient at low prevalence (/yr) [illustrative]
  alpha_h = 4,      # strength of the prevalence-dependent decline in transmission [illustrative]
  n_h = 1,          # shape exponent of that decline [illustrative]
  rho_p = 0.4,      # progression rate between untreated stages (4 stages of 2.5 yr give about 10 yr mean survival) (/yr) [illustrative]
  rho_a = 0.1,      # progression between stages while on ART (slower) (/yr) [illustrative]
  eps_art = 0.01,   # relative infectiousness on ART (99% reduction) [recalled]
  tau = 1,          # testing plus ART initiation rate; 1/yr = annual universal test-and-treat (/yr) [illustrative]
  dropout = 0.015,  # ART discontinuation (/yr) [illustrative]
  mu = 0.02,        # background adult mortality/exit (/yr) [illustrative]
  t_start = 10      # year the strategy starts (yr) [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 0.85,   # susceptible adults aged 15+
  I1 = 0.04,  # untreated HIV, stage 1
  I2 = 0.04,  # untreated HIV, stage 2
  I3 = 0.04,  # untreated HIV, stage 3
  I4 = 0.03,  # untreated HIV, stage 4 (AIDS)
  A1 = 0,     # on ART, entered from stage 1
  A2 = 0,     # on ART, entered from stage 2
  A3 = 0,     # on ART, entered from stage 3
  A4 = 0      # on ART, entered from stage 4
)

# ---- Model equations -------------------------------------------------
hiv_art_tasp_granich_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I1 + I2 + I3 + I4 + A1 + A2 + A3 + A4

    # Derived quantities
    P <- (I1 + I2 + I3 + I4 + A1 + A2 + A3 + A4) / N                          # HIV prevalence
    lam_eff <- lambda0 * exp(-alpha_h * P^n_h)                                # saturating transmission rate (heterogeneity proxy)
    foi <- lam_eff * (I1 + I2 + I3 + I4 + eps_art * (A1 + A2 + A3 + A4)) / N  # force of infection
    tau_t <- ifelse(t >= t_start, tau, 0)                                     # time-switched testing/ART rate

    # Rates of change
    dS <- (mu * N + rho_p * I4 + rho_a * A4) - foi * S - mu * S
    dI1 <- foi * S - rho_p * I1 - tau_t * I1 + dropout * A1 - mu * I1
    dI2 <- rho_p * I1 - rho_p * I2 - tau_t * I2 + dropout * A2 - mu * I2
    dI3 <- rho_p * I2 - rho_p * I3 - tau_t * I3 + dropout * A3 - mu * I3
    dI4 <- rho_p * I3 - rho_p * I4 - tau_t * I4 + dropout * A4 - mu * I4
    dA1 <- tau_t * I1 - rho_a * A1 - dropout * A1 - mu * A1
    dA2 <- tau_t * I2 + rho_a * A1 - rho_a * A2 - dropout * A2 - mu * A2
    dA3 <- tau_t * I3 + rho_a * A2 - rho_a * A3 - dropout * A3 - mu * A3
    dA4 <- tau_t * I4 + rho_a * A3 - rho_a * A4 - dropout * A4 - mu * A4

    list(c(dS, dI1, dI2, dI3, dI4, dA1, dA2, dA3, dA4))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = hiv_art_tasp_granich_model, parms = parameters, hmax = 0.05)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I1 = "#B71C1C", I2 = "#B71C1C", I3 = "#B71C1C", I4 = "#B71C1C", A1 = "#E57373", A2 = "#E57373", A3 = "#E57373", A4 = "#E57373")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "HIV universal test-and-treat (4 untreated + 4 ART stages)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
