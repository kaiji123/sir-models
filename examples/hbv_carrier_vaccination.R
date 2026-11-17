## Hepatitis B with acute infection, chronic carriers, vertical transmission and vaccination
## Hepatitis B.
## Parameter values: 0 of 11 checked against the source; the rest are recalled or illustrative.
## Reference: The conceptual basis is Medley GF, Lindop NA, Edmunds WJ, Nokes DJ. Hepatitis-B virus endemicity: heterogeneity, catastrophic dynamics and control. Nat Med 2001;7(5):619–624
## Source: https://doi.org/10.1038/87953
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, Ac, Cr, R, V
## Flows (time unit: years):
##   (outside) -> S   : (births_ok * (1 - vax) + mu_c * Cr)
##   (outside) -> V   : births_ok * vax
##   (outside) -> Cr  : mu * pv * Cr * pc0
##   (outside) -> Ac  : mu * pv * Cr * (1 - pc0)
##   S -> Ac          : foi * S
##   Ac -> Cr         : q_c * gam_a * Ac
##   Ac -> R          : (1 - q_c) * gam_a * Ac
##   Cr -> R          : gam_c * Cr
##   V -> S           : w_v * V
##   S -> (outside)   : mu * S
##   Ac -> (outside)  : mu * Ac
##   Cr -> (outside)  : (mu + mu_c) * Cr
##   R -> (outside)   : mu * R
##   V -> (outside)   : mu * V

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 100, by = 0.1)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta_h = 0.95,  # horizontal transmission coefficient (/yr) [recalled]
  eps_c = 0.16,   # relative infectiousness of carriers vs acute [recalled]
  gam_a = 4,      # exit from acute stage (about 3 months) (/yr) [recalled]
  q_c = 0.0886,   # probability that horizontal adult infection leads to carriage [recalled]
  gam_c = 0.025,  # carrier clearance (HBsAg loss) (/yr) [recalled]
  pv = 0.11,      # fraction of births to carrier mothers infected perinatally [recalled]
  pc0 = 0.885,    # probability that perinatal infection leads to carriage [recalled]
  mu = 0.00693,   # birth = natural death rate (/yr) [recalled]
  mu_c = 0.002,   # HBV-related excess carrier mortality (/yr) [illustrative]
  vax = 0.0,      # infant vaccination coverage times efficacy (scenario: 0.9) [illustrative]
  w_v = 0.1       # waning of vaccine immunity (/yr) [recalled]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 0.6,    # susceptible
  Ac = 0.01,  # acute infection
  Cr = 0.09,  # chronic carrier
  R = 0.3,    # recovered/immune
  V = 0       # vaccinated
)

# ---- Model equations -------------------------------------------------
hbv_carrier_vaccination_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + Ac + Cr + R + V

    # Derived quantities
    foi <- beta_h * (Ac + eps_c * Cr) / N  # force of infection
    births_ok <- mu * N - mu * pv * Cr     # births not perinatally infected

    # Rates of change
    dS <- (births_ok * (1 - vax) + mu_c * Cr) - foi * S + w_v * V - mu * S
    dAc <- mu * pv * Cr * (1 - pc0) + foi * S - q_c * gam_a * Ac - (1 - q_c) * gam_a * Ac - mu * Ac
    dCr <- mu * pv * Cr * pc0 + q_c * gam_a * Ac - gam_c * Cr - (mu + mu_c) * Cr
    dR <- (1 - q_c) * gam_a * Ac + gam_c * Cr - mu * R
    dV <- births_ok * vax - w_v * V - mu * V

    list(c(dS, dAc, dCr, dR, dV))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = hbv_carrier_vaccination_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", Ac = "#E57373", Cr = "#8E0E0E", R = "#4E8A3E", V = "#2457C5")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Hepatitis B with acute infection, chronic carriers, vertical transmission and vaccination")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
