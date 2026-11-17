## Meningococcal (serogroup C) carriage SCIS model with conjugate vaccine that protects against carriage (Trotter et al. 2005)
## Neisseria meningitidis C (applies to Hib with relabelling).
## Parameter values: 0 of 10 checked against the source; the rest are recalled or illustrative.
## Reference: Trotter CL, Gay NJ, Edmunds WJ. "Dynamic models of meningococcal carriage, disease, and the impact of serogroup C conjugate vaccination." Am J Epidemiol 2005;162:89–100. https://doi.org/10.1093/aje/kwi160 (record: https://researchonline.lshtm.ac.uk/6806/)
## Source: https://doi.org/10.1093/aje/kwi160
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, Ca, Rim, V, Cases
## Flows (time unit: years):
##   (outside) -> S      : (1 - p_vac) * mu * N
##   (outside) -> V      : p_vac * mu * N
##   S -> Ca             : foi * S
##   V -> Ca             : (1 - ve_c) * foi * V
##   Ca -> Rim           : clr * Ca
##   Rim -> S            : w_nat * Rim
##   V -> S              : w_vac * V
##   S -> V              : vac_rate * S
##   Rim -> V            : vac_rate * Rim
##   (outside) -> Cases  : ccr * acq
##   S -> (outside)      : mu * S

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 30, by = 0.05)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  R0_c = 1.36,   # Trotter 2005 abstract
  clr = 2,       # carriage clearance (about 6 months) (1/year) [illustrative]
  w_nat = 0.5,   # waning of immunity after carriage (1/year) [illustrative]
  ve_c = 0.75,   # vaccine efficacy against acquiring carriage [illustrative]
  w_vac = 0.1,   # waning of vaccine protection (1/year) [illustrative]
  vac_rate = 0,  # catch-up vaccination of S and Rim, e.g. 2 during a campaign (1/year) [illustrative]
  p_vac = 0.9,   # fraction of births vaccinated [illustrative]
  ccr = 1e-4,    # invasive cases per carriage acquisition (case:carrier ratio) [illustrative]
  mu = 0.0133,   # birth = death rate (1/year) [illustrative]
  N0 = 1e6       # value taken from the notes text [recalled]
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  S <- 0.6 * N0     # susceptible to carriage
  Ca <- 0.01 * N0   # carrier
  Rim <- 0.39 * N0  # immune after carriage
  V <- 0            # vaccinated, protected
  Cases <- 0        # cumulative invasive disease cases
  c(S = S, Ca = Ca, Rim = Rim, V = V, Cases = Cases)
})

# ---- Model equations -------------------------------------------------
meningo_carriage_trotter_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + Ca + Rim + V

    # Derived quantities
    foi <- R0_c * (clr + mu) * Ca / N
    acq <- foi * S + (1 - ve_c) * foi * V  # total acquisitions

    # Rates of change
    dS <- (1 - p_vac) * mu * N - foi * S + w_nat * Rim + w_vac * V - vac_rate * S - mu * S
    dCa <- foi * S + (1 - ve_c) * foi * V - clr * Ca
    dRim <- clr * Ca - w_nat * Rim - vac_rate * Rim
    dV <- p_vac * mu * N - (1 - ve_c) * foi * V - w_vac * V + vac_rate * S + vac_rate * Rim
    dCases <- ccr * acq

    list(c(dS, dCa, dRim, dV, dCases))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = meningo_carriage_trotter_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", Ca = "#8E0E0E", Rim = "#4E8A3E", V = "#2457C5", Cases = "#8E0E0E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Meningococcal (serogroup C) carriage SCIS model with conjugate vaccine that protects against carriage (Trotter et al. 2005)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
