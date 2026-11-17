## Virus dynamics with a CTL immune response
## Within-host virus–immune-response dynamics (formulated for HIV/HTLV-1 and LCMV).
## Four compartments, and the only nonlinearity beyond standard mass action is the bilinear CTL proliferation term cz Ii Zc, which in the generator's grammar is an absolute inflow to Zc (the rate is not per-capita in a single source compartmen
## Parameter values: 0 of 9 checked against the source; the rest are recalled or illustrative.
## Reference: Nowak MA, Bangham CRM (1996). Population dynamics of immune responses to persistent viruses. Science 272(5258):74–79
## Source: https://doi.org/10.1126/science.272.5258.74
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Tc, Ii, V, Zc
## Flows (time unit: days):
##   (outside) -> Tc  : lam
##   Tc -> (outside)  : dT * Tc
##   Tc -> Ii         : beta * V * Tc
##   Ii -> (outside)  : (aI + pz * Zc) * Ii
##   (outside) -> V   : kv * Ii
##   V -> (outside)   : uV * V
##   (outside) -> Zc  : cz * Ii * Zc
##   Zc -> (outside)  : bz * Zc

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 0.2)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  lam = 10,      # target-cell production [illustrative] (cells day⁻¹)
  dT = 0.1,      # target-cell death rate [illustrative] (per day)
  beta = 0.002,  # infection rate constant [illustrative] (per virus per day)
  aI = 0.5,      # infected-cell death rate in the absence of CTL [illustrative] (per day)
  pz = 1.0,      # rate of CTL-mediated killing of infected cells (pz because p is used for production and c is reserved) [illustrative] (per CTL per day)
  kv = 50,       # virus production per infected cell [illustrative] (virions cell⁻¹ day⁻¹)
  uV = 5,        # virus clearance rate [illustrative] (per day)
  cz = 0.1,      # rate of CTL proliferation driven by antigen [illustrative] (per infected cell per day)
  bz = 0.05      # CTL death rate [illustrative] (per day)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Tc = 1000,  # uninfected target cells
  Ii = 0,     # infected cells
  V = 1e-3,   # free virus
  Zc = 1      # virus-specific cytotoxic T lymphocytes (CTL)
)

# ---- Model equations -------------------------------------------------
hiv_ctl_nowak_bangham_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Derived quantities
    R0v <- beta * kv * (lam / dT) / (uV * aI)  # reproductive ratio without immune control
    Istar_imm <- bz / cz                       # infected-cell abundance at the immune-controlled equilibrium: a key published insight is that this level is set entirely by CTL parameters,
    Zstar <- (beta * Tc * V / Ii - aI) / pz    # equilibrium CTL level (evaluate symbolically, not in the ODE)

    # Rates of change
    dTc <- lam - dT * Tc - beta * V * Tc
    dIi <- beta * V * Tc - (aI + pz * Zc) * Ii
    dV <- kv * Ii - uV * V
    dZc <- cz * Ii * Zc - bz * Zc

    list(c(dTc, dIi, dV, dZc))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = hiv_ctl_nowak_bangham_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Tc = "#333333", Ii = "#B71C1C", V = "#2457C5", Zc = "#333333")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Virus dynamics with a CTL immune response")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
