## Influenza A target-cell-limited model with an eclipse phase (TEIV)
## Influenza A within-host dynamics in experimentally infected humans.
## Four compartments, no target-cell regeneration and no host death terms — over a 7–10 day infection, target-cell depletion alone terminates the infection, which is the model's central claim.
## Parameter values: 0 of 6 checked against the source; the rest are recalled or illustrative.
## Reference: Baccam P, Beauchemin C, Macken CA, Hayden FG, Perelson AS (2006). Kinetics of influenza A virus infection in humans. Journal of Virology 80(15):7590–7599
## Source: https://doi.org/10.1128/JVI.01623-05
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Tc, I1, I2, V
## Flows (time unit: days):
##   Tc -> I1         : beta * V * Tc
##   I1 -> I2         : k * I1
##   I2 -> (outside)  : delta * I2
##   (outside) -> V   : p * I2
##   V -> (outside)   : cV * V

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 10, by = 0.01)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 2.7e-5,  # infection rate constant [unverified (mL TCID50⁻¹ day⁻¹)
  k = 4.0,        # rate of transition from eclipse to productively infected [unverified] (per day)
  delta = 2.4,    # death rate of productively infected cells (a secondary source reports the original estimate as 0.167 per hour ≈ 4.0 per day and notes it implies an implausibly
  p = 1.2e-2,     # virion production rate per productively infected cell [unverified; secondary source reports ≈9.9e4 and ≈5e4 for differently scaled versions] (TCID50 mL⁻¹ day⁻¹
  cV = 3.0,       # viral clearance rate (secondary source: 0.125 per hour = 3.0 per day) [unverified] (per day)
  Tc0 = 4e8       # initial value of Tc
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Tc = 4e8,   # uninfected target (epithelial) cells
  I1 = 0,     # infected cells in the eclipse phase (not yet producing virus)
  I2 = 0,     # productively infected cells
  V = 7.5e-2  # infectious virus titre (TCID50 per mL of nasal wash)
)

# ---- Model equations -------------------------------------------------
influenza_within_host_teiv_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Derived quantities
    R0v <- beta * p * Tc0 / (cV * delta)  # within-host basic reproductive ratio, with Tc0 the initial target-cell count (this expression for R0 is confirmed in the derived literature)
    logV <- log(pmax(V, 1e-12))           # log10-scale output is what the data are fitted on (divide by log(10))

    # Rates of change
    dTc <- -beta * V * Tc
    dI1 <- beta * V * Tc - k * I1
    dI2 <- k * I1 - delta * I2
    dV <- p * I2 - cV * V

    list(c(dTc, dI1, dI2, dV))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = influenza_within_host_teiv_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Tc = "#333333", I1 = "#B71C1C", I2 = "#B71C1C", V = "#2457C5")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Influenza A target-cell-limited model with an eclipse phase (TEIV)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
