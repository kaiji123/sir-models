## Hepatitis C within-host dynamics under interferon-α therapy
## Hepatitis C virus within-host dynamics.
## Structurally identical to the HIV target-cell model — the scientific content is in how therapy enters.
## Parameter values: 0 of 10 checked against the source; the rest are recalled or illustrative.
## Reference: Neumann AU, Lam NP, Dahari H, Gretch DR, Wiley TE, Layden TJ, Perelson AS (1998). Hepatitis C viral dynamics in vivo and the antiviral efficacy of interferon-α therapy. Science 282(5386):103–107
## Source: https://doi.org/10.1126/science.282.5386.103
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Tc, Ii, V
## Flows (time unit: days):
##   (outside) -> Tc  : sh
##   Tc -> (outside)  : dh * Tc
##   Tc -> Ii         : beta * (1 - ei) * V * Tc
##   Ii -> (outside)  : delta * Ii
##   (outside) -> V   : p * (1 - en) * Ii
##   V -> (outside)   : cV * V

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 28, by = 0.05)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  sh = 1e5,      # hepatocyte production rate [illustrative] (cells mL⁻¹ day⁻¹)
  dh = 0.01,     # death rate of uninfected hepatocytes [illustrative] (per day)
  beta = 1e-7,   # infection rate constant [illustrative] (mL virion⁻¹ day⁻¹)
  delta = 0.14,  # death/loss rate of infected hepatocytes; the authors estimated a wide range across patients (≈0.01–0.5 per day), corresponding to infected-cell half-lives of ≈1
  p = 100,       # virion production rate per infected cell [illustrative] (virions cell⁻¹ day⁻¹)
  cV = 6.2,      # virion clearance rate, implying a plasma virion half-life of ≈2.7 hours [unverified; 6.2/day is the headline estimate conventionally attributed to this paper] (
  eta = 0.95,    # interferon efficacy in blocking virion production (dimensionless, 0–1) [recalled]
  eps = 0,       # efficacy in blocking new infections (found to be the minor mechanism) [illustrative] (dimensionless, 0–1)
  tau_rx = 0,    # therapy start time [illustrative] (days)
  tau_lag = 0.3  # pharmacological delay before the drug acts (the authors needed a delay to fit the first hours) [illustrative] (days)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Tc = 1e7,  # uninfected target hepatocytes (cells per mL)
  Ii = 1e6,  # productively infected hepatocytes (cells per mL)
  V = 1e6    # free virions (IU or copies per mL)
)

# ---- Model equations -------------------------------------------------
hcv_neumann_interferon_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Tc + Ii + V

    # Derived quantities
    on_rx <- ifelse(t < tau_rx + tau_lag, 0, 1)  # therapy switch including the delay
    en <- eta * on_rx                            # effective production-blocking efficacy
    ei <- eps * on_rx                            # effective infection-blocking efficacy
    R0v <- beta * p * (sh / dh) / (cV * delta)   # within-host reproductive ratio

    # Rates of change
    dTc <- sh - dh * Tc - beta * (1 - ei) * V * Tc
    dIi <- beta * (1 - ei) * V * Tc - delta * Ii
    dV <- p * (1 - en) * Ii - cV * V

    list(c(dTc, dIi, dV))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = hcv_neumann_interferon_model, parms = parameters, hmax = 0.05)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Tc = "#333333", Ii = "#B71C1C", V = "#2457C5")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Hepatitis C within-host dynamics under interferon-α therapy")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
