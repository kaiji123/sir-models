## HIV-1 target-cell-limited model with antiretroviral therapy
## HIV-1 within-host viral dynamics.
## The canonical three-equation target-cell-limited system.
## Parameter values: 0 of 9 checked against the source; the rest are recalled or illustrative.
## Reference: Perelson AS, Neumann AU, Markowitz M, Leonard JM, Ho DD (1996). HIV-1 dynamics in vivo: virion clearance rate, infected cell life-span, and viral generation time. Science 271(5255):1582–1586
## Source: https://doi.org/10.1126/science.271.5255.1582
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Tc, Ii, V
## Flows (time unit: days):
##   (outside) -> Tc  : lam
##   Tc -> (outside)  : dT * Tc
##   Tc -> Ii         : beta * (1 - er) * V * Tc
##   Ii -> (outside)  : delta * Ii
##   (outside) -> V   : p * (1 - ep) * Ii
##   V -> (outside)   : cV * V

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 60, by = 0.1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  lam = 10,       # target-cell production rate [illustrative] (cells µL⁻¹ day⁻¹)
  dT = 0.01,      # natural death rate of uninfected target cells [illustrative] (per day)
  beta = 2.4e-5,  # infection rate constant [illustrative] (mL copies⁻¹ day⁻¹)
  delta = 0.5,    # death rate of productively infected cells; the authors' estimate implies an infected-cell half-life of ≈1.6 days [unverified; 0.5/day is the value conventionall
  p = 1000,       # virion production rate per infected cell [illustrative] (virions cell⁻¹ day⁻¹)
  cV = 3.07,      # virion clearance rate from plasma (c is reserved, hence cV) [unverified; 3.07/day is the value conventionally cited from this paper] (per day)
  eps_rt = 0,     # efficacy of a reverse-transcriptase inhibitor (blocks new infection) [illustrative] (dimensionless, 0–1)
  eps_pi = 0,     # efficacy of a protease inhibitor (blocks production of infectious virions) [illustrative] (dimensionless, 0–1)
  tau_rx = 30     # time therapy starts [illustrative] (days)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Tc = 1000,  # uninfected target cells (CD4+ T cells, per µL)
  Ii = 0,     # productively infected cells (per µL)
  V = 1e-3    # free virions (copies per mL)
)

# ---- Model equations -------------------------------------------------
hiv_target_cell_therapy_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Derived quantities
    on_rx <- ifelse(t < tau_rx, 0, 1)            # therapy switch
    er <- eps_rt * on_rx                         # effective RT-inhibitor efficacy at time t
    ep <- eps_pi * on_rx                         # effective PI efficacy at time t
    R0v <- beta * p * (lam / dT) / (cV * delta)  # basic reproductive ratio of the virus within host
    Tc_star <- cV * delta / (beta * p)           # infected steady-state target-cell level

    # Rates of change
    dTc <- lam - dT * Tc - beta * (1 - er) * V * Tc
    dIi <- beta * (1 - er) * V * Tc - delta * Ii
    dV <- p * (1 - ep) * Ii - cV * V

    list(c(dTc, dIi, dV))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = hiv_target_cell_therapy_model, parms = parameters, hmax = 0.1)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Tc = "#333333", Ii = "#B71C1C", V = "#2457C5")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "HIV-1 target-cell-limited model with antiretroviral therapy")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
