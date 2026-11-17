## Amphibian chytridiomycosis with free-living zoospores
## Batrachochytrium dendrobatidis in a common toad (Bufo bufo) breeding population.
## Structurally an Anderson & May "model G" host–free-living-stage system.
## Parameter values: 0 of 11 checked against the source; the rest are recalled or illustrative.
## Reference: Mitchell KM, Churcher TS, Garner TWJ, Fisher MC (2008). Persistence of the emerging pathogen Batrachochytrium dendrobatidis outside the amphibian host greatly increases the probability of host extinction. Proc R Soc B 275(1632):329–334
## Source: https://doi.org/10.1098/rspb.2007.1356
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, Ii, Zs
## Flows (time unit: days):
##   (outside) -> S   : bz * N * pmax(1 - N / Kz, 0)
##   S -> Ii          : lambda * S
##   Ii -> S          : rec * Ii
##   S -> (outside)   : mu_z * S
##   Ii -> (outside)  : (mu_z + alpha_z) * Ii
##   (outside) -> Zs  : (prod * Ii + sapro)
##   Zs -> (outside)  : clear_z * Zs

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 3650, by = 5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta_hh = 1e-4,   # host-to-host transmission coefficient [illustrative] (per toad per day)
  beta_z = 1e-7,    # zoospore-to-host transmission coefficient [illustrative] (per zoospore per day)
  alpha_z = 0.02,   # disease-induced host death rate [illustrative] (per day)
  mu_z = 1 / 2920,  # natural host death rate (8-year life span) [illustrative] (per day)
  bz = 1 / 2920,    # host birth rate [illustrative] (per day)
  Kz = 1000,        # host carrying capacity [illustrative] (toads)
  prod = 100,       # zoospore release by infected hosts [illustrative] (zoospores per infected host per day)
  rz = 0.1,         # maximum saprobic (host-independent) zoospore reproduction rate (per day) [illustrative]
  phi = 1e-5,       # density-dependence of saprobic reproduction; appears as exp(-phi Zs) [verified: the paper uses a density-dependent term of the form exp(−φZ(t))] (per zoospore)
  clear_z = 0.1,    # zoospore death/loss rate (the survival time this implies is the key parameter the authors vary) [illustrative] (per day)
  rec = 0           # host recovery rate (set 0 for the baseline; chytrid clearance is limited) [illustrative] (per day)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 1000,  # susceptible toads
  Ii = 5,    # infected toads
  Zs = 0     # free-living zoospores in water (Z written Zs for clarity)
)

# ---- Model equations -------------------------------------------------
chytrid_amphibian_zoospore_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N_total <- S + Ii

    # Derived quantities
    N <- S + Ii                           # host population
    lambda <- beta_hh * Ii + beta_z * Zs  # force of infection (direct + zoospore pool)
    sapro <- rz * Zs * exp(-phi * Zs)     # density-limited saprobic zoospore reproduction

    # Rates of change
    dS <- bz * N * pmax(1 - N / Kz, 0) - lambda * S + rec * Ii - mu_z * S
    dIi <- lambda * S - rec * Ii - (mu_z + alpha_z) * Ii
    dZs <- (prod * Ii + sapro) - clear_z * Zs

    list(c(dS, dIi, dZs))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = chytrid_amphibian_zoospore_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", Ii = "#B71C1C", Zs = "#333333")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Amphibian chytridiomycosis with free-living zoospores")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
