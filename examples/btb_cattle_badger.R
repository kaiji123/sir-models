## Bovine tuberculosis in cattle and badgers with test-and-slaughter
## Mycobacterium bovis, cattle ↔ European badger (Meles meles).
## The two-host coupling is the whole point: the paper examines how much of cattle infection is attributable to badgers versus within-herd spread, and how badger culling changes that.
## Parameter values: 0 of 12 checked against the source; the rest are recalled or illustrative.
## Reference: Cox DR, Donnelly CA, Bourne FJ, Gettinby G, McInerney JP, Morrison WI, Woodroffe R (2005). Simple model for tuberculosis in cattle and badgers. PNAS 102(49):17588–17593
## Source: https://doi.org/10.1073/pnas.0509003102
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sc, Ec, Ic, Sb, Eb, Ib, CRc
## Flows (time unit: years):
##   (outside) -> Sc  : (mu_c * N + test * Ic)
##   Sc -> Ec         : lambda_c * Sc
##   Ec -> Ic         : sigma_c * Ec
##   Ic -> CRc        : test * Ic
##   Sc -> (outside)  : mu_c * Sc
##   (outside) -> Sb  : (mu_b * Nb + alpha_b * Ib + cull_b * Nb)
##   Sb -> Eb         : lambda_b * Sb
##   Eb -> Ib         : sigma_b * Eb
##   Sb -> (outside)  : (mu_b + cull_b) * Sb

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 40, by = 0.05)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta_cc = 0.4,    # cattle-to-cattle transmission rate [illustrative] (per year)
  beta_bc = 0.02,   # badger-to-cattle transmission coefficient [illustrative] (per badger per year)
  beta_bb = 0.6,    # badger-to-badger transmission rate [illustrative] (per year)
  beta_cb = 0.002,  # cattle-to-badger transmission coefficient [illustrative] (per cattle per year)
  sigma_c = 2.0,    # rate cattle progress from test-negative to test-positive [illustrative] (per year)
  test = 1.0,       # herd-testing frequency × test sensitivity, i.e. rate at which Ic animals are detected and slaughtered [illustrative] (per year)
  mu_c = 0.25,      # cattle turnover (culling/sale, 4-year productive life) [illustrative] (per year)
  sigma_b = 1.5,    # badger progression rate to excreting state [illustrative] (per year)
  mu_b = 0.3,       # badger natural death rate [illustrative] (per year)
  alpha_b = 0.2,    # additional disease-induced badger death rate [illustrative] (per year)
  cull_b = 0,       # badger culling rate (the intervention the trial was designed to evaluate) [illustrative] (per year)
  perturb = 0       # proportional increase in beta_bc caused by culling-induced badger ranging ("perturbation effect"); e.g. 0.3 [illustrative] (dimensionless)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sc = 10000,  # susceptible cattle
  Ec = 0,      # cattle infected but test-negative (occult/early infection)
  Ic = 10,     # cattle infected and test-positive / infectious
  Sb = 1000,   # susceptible badgers
  Eb = 0,      # badgers infected, not yet excreting
  Ib = 20,     # badgers excreting M. bovis
  CRc = 0      # cumulative cattle reactors removed (herd-breakdown counter)
)

# ---- Model equations -------------------------------------------------
btb_cattle_badger_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N_total <- Sc + Ec + Ic

    # Derived quantities
    N <- Sc + Ec + Ic                                                           # cattle herd size
    Nb <- Sb + Eb + Ib                                                          # badger population
    lambda_c <- beta_cc * Ic / N + beta_bc * (1 + perturb * (cull_b > 0)) * Ib  # cattle force of infection
    lambda_b <- beta_bb * Ib / Nb + beta_cb * Ic                                # badger force of infection

    # Rates of change
    dSc <- (mu_c * N + test * Ic) - lambda_c * Sc - mu_c * Sc
    dEc <- lambda_c * Sc - sigma_c * Ec
    dIc <- sigma_c * Ec - test * Ic
    dSb <- (mu_b * Nb + alpha_b * Ib + cull_b * Nb) - lambda_b * Sb - (mu_b + cull_b) * Sb
    dEb <- lambda_b * Sb - sigma_b * Eb
    dIb <- sigma_b * Eb
    dCRc <- test * Ic

    list(c(dSc, dEc, dIc, dSb, dEb, dIb, dCRc))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = btb_cattle_badger_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sc = "#333333", Ec = "#E8C468", Ic = "#B71C1C", Sb = "#333333", Eb = "#E8C468", Ib = "#B71C1C", CRc = "#8E0E0E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Bovine tuberculosis in cattle and badgers with test-and-slaughter")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
