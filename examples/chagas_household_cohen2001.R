## Domestic Chagas transmission: Triatoma infestans with humans, dogs and non-competent chickens (after Cohen & Gürtler 2001) [ILLUSTRATIVE-COMPOSITE]
## Chagas disease (Trypanosoma cruzi), rural households, Santiago del Estero, Argentina.
## Parameter values: 0 of 14 checked against the source; the rest are recalled or illustrative.
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Hs, Hi, Gs, Gi, Vs, Vi
## Flows (time unit: days):
##   (outside) -> Hs  : mu_hum * Nh
##   Hs -> Hi         : lam_hum * Hs
##   Hs -> (outside)  : mu_hum * Hs
##   Hi -> (outside)  : mu_hum * Hi
##   (outside) -> Gs  : mu_dog * Ng
##   Gs -> Gi         : lam_dog * Gs
##   Gs -> (outside)  : mu_dog * Gs
##   Gi -> (outside)  : mu_dog * Gi
##   (outside) -> Vs  : r_bug * Nv
##   Vs -> Vi         : lam_bug * Vs
##   Vs -> (outside)  : (mu_bug + (r_bug - mu_bug) * Nv / K_bug) * Vs
##   Vi -> (outside)  : (mu_bug + (r_bug - mu_bug) * Nv / K_bug) * Vi

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 7300, by = 10)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  a_b = 0.1429,      # blood meals per bug per day (weekly) (1/day)
  w_hum = 1,         # feeding weight on humans
  w_dog = 2.5,       # feeding weight on dogs (dogs preferred)
  w_chk = 1,         # feeding weight on chickens
  N_chk = 5,         # constant (chickens)
  p_vh = 5.8e-4,     # bug → human per-contact probability (Nouvellet 2013)
  p_vg = 0.001,      # bug → dog per-contact probability
  p_hv = 0.05,       # infected human → bug per-meal probability
  p_gv = 0.49,       # infected dog → bug per-meal probability
  mu_hum = 4.57e-5,  # human birth/death (1/day)
  mu_dog = 9.13e-4,  # dog birth/death (3-yr life) (1/day)
  r_bug = 0.02,      # bug per-capita recruitment (1/day)
  mu_bug = 0.005,    # bug death (1/day)
  K_bug = 1000       # (bugs/house)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Hs = 5,    # uninfected humans
  Hi = 0,    # infected humans (chronic, lifelong)
  Gs = 1.9,  # uninfected dogs
  Gi = 0.1,  # infected dogs
  Vs = 990,  # uninfected triatomine bugs
  Vi = 10    # infected bugs (lifelong)
)

# ---- Model equations -------------------------------------------------
chagas_household_cohen2001_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Hs + Hi + Gs + Gi

    # Derived quantities
    Nh <- Hs + Hi
    Ng <- Gs + Gi
    Nv <- Vs + Vi
    wsum <- w_hum * Nh + w_dog * Ng + w_chk * N_chk
    f_hum <- w_hum * Nh / wsum                                          # share of bites on humans
    f_dog <- w_dog * Ng / wsum                                          # share of bites on dogs
    lam_hum <- a_b * f_hum * p_vh * Vi / Nh                             # human force of infection
    lam_dog <- a_b * f_dog * p_vg * Vi / Ng                             # dog force of infection
    lam_bug <- a_b * (f_hum * p_hv * Hi / Nh + f_dog * p_gv * Gi / Ng)  # bug force of infection

    # Rates of change
    dHs <- mu_hum * Nh - lam_hum * Hs - mu_hum * Hs
    dHi <- lam_hum * Hs - mu_hum * Hi
    dGs <- mu_dog * Ng - lam_dog * Gs - mu_dog * Gs
    dGi <- lam_dog * Gs - mu_dog * Gi
    dVs <- r_bug * Nv - lam_bug * Vs - (mu_bug + (r_bug - mu_bug) * Nv / K_bug) * Vs
    dVi <- lam_bug * Vs - (mu_bug + (r_bug - mu_bug) * Nv / K_bug) * Vi

    list(c(dHs, dHi, dGs, dGi, dVs, dVi))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = chagas_household_cohen2001_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Hs = "#6A2C9C", Hi = "#6A2C9C", Gs = "#333333", Gi = "#333333", Vs = "#2457C5", Vi = "#2457C5")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Domestic Chagas transmission: Triatoma infestans with humans, dogs and non-competent chickens (after Cohen & Gürtler 2001) [ILLUSTRATIVE-COMPOSITE]")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
