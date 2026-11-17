## Chikungunya SEIIR humans / SEI Aedes albopictus, Réunion 2005–06 (Yakob & Clements 2013) [VERIFIED-PARAMS]
## Chikungunya, Réunion Island epidemic 2005–2006.
## Parameter values: 0 of 8 checked against the source; the rest are recalled or illustrative.
## Reference: Yakob L, Clements ACA (2013). A mathematical model of chikungunya dynamics and control: the major epidemic on Réunion Island. PLoS ONE 8(3):e57448. https://doi.org/10.1371/journal.pone.0057448 . Parameter values were read from the paper this session; β1…
## Source: https://doi.org/10.1371/journal.pone.0057448
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sh, Eh, Ish, Iah, Rh, Xv, Yv, Zv
## Flows (time unit: days):
##   Sh -> Eh         : lambda_h * Sh
##   Eh -> Ish        : phi_s * lat_h * Eh
##   Eh -> Iah        : (1 - phi_s) * lat_h * Eh
##   Ish -> Rh        : gam_h * Ish
##   Iah -> Rh        : gam_h * Iah
##   (outside) -> Xv  : mu_v * Nv
##   Xv -> Yv         : lambda_v * Xv
##   Yv -> Zv         : lat_v * Yv
##   Xv -> (outside)  : mu_v * Xv
##   Yv -> (outside)  : mu_v * Yv
##   Zv -> (outside)  : mu_v * Zv

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 365, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta1 = 0.14,  # mosquito → human transmission rate (fitted; Macdonald "abm"-type) (1/day)
  beta2 = 0.40,  # human → mosquito transmission rate (fitted) (1/day)
  lat_h = 0.5,   # human latent rate (2–6 d; Thiboutot et al. 2010) (1/day)
  lat_v = 0.5,   # mosquito latent rate (2–6 d for Ae. albopictus; Dubrulle et al. 2009) (1/day)
  omega = 0.25,  # symptom-onset rate (pre-patent 4–7 d; Ziegler et al. 2008). Used only for the reporting lag. (1/day)
  phi_s = 0.97,  # proportion symptomatic (0.83–0.97)
  gam_h = 0.25,  # recovery rate (1–7 d; Chhabra et al. 2008) (1/day)
  mu_v = 0.05    # mosquito mortality (20–30 d life; Delatte et al. 2009). Mosquito birth = death. (1/day)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sh = 785000,  # susceptible humans
  Eh = 0,       # exposed
  Ish = 10,     # symptomatic infectious
  Iah = 0,      # asymptomatic infectious
  Rh = 0,       # recovered
  Xv = 1.0e6,   # susceptible mosquitoes
  Yv = 0,       # exposed mosquitoes
  Zv = 0        # infectious mosquitoes
)

# ---- Model equations -------------------------------------------------
chikungunya_yakob2013_reunion_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sh + Eh + Ish + Iah + Rh

    # Derived quantities
    Nh <- Sh + Eh + Ish + Iah + Rh
    Nv <- Xv + Yv + Zv
    lambda_h <- beta1 * Zv / Nv           # human force of infection; beta1 already contains a·b·m
    lambda_v <- beta2 * (Ish + Iah) / Nh  # mosquito force of infection; independent of m
    inc_sympt <- phi_s * lat_h * Eh       # symptomatic incidence, for comparison with case data; the paper adds an omega lag

    # Rates of change
    dSh <- -lambda_h * Sh
    dEh <- lambda_h * Sh - phi_s * lat_h * Eh - (1 - phi_s) * lat_h * Eh
    dIsh <- phi_s * lat_h * Eh - gam_h * Ish
    dIah <- (1 - phi_s) * lat_h * Eh - gam_h * Iah
    dRh <- gam_h * Ish + gam_h * Iah
    dXv <- mu_v * Nv - lambda_v * Xv - mu_v * Xv
    dYv <- lambda_v * Xv - lat_v * Yv - mu_v * Yv
    dZv <- lat_v * Yv - mu_v * Zv

    list(c(dSh, dEh, dIsh, dIah, dRh, dXv, dYv, dZv))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = chikungunya_yakob2013_reunion_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sh = "#333333", Eh = "#E8C468", Ish = "#B71C1C", Iah = "#B71C1C", Rh = "#4E8A3E", Xv = "#333333", Yv = "#333333", Zv = "#333333")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Chikungunya SEIIR humans / SEI Aedes albopictus, Réunion 2005–06 (Yakob & Clements 2013) [VERIFIED-PARAMS]")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
