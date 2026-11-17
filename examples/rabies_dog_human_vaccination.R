## Dog rabies with mass vaccination and dog-to-human transmission
## Rabies, domestic dogs → humans, urban African setting (N'Djaména, Chad).
## Transmission is density-dependent with an explicit area term — beta_d has units km²·dog⁻¹·week⁻¹, so the per-capita hazard is beta_d Id / area.
## Parameter values: 0 of 15 checked against the source; the rest are recalled or illustrative.
## Reference: Zinsstag J, Dürr S, Penny MA, Mindekem R, Roth F, Menendez Gonzalez S, Naissengar S, Hattendorf J (2009). Transmission dynamics and economics of rabies control in dogs and humans in an African city. PNAS 106(35):14996–15001
## Source: https://doi.org/10.1073/pnas.0904740106
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sd, Ed, Id, Vd, Sh, Eh, Ih, Dh, Dd
## Flows (time unit: weeks):
##   (outside) -> Sd  : bd * Nd
##   Sd -> Ed         : lambda_d * Sd
##   Sd -> Vd         : nu_d * Sd
##   Vd -> Sd         : wane_d * Vd
##   Ed -> Id         : sigma_d * rd * Ed
##   Ed -> Sd         : sigma_d * (1 - rd) * Ed
##   Id -> Dd         : mu_d * Id
##   Sd -> (outside)  : (md + cull_d) * Sd
##   Sh -> Eh         : lambda_h * Sh
##   Eh -> Ih         : sigma_h * rh * (1 - pep) * Eh
##   Eh -> Sh         : sigma_h * (1 - rh * (1 - pep)) * Eh
##   Ih -> Dh         : mu_h * Ih

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 520, by = 1)  # weeks

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta_d = 0.0807,    # dog-to-dog transmission rate [verified] (estimated by the authors) (km² per dog per week)
  beta_dh = 0.0002,   # dog-to-human transmission rate [verified] (km² per dog per week)
  area = 100,         # area of the study city, converting counts to densities [illustrative] (km²)
  sigma_d = 1 / 3.1,  # rate of leaving dog incubation (≈22 d mean incubation) [illustrative] (per week)
  mu_d = 1 / 0.6,     # rabies-induced death rate of rabid dogs (≈4 d survival); the sensitivity analysis identifies this as one of the three most influential parameters [verified as i
  rd = 0.49,          # probability that an exposed dog develops clinical rabies; also identified as highly influential [verified as influential; value illustrative] (dimensionless)
  bd = 0.0052,        # dog birth rate [illustrative] (per week)
  md = 0.0048,        # dog non-rabies death rate [illustrative] (per week)
  nu_d = 0,           # dog vaccination rate (campaign lever) [illustrative] (per week)
  wane_d = 1 / 52,    # waning of dog vaccine immunity + turnover of vaccinated dogs [illustrative] (per week)
  cull_d = 0,         # dog culling rate (alternative control the authors simulate) [illustrative] (per week)
  rh = 0.19,          # probability a bitten human develops rabies absent prophylaxis [illustrative] (dimensionless)
  pep = 0,            # proportion of exposed humans receiving effective post-exposure prophylaxis [illustrative] (dimensionless)
  sigma_h = 1 / 6,    # rate of leaving human incubation (≈6 weeks) [illustrative] (per week)
  mu_h = 1            # death rate of clinical human rabies (uniformly fatal within days) [illustrative] (per week)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sd = 23000,  # susceptible dogs
  Ed = 10,     # incubating (exposed) dogs
  Id = 2,      # rabid, infectious dogs
  Vd = 0,      # vaccinated, immune dogs
  Sh = 6e5,    # susceptible humans
  Eh = 0,      # exposed humans (bitten, incubating)
  Ih = 0,      # humans with clinical rabies
  Dh = 0,      # cumulative human rabies deaths
  Dd = 0       # cumulative rabid-dog deaths (for fitting to weekly rabid-dog reports)
)

# ---- Model equations -------------------------------------------------
rabies_dog_human_vaccination_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N_total <- Sh + Eh + Ih

    # Derived quantities
    Nd <- Sd + Ed + Id + Vd            # dog population
    N <- Sh + Eh + Ih                  # human population
    lambda_d <- beta_d * Id / area     # dog force of infection (density-dependent: the transmission rate is per unit area)
    lambda_h <- beta_dh * Id / area    # human exposure hazard
    Re <- beta_d * Sd / (area * mu_d)  # effective reproduction ratio; the authors estimate Re ≈ 1.01 at the onset of observation, consistent with stable endemic dynamics [verified]
    cover <- Vd / Nd                   # vaccination coverage; a single campaign reaching ≥70% coverage was reported sufficient to interrupt transmission to humans for at least 6 ye

    # Rates of change
    dSd <- bd * Nd - lambda_d * Sd - nu_d * Sd + wane_d * Vd + sigma_d * (1 - rd) * Ed - (md + cull_d) * Sd
    dEd <- lambda_d * Sd - sigma_d * rd * Ed - sigma_d * (1 - rd) * Ed
    dId <- sigma_d * rd * Ed - mu_d * Id
    dVd <- nu_d * Sd - wane_d * Vd
    dSh <- -lambda_h * Sh + sigma_h * (1 - rh * (1 - pep)) * Eh
    dEh <- lambda_h * Sh - sigma_h * rh * (1 - pep) * Eh - sigma_h * (1 - rh * (1 - pep)) * Eh
    dIh <- sigma_h * rh * (1 - pep) * Eh - mu_h * Ih
    dDh <- mu_h * Ih
    dDd <- mu_d * Id

    list(c(dSd, dEd, dId, dVd, dSh, dEh, dIh, dDh, dDd))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = rabies_dog_human_vaccination_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sd = "#333333", Ed = "#E8C468", Id = "#B71C1C", Vd = "#2457C5", Sh = "#333333", Eh = "#E8C468", Ih = "#B71C1C", Dh = "#111111", Dd = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (weeks)", ylab = "Number of individuals",
        main = "Dog rabies with mass vaccination and dog-to-human transmission")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
