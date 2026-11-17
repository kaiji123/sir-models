## Gonorrhoea SIS with core/non-core groups, two sexes
## Gonorrhoea, heterosexual.
## Parameter values: 0 of 7 checked against the source; the rest are recalled or illustrative.
## Reference: Hethcote HW, Yorke JA. Gonorrhea Transmission Dynamics and Control. Lecture Notes in Biomathematics vol. 56. Springer-Verlag, 1984
## Source: https://doi.org/10.1007/978-3-662-07544-9
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sm_c, Im_c, Sm_n, Im_n, Sf_c, If_c, Sf_n, If_n
## Flows (time unit: years):
##   Sm_c -> Im_c  : foi_m_c * Sm_c
##   Im_c -> Sm_c  : 1 / dur_m * Im_c
##   Sm_n -> Im_n  : foi_m_n * Sm_n
##   Im_n -> Sm_n  : 1 / dur_m * Im_n
##   Sf_c -> If_c  : foi_f_c * Sf_c
##   If_c -> Sf_c  : 1 / dur_f * If_c
##   Sf_n -> If_n  : foi_f_n * Sf_n
##   If_n -> Sf_n  : 1 / dur_f * If_n

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 20, by = 0.02)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  cp_c = 20,     # core partner rate (partners/yr) [illustrative]
  cp_n = 1,      # non-core partner rate (partners/yr) [illustrative]
  b_mf = 0.5,    # male-to-female transmission probability (per partnership) [illustrative]
  b_fm = 0.25,   # female-to-male transmission probability (per partnership) [illustrative]
  dur_m = 0.05,  # mean infectious duration in men (most symptomatic, treated in about 2–3 weeks) (yr) [illustrative]
  dur_f = 0.25,  # mean infectious duration in women (many asymptomatic) (yr) [illustrative]
  eps_mix = 0.6  # assortativity [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sm_c = 0.019,  # susceptible core men
  Im_c = 0.001,  # infected core men
  Sm_n = 0.98,   # susceptible non-core men
  Im_n = 0,      # infected non-core men
  Sf_c = 0.019,  # susceptible core women
  If_c = 0.001,  # infected core women
  Sf_n = 0.98,   # susceptible non-core women
  If_n = 0       # infected non-core women
)

# ---- Model equations -------------------------------------------------
gonorrhoea_core_group_hethcote_yorke_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sm_c + Im_c + Sm_n + Im_n + Sf_c + If_c + Sf_n + If_n

    # Derived quantities
    Nm_c <- Sm_c + Im_c                                                        # core men total
    Nm_n <- Sm_n + Im_n                                                        # non-core men total
    Nf_c <- Sf_c + If_c                                                        # core women total
    Nf_n <- Sf_n + If_n                                                        # non-core women total
    pf_mix <- (cp_c * If_c + cp_n * If_n) / (cp_c * Nf_c + cp_n * Nf_n)        # prevalence among female partnerships (proportionate mixing)
    pm_mix <- (cp_c * Im_c + cp_n * Im_n) / (cp_c * Nm_c + cp_n * Nm_n)        # prevalence among male partnerships (proportionate mixing)
    foi_m_c <- cp_c * b_fm * (eps_mix * If_c / Nf_c + (1 - eps_mix) * pf_mix)  # force of infection, core men
    foi_m_n <- cp_n * b_fm * (eps_mix * If_n / Nf_n + (1 - eps_mix) * pf_mix)  # force of infection, non-core men
    foi_f_c <- cp_c * b_mf * (eps_mix * Im_c / Nm_c + (1 - eps_mix) * pm_mix)  # force of infection, core women
    foi_f_n <- cp_n * b_mf * (eps_mix * Im_n / Nm_n + (1 - eps_mix) * pm_mix)  # force of infection, non-core women

    # Rates of change
    dSm_c <- -foi_m_c * Sm_c + 1 / dur_m * Im_c
    dIm_c <- foi_m_c * Sm_c - 1 / dur_m * Im_c
    dSm_n <- -foi_m_n * Sm_n + 1 / dur_m * Im_n
    dIm_n <- foi_m_n * Sm_n - 1 / dur_m * Im_n
    dSf_c <- -foi_f_c * Sf_c + 1 / dur_f * If_c
    dIf_c <- foi_f_c * Sf_c - 1 / dur_f * If_c
    dSf_n <- -foi_f_n * Sf_n + 1 / dur_f * If_n
    dIf_n <- foi_f_n * Sf_n - 1 / dur_f * If_n

    list(c(dSm_c, dIm_c, dSm_n, dIm_n, dSf_c, dIf_c, dSf_n, dIf_n))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = gonorrhoea_core_group_hethcote_yorke_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sm_c = "#333333", Im_c = "#B71C1C", Sm_n = "#333333", Im_n = "#B71C1C", Sf_c = "#333333", If_c = "#B71C1C", Sf_n = "#333333", If_n = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Gonorrhoea SIS with core/non-core groups, two sexes")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
