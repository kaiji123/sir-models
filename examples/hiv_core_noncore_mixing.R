## HIV in core and non-core risk groups with preferential mixing and ART
## HIV, heterogeneous sexual activity.
## Parameter values: 0 of 10 checked against the source; the rest are recalled or illustrative.
## Reference: The core-group concept and the mixing structure come from Hethcote HW, Yorke JA. Gonorrhea Transmission Dynamics and Control. Lecture Notes in Biomathematics 56, Springer, 1984
## Source: https://doi.org/10.1007/978-3-662-07544-9
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S_h, I_h, Y_h, A_h, S_l, I_l, Y_l, A_l
## Flows (time unit: years):
##   (outside) -> S_h  : f_h * (mu * N + delta_y * (Y_h + Y_l))
##   (outside) -> S_l  : (1 - f_h) * (mu * N + delta_y * (Y_h + Y_l))
##   S_h -> I_h        : foi_h * S_h
##   S_l -> I_l        : foi_l * S_l
##   I_h -> Y_h        : nu_i * I_h
##   I_l -> Y_l        : nu_i * I_l
##   I_h -> A_h        : tau * I_h
##   I_l -> A_l        : tau * I_l
##   Y_h -> (outside)  : delta_y * Y_h
##   Y_l -> (outside)  : delta_y * Y_l
##   S_h -> (outside)  : mu * S_h
##   I_h -> (outside)  : mu * I_h
##   Y_h -> (outside)  : mu * Y_h
##   A_h -> (outside)  : mu * A_h
##   S_l -> (outside)  : mu * S_l
##   I_l -> (outside)  : mu * I_l
##   Y_l -> (outside)  : mu * Y_l
##   A_l -> (outside)  : mu * A_l

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 60, by = 0.1)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  cp_h = 10,       # partner change rate, core (partners/yr) [illustrative]
  cp_l = 1,        # partner change rate, non-core (partners/yr) [illustrative]
  beta_p = 0.05,   # transmission probability per partnership (per partnership) [illustrative]
  eps_mix = 0.5,   # assortativity: 0 = proportionate mixing, 1 = fully assortative [illustrative]
  nu_i = 0.1,      # progression from chronic infection to AIDS (/yr) [illustrative]
  delta_y = 1,     # AIDS mortality (/yr) [illustrative]
  tau = 0.2,       # ART initiation (/yr) [illustrative]
  eps_art = 0.04,  # relative infectiousness on ART (96% reduction; HPTN 052, Cohen et al. NEJM 2011, DOI 10.1056/NEJMoa1105243) [recalled]
  mu = 0.03,       # exit from the sexually active population (/yr) [illustrative]
  f_h = 0.05       # fraction of recruits entering the core [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S_h = 0.049,  # susceptible, high activity (core)
  I_h = 0.001,  # infected untreated, high activity
  Y_h = 0,      # AIDS, high activity
  A_h = 0,      # on ART, high activity
  S_l = 0.95,   # susceptible, low activity
  I_l = 0,      # infected untreated, low activity
  Y_l = 0,      # AIDS, low activity
  A_l = 0       # on ART, low activity
)

# ---- Model equations -------------------------------------------------
hiv_core_noncore_mixing_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S_h + I_h + Y_h + A_h + S_l + I_l + Y_l + A_l

    # Derived quantities
    N_h <- S_h + I_h + Y_h + A_h                                                          # core group total
    N_l <- S_l + I_l + Y_l + A_l                                                          # non-core group total
    Ninf_h <- I_h + eps_art * A_h                                                         # effective infectious count, core
    Ninf_l <- I_l + eps_art * A_l                                                         # effective infectious count, non-core
    act_h <- cp_h * (N_h - Y_h)                                                           # partnerships offered by the core (AIDS cases assumed inactive)
    act_l <- cp_l * (N_l - Y_l)                                                           # partnerships offered by the non-core
    prop_inf <- (cp_h * Ninf_h + cp_l * Ninf_l) / (act_h + act_l)                         # infectious fraction among partnerships under proportionate mixing
    foi_h <- cp_h * beta_p * (eps_mix * Ninf_h / (N_h - Y_h) + (1 - eps_mix) * prop_inf)  # force of infection, core
    foi_l <- cp_l * beta_p * (eps_mix * Ninf_l / (N_l - Y_l) + (1 - eps_mix) * prop_inf)  # force of infection, non-core

    # Rates of change
    dS_h <- f_h * (mu * N + delta_y * (Y_h + Y_l)) - foi_h * S_h - mu * S_h
    dI_h <- foi_h * S_h - nu_i * I_h - tau * I_h - mu * I_h
    dY_h <- nu_i * I_h - delta_y * Y_h - mu * Y_h
    dA_h <- tau * I_h - mu * A_h
    dS_l <- (1 - f_h) * (mu * N + delta_y * (Y_h + Y_l)) - foi_l * S_l - mu * S_l
    dI_l <- foi_l * S_l - nu_i * I_l - tau * I_l - mu * I_l
    dY_l <- nu_i * I_l - delta_y * Y_l - mu * Y_l
    dA_l <- tau * I_l - mu * A_l

    list(c(dS_h, dI_h, dY_h, dA_h, dS_l, dI_l, dY_l, dA_l))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = hiv_core_noncore_mixing_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S_h = "#333333", I_h = "#B71C1C", Y_h = "#333333", A_h = "#E57373", S_l = "#333333", I_l = "#B71C1C", Y_l = "#333333", A_l = "#E57373")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "HIV in core and non-core risk groups with preferential mixing and ART")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
