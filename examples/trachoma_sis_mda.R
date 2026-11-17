## Trachoma SIS with periodic mass azithromycin (Lietman 1999)
## Ocular Chlamydia trachomatis.
## Parameter values: 0 of 7 checked against the source; the rest are recalled or illustrative.
## Reference: Lietman T, Porco T, Dawson C, Blower S. Global elimination of trachoma: how frequently should we administer mass chemotherapy? Nat Med 1999;5(5):572–576
## Source: https://doi.org/10.1038/8451
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I
## Flows (time unit: weeks):
##   S -> I  : b_tr * I / N * S
##   I -> S  : (1 / dur_tr + mda_rate) * I

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 10, by = 0.01)  # weeks

# ---- Parameters ------------------------------------------------------
parameters <- c(
  dur_tr = 26,    # mean duration of infection in children [I] (week)
  R0_tr = 2,      # gives about 50% endemic prevalence (1 − 1/R0) [I]
  cover_m = 0.8,  # MDA coverage [I]
  eff_m = 0.95,   # azithromycin efficacy [M]
  int_m = 52,     # treatment interval: 52 annual, 26 biannual [V that these intervals were compared] (week)
  t_start = 52,   # first MDA round [I] (week)
  win_m = 1       # length of the treatment pulse [I] (week)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 500,  # uninfected
  I = 500   # infected with ocular chlamydia
)

# ---- Model equations -------------------------------------------------
trachoma_sis_mda_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I

    # Derived quantities
    b_tr <- R0_tr / dur_tr
    mda_on <- ifelse(t >= t_start & (t - t_start - int_m * floor((t - t_start) / int_m)) < win_m, 1, 0)
    mda_rate <- mda_on * (-log(1 - cover_m * eff_m) / win_m)  # a pulse giving the cover×eff kill in each round

    # Rates of change
    dS <- -b_tr * I / N * S + (1 / dur_tr + mda_rate) * I
    dI <- b_tr * I / N * S - (1 / dur_tr + mda_rate) * I

    list(c(dS, dI))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = trachoma_sis_mda_model, parms = parameters, hmax = 0.01)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (weeks)", ylab = "Number of individuals",
        main = "Trachoma SIS with periodic mass azithromycin (Lietman 1999)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
