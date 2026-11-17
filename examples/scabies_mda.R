## Scabies SIS with egg-stage reinfection and two-dose ivermectin MDA
## Scabies (Sarcoptes scabiei).
## Parameter values: 0 of 9 checked against the source; the rest are recalled or illustrative.
## Reference: Lydeamore MJ, Campbell PT, Regan DG, Tong SYC, Andrews RM, Steer AC, Romani L, Kaldor JM, McVernon J, McCaw JM. A biological model of scabies infection dynamics and treatment informs mass drug administration strategies to increase the likelihood of…
## Source: https://doi.org/10.1016/j.mbs.2018.08.007
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, Ie
## Flows (time unit: days):
##   S -> I   : b_sc * I / N * S
##   I -> S   : 1 / dur_sc * I
##   I -> Ie  : mda_k * I
##   Ie -> I  : hatch * Ie
##   Ie -> S  : mda_k * Ie

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 3, by = 0.005)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  R0_sc = 1.5,     # R0 [I]
  dur_sc = 365,    # untreated infestation duration, chronic [I] (day)
  hatch = 1 / 7,   # egg hatch/maturation back to active infestation (ivermectin is not ovicidal) [M] (day^-1)
  cover_s = 0.85,  # MDA coverage [I]
  eff_s = 0.95,    # ivermectin kill of active mites [I]
  t1 = 100,        # first MDA dose [I] (day)
  t2 = 110,        # second dose, 7–14 days later [M on interval] (day)
  win_s = 1,       # pulse width [I] (day)
  N0 = 100000      # total initial population [illustrative default]
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  S <- 0.7 * N0  # uninfested
  I <- 0.3 * N0  # infested (active mites)
  Ie <- 0        # post-treatment carriage of surviving eggs, which re-hatch unless a second dose is given
  c(S = S, I = I, Ie = Ie)
})

# ---- Model equations -------------------------------------------------
scabies_mda_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + Ie

    # Derived quantities
    b_sc <- R0_sc / dur_sc
    dose_on <- ifelse((t >= t1 & t < t1 + win_s) | (t >= t2 & t < t2 + win_s), 1, 0)
    mda_k <- dose_on * (-log(1 - cover_s * eff_s) / win_s)

    # Rates of change
    dS <- -b_sc * I / N * S + 1 / dur_sc * I + mda_k * Ie
    dI <- b_sc * I / N * S - 1 / dur_sc * I - mda_k * I + hatch * Ie
    dIe <- mda_k * I - hatch * Ie - mda_k * Ie

    list(c(dS, dI, dIe))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = scabies_mda_model, parms = parameters, hmax = 0.005)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", Ie = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Scabies SIS with egg-stage reinfection and two-dose ivermectin MDA")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
