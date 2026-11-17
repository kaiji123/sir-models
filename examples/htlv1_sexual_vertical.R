## HTLV-1 lifelong infection with sexual and mother-to-child transmission and disease (ATL/HAM)
## HTLV-1.
## Parameter values: 0 of 7 checked against the source; the rest are recalled or illustrative.
## Reference: I found no canonical population-level deterministic HTLV-1 transmission ODE paper in this session; published HTLV ODE models are mostly within-host. This is a composite model with illustrative values. The epidemiological facts behind it (about 2–5%…
## Source: https://doi.org/10.3389/fmicb.2012.00388
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sf, Hf, Sm, Hm, Dz
## Flows (time unit: years):
##   (outside) -> Sf  : 0.5 * births * (1 - pinf_birth)
##   (outside) -> Hf  : 0.5 * births * pinf_birth
##   (outside) -> Sm  : 0.5 * births * (1 - pinf_birth)
##   (outside) -> Hm  : 0.5 * births * pinf_birth
##   Sf -> Hf         : foi_f * Sf
##   Sm -> Hm         : foi_m * Sm
##   Hf -> Dz         : dis * Hf
##   Hm -> Dz         : dis * Hm
##   Dz -> (outside)  : (mu_d + mu) * Dz
##   Sf -> (outside)  : mu * Sf
##   Hf -> (outside)  : mu * Hf
##   Sm -> (outside)  : mu * Sm
##   Hm -> (outside)  : mu * Hm

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 0.2)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  cp = 1,        # partner change rate (partners/yr) [illustrative]
  b_mf = 0.06,   # male-to-female transmission (higher than reverse, per the literature) (per partnership/yr) [illustrative]
  b_fm = 0.01,   # female-to-male transmission (per partnership/yr) [illustrative]
  p_mtct = 0.2,  # infant infection probability from a carrier mother (prolonged breastfeeding) [illustrative]
  dis = 0.001,   # progression carrier to disease (about 4% lifetime over a 40-yr exposure) (/yr) [illustrative]
  mu = 1 / 50,   # natural death = birth (/yr) [illustrative]
  mu_d = 0.5     # disease-related death (ATL is rapidly fatal; HAM is not) (/yr) [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sf = 0.49,   # susceptible women
  Hf = 0.01,   # HTLV-1 carrier women
  Sm = 0.495,  # susceptible men
  Hm = 0.005,  # carrier men
  Dz = 0       # ATL/HAM disease (both sexes)
)

# ---- Model equations -------------------------------------------------
htlv1_sexual_vertical_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sf + Hf + Sm + Hm + Dz

    # Derived quantities
    Nf <- Sf + Hf + Dz / 2          # women total (assumes disease split evenly by sex)
    Nm <- Sm + Hm + Dz / 2          # men total
    foi_f <- cp * b_mf * Hm / Nm    # force of infection on women
    foi_m <- cp * b_fm * Hf / Nf    # force of infection on men
    births <- mu * N + mu_d * Dz    # births keeping N constant
    pinf_birth <- p_mtct * Hf / Nf  # fraction of births infected vertically

    # Rates of change
    dSf <- 0.5 * births * (1 - pinf_birth) - foi_f * Sf - mu * Sf
    dHf <- 0.5 * births * pinf_birth + foi_f * Sf - dis * Hf - mu * Hf
    dSm <- 0.5 * births * (1 - pinf_birth) - foi_m * Sm - mu * Sm
    dHm <- 0.5 * births * pinf_birth + foi_m * Sm - dis * Hm - mu * Hm
    dDz <- dis * Hf + dis * Hm - (mu_d + mu) * Dz

    list(c(dSf, dHf, dSm, dHm, dDz))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = htlv1_sexual_vertical_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sf = "#333333", Hf = "#6A2C9C", Sm = "#333333", Hm = "#6A2C9C", Dz = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "HTLV-1 lifelong infection with sexual and mother-to-child transmission and disease (ATL/HAM)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
