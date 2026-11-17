## Two-serotype dengue with secondary infection (after Esteva & Vargas 2003) [ILLUSTRATIVE-COMPOSITE]
## Dengue (DENV, two co-circulating serotypes), tropical urban settings.
## Parameter values: 0 of 8 checked against the source; the rest are recalled or illustrative.
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sh, I1, I2, R1, R2, I12, I21, Rh, Sv, Iv1, Iv2
## Flows (time unit: days):
##   (outside) -> Sh   : mu_h * Nh
##   Sh -> I1          : lam1 * Sh
##   Sh -> I2          : lam2 * Sh
##   I1 -> R1          : gam_h * I1
##   I2 -> R2          : gam_h * I2
##   R1 -> I12         : x_sec * lam2 * R1
##   R2 -> I21         : x_sec * lam1 * R2
##   I12 -> Rh         : gam_h * I12
##   I21 -> Rh         : gam_h * I21
##   Sh -> (outside)   : mu_h * Sh
##   I1 -> (outside)   : mu_h * I1
##   I2 -> (outside)   : mu_h * I2
##   R1 -> (outside)   : mu_h * R1
##   R2 -> (outside)   : mu_h * R2
##   I12 -> (outside)  : mu_h * I12
##   I21 -> (outside)  : mu_h * I21
##   Rh -> (outside)   : mu_h * Rh
##   (outside) -> Sv   : A_v
##   Sv -> Iv1         : lamv1 * Sv
##   Sv -> Iv2         : lamv2 * Sv
##   Sv -> (outside)   : mu_v * Sv
##   Iv1 -> (outside)  : mu_v * Iv1
##   Iv2 -> (outside)  : mu_v * Iv2

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 3650, by = 5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  a_bite = 0.5,    # biting rate (1/day)
  p_vh = 0.5,      # mosquito → human per-bite transmission probability
  p_hv = 0.5,      # human → mosquito per-bite transmission probability
  gam_h = 0.1429,  # recovery rate (7-day viraemia) (1/day)
  mu_h = 4.2e-5,   # human birth/death rate (65-year life) (1/day)
  mu_v = 0.0714,   # mosquito death rate (14-day life) (1/day)
  A_v = 1428,      # mosquito recruitment, = mu_v × 20000 (mosquitoes/day)
  x_sec = 1.0      # susceptibility of R1/R2 to the heterologous serotype relative to Sh. Esteva & Vargas's enhancement factor; values >1 model ADE.
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sh = 9989,   # susceptible to both serotypes
  I1 = 10,     # primary infection with serotype 1
  I2 = 1,      # primary infection with serotype 2
  R1 = 0,      # recovered from serotype 1, susceptible to serotype 2
  R2 = 0,      # recovered from serotype 2, susceptible to serotype 1
  I12 = 0,     # secondary infection with serotype 2 after serotype 1
  I21 = 0,     # secondary infection with serotype 1 after serotype 2
  Rh = 0,      # immune to both
  Sv = 20000,  # susceptible mosquitoes
  Iv1 = 0,     # mosquitoes infectious with serotype 1
  Iv2 = 0      # mosquitoes infectious with serotype 2
)

# ---- Model equations -------------------------------------------------
dengue_two_serotype_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sh + I1 + I2 + R1 + R2 + I12 + I21 + Rh

    # Derived quantities
    Nh <- Sh + I1 + I2 + R1 + R2 + I12 + I21 + Rh
    Nv <- Sv + Iv1 + Iv2
    lam1 <- a_bite * p_vh * Iv1 / Nh          # human force of infection, serotype 1
    lam2 <- a_bite * p_vh * Iv2 / Nh          # human force of infection, serotype 2
    lamv1 <- a_bite * p_hv * (I1 + I21) / Nh  # mosquito force of infection, serotype 1
    lamv2 <- a_bite * p_hv * (I2 + I12) / Nh  # mosquito force of infection, serotype 2

    # Rates of change
    dSh <- mu_h * Nh - lam1 * Sh - lam2 * Sh - mu_h * Sh
    dI1 <- lam1 * Sh - gam_h * I1 - mu_h * I1
    dI2 <- lam2 * Sh - gam_h * I2 - mu_h * I2
    dR1 <- gam_h * I1 - x_sec * lam2 * R1 - mu_h * R1
    dR2 <- gam_h * I2 - x_sec * lam1 * R2 - mu_h * R2
    dI12 <- x_sec * lam2 * R1 - gam_h * I12 - mu_h * I12
    dI21 <- x_sec * lam1 * R2 - gam_h * I21 - mu_h * I21
    dRh <- gam_h * I12 + gam_h * I21 - mu_h * Rh
    dSv <- A_v - lamv1 * Sv - lamv2 * Sv - mu_v * Sv
    dIv1 <- lamv1 * Sv - mu_v * Iv1
    dIv2 <- lamv2 * Sv - mu_v * Iv2

    list(c(dSh, dI1, dI2, dR1, dR2, dI12, dI21, dRh, dSv, dIv1, dIv2))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = dengue_two_serotype_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sh = "#333333", I1 = "#B71C1C", I2 = "#B71C1C", R1 = "#4E8A3E", R2 = "#4E8A3E", I12 = "#B71C1C", I21 = "#B71C1C", Rh = "#4E8A3E", Sv = "#333333", Iv1 = "#B71C1C", Iv2 = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Two-serotype dengue with secondary infection (after Esteva & Vargas 2003) [ILLUSTRATIVE-COMPOSITE]")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
