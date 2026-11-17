## Typhoid with chronic carriers, short- and long-cycle (water) transmission
## Typhoid fever (S. Typhi), South Asia (Kolkata, Delhi, Vellore).
## Parameter values: 3 of 11 checked against the source; the rest are recalled or illustrative.
## Reference: Pitzer VE, Bowles CC, Baker S, Kang G, Balaji V, Farrar JJ, Grenfell BT. Predicting the impact of vaccination on the transmission dynamics of typhoid in South Asia: a mathematical modeling study. PLoS Negl Trop Dis 2014;8(1):e2642
## Source: https://doi.org/10.1371/journal.pntd.0002642
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S1, I1, Car, R, S2, I2, W
## Flows (time unit: weeks):
##   (outside) -> S1   : mu_h * N
##   S1 -> I1          : lam * S1
##   I1 -> Car         : theta_c * delta_i * I1
##   I1 -> R           : (1 - theta_c) * delta_i * I1
##   R -> S2           : omega_r * R
##   S2 -> I2          : lam * S2
##   I2 -> R           : delta_i * I2
##   S1 -> (outside)   : mu_h * S1
##   I1 -> (outside)   : mu_h * I1
##   Car -> (outside)  : mu_h * Car
##   R -> (outside)    : mu_h * R
##   S2 -> (outside)   : mu_h * S2
##   I2 -> (outside)   : mu_h * I2
##   (outside) -> W    : xi_w * prev_w
##   W -> (outside)    : xi_w * W

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 50, by = 0.05)  # weeks

# ---- Parameters ------------------------------------------------------
parameters <- c(
  R0_ty = 2.8,       # basic reproduction number, Kolkata fit [verified]
  frac_short = 0.9,  # share of R0 from short-cycle transmission [verified]
  r_car = 0.01,      # relative infectiousness of carriers (best fit). The model also applies it to I2 [V for carriers; M/I for applying it to I2].
  theta_c = 0.03,    # probability that a primary infection becomes a lifelong carrier. The paper uses age-specific values of 0.3–10.1%, so 0.03 is an [I] average [V range].
  p_clin = 0.1,      # fraction of primary infections that are clinical typhoid. Used only for reporting [V].
  rep_f = 0.0052,    # reporting fraction [verified]
  delta_i = 0.25,    # recovery rate (4-week infectiousness) [M] (week^-1)
  omega_r = 0.0096,  # waning R → S2 (about 104 weeks) [M] (week^-1)
  xi_w = 0.333,      # decay of water contamination (about 3 weeks) [M] (week^-1)
  mu_h = 0.000385,   # birth = death rate (50-year life) [I]. The paper used observed demography. (week^-1)
  amp_w = 0.3        # seasonal amplitude of long-cycle transmission [I]. The paper fitted seasonality.
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S1 = 0.3,     # fully susceptible
  I1 = 100,     # primary (clinical-capable) infection
  Car = 0.005,  # chronic carrier. C is avoided because it is an R function
  R = 0.4,      # immune
  S2 = 0.3,     # partially susceptible after waning
  I2 = 0,       # subclinical reinfection
  W = 0         # water-supply contamination, scaled to prevalence units
)

# ---- Model equations -------------------------------------------------
typhoid_pitzer2014_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S1 + I1 + Car + R + S2 + I2

    # Derived quantities
    D_eff <- 1 / (delta_i + mu_h) + theta_c * delta_i / (delta_i + mu_h) * r_car / mu_h  # effective infectious duration per primary case (about 4.78 weeks)
    beta_s <- frac_short * R0_ty / D_eff                                                 # short-cycle rate (week^-1)
    beta_l <- (1 - frac_short) * R0_ty / D_eff                                           # long-cycle rate. W is scaled so that its equilibrium equals the weighted prevalence.
    prev_w <- (I1 + r_car * Car + r_car * I2) / N
    lam <- beta_s * prev_w + beta_l * (1 + amp_w * cos(6.2831853 * t / 52)) * W

    # Rates of change
    dS1 <- mu_h * N - lam * S1 - mu_h * S1
    dI1 <- lam * S1 - theta_c * delta_i * I1 - (1 - theta_c) * delta_i * I1 - mu_h * I1
    dCar <- theta_c * delta_i * I1 - mu_h * Car
    dR <- (1 - theta_c) * delta_i * I1 - omega_r * R + delta_i * I2 - mu_h * R
    dS2 <- omega_r * R - lam * S2 - mu_h * S2
    dI2 <- lam * S2 - delta_i * I2 - mu_h * I2
    dW <- xi_w * prev_w - xi_w * W

    list(c(dS1, dI1, dCar, dR, dS2, dI2, dW))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = typhoid_pitzer2014_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S1 = "#333333", I1 = "#B71C1C", Car = "#8E0E0E", R = "#4E8A3E", S2 = "#333333", I2 = "#B71C1C", W = "#A5D6A7")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (weeks)", ylab = "Number of individuals",
        main = "Typhoid with chronic carriers, short- and long-cycle (water) transmission")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
