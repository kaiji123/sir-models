## sensitive and resistant strains with treatment and de novo resistance (single drug)
## bacterial infection with antibiotic treatment.
## Resistance takes over when f·h > Δr, as the ρ equation shows.
## Parameter values: 0 of 9 checked against the source; the rest are recalled or illustrative.
## Reference: Bonhoeffer S., Lipsitch M., Levin B.R. (1997) "Evaluating treatment protocols to prevent antibiotic resistance", PNAS 94(22):12106–12111
## Source: https://doi.org/10.1073/pnas.94.22.12106
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: X, Yw, Yr
## Flows (time unit: weeks):
##   (outside) -> X   : lambda_in
##   X -> (outside)   : d_x * X
##   X -> Yw          : b_tr * Yw * X
##   X -> Yr          : b_tr * Yr * X
##   Yw -> X          : (r_w + f_tr * h_max * (1 - s_res)) * Yw
##   Yw -> Yr         : f_tr * h_max * s_res * Yw
##   Yr -> X          : r_r * Yr
##   Yw -> (outside)  : c_inf * Yw
##   Yr -> (outside)  : c_inf * Yr

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 100, by = 0.1)  # weeks

# ---- Parameters ------------------------------------------------------
parameters <- c(
  lambda_in = 100,  # entry of uninfected hosts (hosts/time)
  d_x = 1,          # death/removal of uninfected (1/time)
  c_inf = 1.5,      # death/removal of infected hosts (1/time)
  b_tr = 0.05,      # mass-action transmission [illustrative]
  h_max = 1,        # cure rate under treatment (1/time)
  f_tr = 0.5,       # fraction of infected treated [illustrative]
  s_res = 0.001,    # fraction of treated wild-type infections that acquire resistance
  r_w = 0,          # spontaneous recovery, wild-type (1/time)
  r_r = 0.1         # spontaneous recovery, resistant (cost of resistance Δr = r_r − r_w) (1/time)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  X = 99,  # uninfected hosts
  Yw = 1,  # infected with wild-type (sensitive)
  Yr = 0   # infected with resistant strain
)

# ---- Model equations -------------------------------------------------
amr_bonhoeffer_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- X + Yw + Yr

    # Derived quantities
    res_frac <- Yr / (Yw + Yr)  # proportion of infections that are resistant (ρ)

    # Rates of change
    dX <- lambda_in - d_x * X - b_tr * Yw * X - b_tr * Yr * X + (r_w + f_tr * h_max * (1 - s_res)) * Yw + r_r * Yr
    dYw <- b_tr * Yw * X - (r_w + f_tr * h_max * (1 - s_res)) * Yw - f_tr * h_max * s_res * Yw - c_inf * Yw
    dYr <- b_tr * Yr * X + f_tr * h_max * s_res * Yw - r_r * Yr - c_inf * Yr

    list(c(dX, dYw, dYr))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = amr_bonhoeffer_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(X = "#333333", Yw = "#333333", Yr = "#333333")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (weeks)", ylab = "Number of individuals",
        main = "sensitive and resistant strains with treatment and de novo resistance (single drug)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
