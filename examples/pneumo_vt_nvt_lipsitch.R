## Pneumococcal carriage: vaccine-type (VT) and non-vaccine-type (NVT) serotypes with competition (Lipsitch 1997)
## Streptococcus pneumoniae carriage.
## Parameter values: 0 of 8 checked against the source; the rest are recalled or illustrative.
## Reference: Lipsitch M. "Vaccination against colonizing bacteria with multiple serotypes." PNAS 1997;94:6571–6576. https://doi.org/10.1073/pnas.94.12.6571 (PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC21091)
## Source: https://doi.org/10.1073/pnas.94.12.6571
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: X, Y1, Y2, Y12, Xv, Y1v, Y2v, Y12v
## Flows (time unit: days):
##   (outside) -> X   : (1 - p_vac) * mu * N
##   (outside) -> Xv  : p_vac * mu * N
##   X -> Y1          : foi1 * X
##   X -> Y2          : foi2 * X
##   Y1 -> Y12        : comp * foi2 * Y1
##   Y2 -> Y12        : comp * foi1 * Y2
##   Y1 -> X          : clr * Y1
##   Y2 -> X          : clr * Y2
##   Y12 -> Y2        : clr * Y12
##   Y12 -> Y1        : clr * Y12
##   Xv -> Y1v        : (1 - ve) * foi1 * Xv
##   Xv -> Y2v        : foi2 * Xv
##   Y1v -> Y12v      : comp * foi2 * Y1v
##   Y2v -> Y12v      : comp * (1 - ve) * foi1 * Y2v
##   Y1v -> Xv        : clr * Y1v
##   Y2v -> Xv        : clr * Y2v
##   Y12v -> Y2v      : clr * Y12v
##   Y12v -> Y1v      : clr * Y12v
##   X -> (outside)   : mu * X

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 3650, by = 5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  b_vt = 0.06,         # acquisition rate per carrier contact [illustrative]
  b_nvt = 0.055,       # acquisition rate per carrier contact (1/day) [illustrative]
  clr = 1 / 42,        # clearance of each serotype (carriage of about 6 weeks) (1/day) [illustrative]
  comp = 0.3,          # relative rate at which a carrier acquires the second serotype (1 = no competition) [illustrative]
  ve = 0.9,            # vaccine reduction in VT acquisition [illustrative]
  p_vac = 0.9,         # fraction of births vaccinated [illustrative]
  mu = 1 / (5 * 365),  # turnover of the at-risk population (children under 5) (1/day) [illustrative]
  N0 = 1e5             # value taken from the notes text [recalled]
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  X <- 0.6 * N0    # unvaccinated, uncolonised
  Y1 <- 0.15 * N0  # unvaccinated, carrying VT
  Y2 <- 0.15 * N0  # unvaccinated, carrying NVT
  Y12 <- 0.1 * N0  # unvaccinated, carrying both
  Xv <- 0          # the same four states in vaccinated people
  Y1v <- 0         # the same four states in vaccinated people
  Y2v <- 0         # the same four states in vaccinated people
  Y12v <- 0        # the same four states in vaccinated people
  c(X = X, Y1 = Y1, Y2 = Y2, Y12 = Y12, Xv = Xv, Y1v = Y1v, Y2v = Y2v, Y12v = Y12v)
})

# ---- Model equations -------------------------------------------------
pneumo_vt_nvt_lipsitch_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- X + Y1 + Y2 + Y12 + Xv + Y1v + Y2v + Y12v

    # Derived quantities
    foi1 <- b_vt * (Y1 + Y12 + Y1v + Y12v) / N
    foi2 <- b_nvt * (Y2 + Y12 + Y2v + Y12v) / N

    # Rates of change
    dX <- (1 - p_vac) * mu * N - foi1 * X - foi2 * X + clr * Y1 + clr * Y2 - mu * X
    dY1 <- foi1 * X - comp * foi2 * Y1 - clr * Y1 + clr * Y12
    dY2 <- foi2 * X - comp * foi1 * Y2 - clr * Y2 + clr * Y12
    dY12 <- comp * foi2 * Y1 + comp * foi1 * Y2 - clr * Y12 - clr * Y12
    dXv <- p_vac * mu * N - (1 - ve) * foi1 * Xv - foi2 * Xv + clr * Y1v + clr * Y2v
    dY1v <- (1 - ve) * foi1 * Xv - comp * foi2 * Y1v - clr * Y1v + clr * Y12v
    dY2v <- foi2 * Xv - comp * (1 - ve) * foi1 * Y2v - clr * Y2v + clr * Y12v
    dY12v <- comp * foi2 * Y1v + comp * (1 - ve) * foi1 * Y2v - clr * Y12v - clr * Y12v

    list(c(dX, dY1, dY2, dY12, dXv, dY1v, dY2v, dY12v))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = pneumo_vt_nvt_lipsitch_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(X = "#333333", Y1 = "#333333", Y2 = "#333333", Y12 = "#333333", Xv = "#333333", Y1v = "#333333", Y2v = "#333333", Y12v = "#333333")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Pneumococcal carriage: vaccine-type (VT) and non-vaccine-type (NVT) serotypes with competition (Lipsitch 1997)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
