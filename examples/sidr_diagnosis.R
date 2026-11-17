## acute infection with diagnosis, isolation-at-home and treatment (no carriers)
## generic acute bacterial or viral infection with case finding (e.g. diphtheria, scabies, early COVID testing).
## R0 is about 1.53 with these values, against 4 without diagnosis.
## Parameter values: 0 of 5 checked against the source; the rest are recalled or illustrative.
## Reference: Illustrative structure in the spirit of the treatment and control chapters of Vynnycky E. & White R.G. (2010) An Introduction to Infectious Disease Modelling, Oxford University Press, ISBN 978-0-19-856576-5, and Martcheva M. (2015) An Introduction to…
## Source: https://doi.org/10.1007/978-1-4899-7612-3
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, Iu, Id, R
## Flows (time unit: days):
##   S -> Iu   : lambda * S
##   Iu -> Id  : delta_d * Iu
##   Iu -> R   : gamma_u * Iu
##   Id -> R   : gamma_d * Id

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 0.2)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 0.4,        # transmission (1/day) [illustrative]
  delta_d = 0.2,     # diagnosis rate (mean 5 days to diagnosis) (1/day) [illustrative]
  kappa_d = 0.2,     # relative infectiousness after diagnosis (isolation plus treatment) [illustrative]
  gamma_u = 1 / 10,  # recovery, undiagnosed (1/day) [illustrative]
  gamma_d = 1 / 3    # recovery, diagnosed and treated (1/day) [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9990,  # susceptible
  Iu = 10,   # infectious, undiagnosed
  Id = 0,    # diagnosed and treated, reduced contacts
  R = 0      # recovered
)

# ---- Model equations -------------------------------------------------
sidr_diagnosis_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + Iu + Id + R

    # Derived quantities
    lambda <- beta * (Iu + kappa_d * Id) / N
    R0 <- beta / (gamma_u + delta_d) + beta * delta_d * kappa_d / ((gamma_u + delta_d) * gamma_d)

    # Rates of change
    dS <- -lambda * S
    dIu <- lambda * S - delta_d * Iu - gamma_u * Iu
    dId <- delta_d * Iu - gamma_d * Id
    dR <- gamma_u * Iu + gamma_d * Id

    list(c(dS, dIu, dId, dR))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sidr_diagnosis_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", Iu = "#B71C1C", Id = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "acute infection with diagnosis, isolation-at-home and treatment (no carriers)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
