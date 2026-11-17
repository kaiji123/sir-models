## MSIR (maternal immunity)
## Newborns are protected by maternal antibodies (M) before becoming susceptible.
## Reference: Hethcote (2000), SIAM Review 42:599-653
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: M, S, I, R
## Flows (time unit: years):
##   (outside) -> M  : mu * N
##   M -> S          : delta * M
##   S -> I          : lam * S
##   I -> R          : gamma * I
##   M -> (outside)  : mu * M
##   S -> (outside)  : mu * S
##   I -> (outside)  : mu * I
##   R -> (outside)  : mu * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 60, by = 0.02)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta = 700,       # transmission rate per year
  delta = 2,        # rate of losing maternal immunity per year (6 months)
  gamma = 365 / 7,  # recovery rate per year
  mu = 1 / 70       # birth = death rate per year
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  M = 500,   # Maternally immune
  S = 9500,  # Susceptible
  I = 10,    # Infectious
  R = 89990  # Recovered
)

# ---- Model equations -------------------------------------------------
msir_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- M + S + I + R

    # Derived quantities
    lam <- beta * I / N  # force of infection

    # Flows recorded as outputs
    incidence <- lam * S

    # Rates of change
    dM <- mu * N - delta * M - mu * M
    dS <- delta * M - lam * S - mu * S
    dI <- lam * S - gamma * I - mu * I
    dR <- gamma * I - mu * R

    list(c(dM, dS, dI, dR), c(incidence = incidence))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = msir_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(M = "#F3C6D3", S = "#333333", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "MSIR (maternal immunity)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("incidence")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (years)", ylab = "Value",
        main = "MSIR (maternal immunity): outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
