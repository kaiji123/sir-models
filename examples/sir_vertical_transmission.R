## SIR with vertical (mother-to-child) transmission
## hepatitis B, rubella (congenital), HTLV, Chagas.
## Vertical transmission adds q_v·b_birth to the infected's "reproductive output", so R0 = (beta + q_vb_birth)/(gamma + mu).
## Parameter values: 0 of 5 checked against the source; the rest are recalled or illustrative.
## Reference: Busenberg S. & Cooke K. (1993) Vertically Transmitted Diseases: Models and Dynamics, Biomathematics 23, Springer
## Source: https://doi.org/10.1007/978-3-642-75301-5
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, R
## Flows (time unit: years):
##   (outside) -> S  : b_birth * (S + R + (1 - q_v) * I)
##   (outside) -> I  : b_birth * q_v * I
##   S -> I          : lambda * S
##   I -> R          : gamma * I
##   S -> (outside)  : mu * S
##   I -> (outside)  : mu * I
##   R -> (outside)  : mu * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 0.2)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  b_birth = 1 / 20,  # per-capita birth rate (1/yr) [illustrative, adjusted for numerical stability]
  mu = 1 / 20,       # natural death (1/yr) [illustrative, adjusted for numerical stability]
  q_v = 0.3,         # probability a child of an infected mother is born infected [illustrative]
  beta = 1.5,        # horizontal transmission (1/yr) [illustrative, adjusted for numerical stability]
  gamma = 0.5        # recovery (1/yr) [illustrative, adjusted for numerical stability]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9900,  # susceptible
  I = 100,   # infectious
  R = 0      # recovered/immune
)

# ---- Model equations -------------------------------------------------
sir_vertical_transmission_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + R

    # Derived quantities
    lambda <- beta * I / N

    # Rates of change
    dS <- b_birth * (S + R + (1 - q_v) * I) - lambda * S - mu * S
    dI <- b_birth * q_v * I + lambda * S - gamma * I - mu * I
    dR <- gamma * I - mu * R

    list(c(dS, dI, dR))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sir_vertical_transmission_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "SIR with vertical (mother-to-child) transmission")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
