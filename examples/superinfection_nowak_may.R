## two-strain superinfection (virulence evolution)
## evolution of virulence; generic parasite.
## Without superinfection (s_sup = 0) the strain with the highest R0 = b·(Λ/μ)/(μ+v) excludes the other.
## Parameter values: 0 of 7 checked against the source; the rest are recalled or illustrative.
## Reference: Nowak M.A. & May R.M. (1994) "Superinfection and the evolution of parasite virulence", Proc. R. Soc. Lond. B 255:81–89
## Source: https://doi.org/10.1098/rspb.1994.0012
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: X, Y1, Y2
## Flows (time unit: years):
##   (outside) -> X   : Lambda
##   X -> Y1          : b1 * Y1 * X
##   X -> Y2          : b2 * Y2 * X
##   Y2 -> Y1         : s_sup * b1 * Y1 * Y2
##   X -> (outside)   : mu * X
##   Y1 -> (outside)  : (mu + v1) * Y1
##   Y2 -> (outside)  : (mu + v2) * Y2

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 0.2)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  Lambda = 100,  # host recruitment (hosts/yr) [illustrative]
  mu = 0.1,      # natural host death (1/yr) [illustrative]
  b1 = 0.005,    # transmission, strain 1 [illustrative]
  b2 = 0.004,    # transmission, strain 2 [illustrative]
  v1 = 0.5,      # virulence (extra death), strain 1 (1/yr) [illustrative]
  v2 = 0.2,      # virulence, strain 2 (1/yr) [illustrative]
  s_sup = 1.0    # superinfection coefficient (strain 1 superinfecting Y2 hosts) [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  X = 990,  # uninfected hosts
  Y1 = 5,   # infected with the more virulent strain 1
  Y2 = 5    # infected with the less virulent strain 2
)

# ---- Model equations -------------------------------------------------
superinfection_nowak_may_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- X + Y1 + Y2

    # Rates of change
    dX <- Lambda - b1 * Y1 * X - b2 * Y2 * X - mu * X
    dY1 <- b1 * Y1 * X + s_sup * b1 * Y1 * Y2 - (mu + v1) * Y1
    dY2 <- b2 * Y2 * X - s_sup * b1 * Y1 * Y2 - (mu + v2) * Y2

    list(c(dX, dY1, dY2))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = superinfection_nowak_may_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(X = "#333333", Y1 = "#333333", Y2 = "#333333")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "two-strain superinfection (virulence evolution)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
