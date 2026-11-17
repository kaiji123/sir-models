## Codeço endemic/epidemic cholera with aquatic reservoir and dose–response
## Cholera (Vibrio cholerae).
## Parameter values: 6 of 6 checked against the source; the rest are recalled or illustrative.
## Reference: Codeço CT. Endemic and epidemic dynamics of cholera: the role of the aquatic reservoir. BMC Infect Dis 2001;1:1
## Source: https://doi.org/10.1186/1471-2334-1-1
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, I, R, B
## Flows (time unit: days):
##   (outside) -> S  : mu_h * N
##   S -> (outside)  : mu_h * S
##   S -> I          : a * B / (K + B) * S
##   I -> R          : r * I
##   I -> (outside)  : mu_h * I
##   R -> (outside)  : mu_h * R
##   (outside) -> B  : e * I
##   B -> (outside)  : (-nb_mb) * B

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 730, by = 1)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  mu_h = 0.0001,  # human birth = death rate (the paper's n) (day^-1) [verified]
  a = 1,          # rate of exposure to contaminated water (day^-1) [verified]
  K = 1e6,        # concentration giving a 50% chance of infection (cells/ml) [verified]
  r = 0.2,        # recovery rate (5-day infection) (day^-1) [verified]
  nb_mb = -0.33,  # net growth rate of vibrios, nb − mb. It is negative, so the population cannot persist without human input (day^-1) [verified]
  e = 10          # contribution of each infected person to the aquatic population (cells ml^-1 day^-1 person^-1) [verified]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 9999,  # susceptible
  I = 1,     # infected
  R = 0,     # recovered/removed. Added for bookkeeping: the original has only S and I with constant total H, and recovered people are
  B = 0      # toxigenic V. cholerae concentration in water (cells/ml)
)

# ---- Model equations -------------------------------------------------
cholera_codeco2001_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + I + R

    # Derived quantities
    lam_B <- B / (K + B)                  # probability of infection on exposure (saturating dose–response)
    S_crit <- K * r * (-nb_mb) / (e * a)  # critical susceptible number for invasion: 6600 for Community 2 and 13,200 for Community 1. These match the paper's Table 2 S_C values, which

    # Rates of change
    dS <- mu_h * N - mu_h * S - a * B / (K + B) * S
    dI <- a * B / (K + B) * S - r * I - mu_h * I
    dR <- r * I - mu_h * R
    dB <- e * I - (-nb_mb) * B

    list(c(dS, dI, dR, dB))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = cholera_codeco2001_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", I = "#B71C1C", R = "#4E8A3E", B = "#333333")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Codeço endemic/epidemic cholera with aquatic reservoir and dose–response")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
