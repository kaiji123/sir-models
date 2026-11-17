## Ebola SEIR with intervention-time switch (Congo 1995 / Uganda 2000)
## Ebola virus disease, human.
## The single structural feature that distinguishes this model is the step-then-exponential-decay transmission rate, implemented with ifelse on t.
## Parameter values: 0 of 6 checked against the source; the rest are recalled or illustrative.
## Reference: Chowell G, Hengartner NW, Castillo-Chavez C, Fenimore PW, Hyman JM (2004). The basic reproductive number of Ebola and the effects of public health measures: the cases of Congo and Uganda. J Theor Biol 229(1):119–126
## Source: https://doi.org/10.1016/j.jtbi.2004.03.006
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, I, R, D
## Flows (time unit: days):
##   S -> E  : lambda * S
##   E -> I  : k * E
##   I -> D  : gamma * cfr * I
##   I -> R  : gamma * (1 - cfr) * I

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 250, by = 0.5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta0 = 0.33,    # transmission rate before interventions (Congo 1995; 0.30 commonly quoted for Uganda 2000) [unverified] (per day)
  q = 0.98,        # exponential decay rate of transmission after interventions begin [unverified] (per day)
  tau = 130,       # day interventions begin (relative to epidemic start used in the fit) [unverified] (days)
  k = 0.1887,      # rate of leaving latency (mean incubation 5.3 d) [unverified] (per day)
  gamma = 0.1783,  # rate of leaving infectiousness (mean infectious period 5.61 d) [unverified] (per day)
  cfr = 0.67       # proportion of infectious cases that die (Congo 1995; ≈0.53 Uganda 2000) [unverified] (dimensionless)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 5364500,  # susceptible
  E = 0,        # infected, not yet infectious (latent)
  I = 1,        # infectious (symptomatic, in community)
  R = 0,        # recovered, immune
  D = 0         # died of disease (absorbing; cumulative deaths)
)

# ---- Model equations -------------------------------------------------
ebola_chowell_2004_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N_total <- S + E + I + R

    # Derived quantities
    N <- S + E + I + R                                             # living population
    beta_t <- ifelse(t < tau, beta0, beta0 * exp(-q * (t - tau)))  # time-varying transmission rate
    lambda <- beta_t * I / N                                       # force of infection
    R0 <- beta0 / gamma                                            # basic reproduction number (authors report ≈1.83 Congo 1995, ≈1.34 Uganda 2000) [unverified]
    CC_rate <- (k) * E                                             # rate of new cumulative cases (bookkeeping for fitting to cumulative case curves) (integrate with cumsum for the cumulative count)

    # Rates of change
    dS <- -lambda * S
    dE <- lambda * S - k * E
    dI <- k * E - gamma * cfr * I - gamma * (1 - cfr) * I
    dR <- gamma * (1 - cfr) * I
    dD <- gamma * cfr * I

    list(c(dS, dE, dI, dR, dD), c(CC_rate = CC_rate))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = ebola_chowell_2004_model, parms = parameters, hmax = 0.5)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", I = "#B71C1C", R = "#4E8A3E", D = "#111111")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Ebola SEIR with intervention-time switch (Congo 1995 / Uganda 2000)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("CC_rate")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (days)", ylab = "Value",
        main = "Ebola SEIR with intervention-time switch (Congo 1995 / Uganda 2000): outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
