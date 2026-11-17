## Norovirus SEIAR with asymptomatic reinfection and waning disease immunity
## Norovirus gastroenteritis.
## Parameter values: 0 of 8 checked against the source; the rest are recalled or illustrative.
## Reference: Simmons K, Gambhir M, Leon J, Lopman B. Duration of immunity to norovirus gastroenteritis. Emerg Infect Dis 2013;19(8):1260–1267
## Source: https://doi.org/10.3201/eid1908.130472
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, E, I, A, R
## Flows (time unit: days):
##   (outside) -> S  : mu_h * N
##   S -> E          : lam * S
##   E -> I          : 1 / inc_d * E
##   I -> A          : 1 / sym_d * I
##   A -> R          : 1 / shed_d * A
##   R -> A          : lam * R
##   R -> S          : 1 / theta_d * R

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 30, by = 0.05)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  inc_d = 1,              # incubation, E → I [V, "about 1 day"] (day)
  sym_d = 2,              # symptomatic duration, I → A [V, "1 day, half the symptomatic phase"] (day)
  shed_d = 10,            # asymptomatic shedding, A → R [V, Model B] (day)
  eps_ea = 0.05,          # relative infectiousness of E and A: 0.05 in Model B, 0.25 in Model C, 0 in Model A [V]
  theta_d = 5.1 * 365,    # duration of disease immunity, R → S. Fitted values: 5.1 years (Model A/B), 8.7 (C), 4.1 (D, 20% genetically resistant) [V] (day)
  R0_nv = 1.7,            # R0. Reported values: 1.64–1.88 for most models, 7.16 for C [V]
  mu_h = 1 / (70 * 365),  # birth = death [I] (day^-1)
  N0 = 100000             # total initial population [illustrative default]
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  S <- 0.3 * N0  # susceptible to infection and disease
  E <- 0         # exposed (presymptomatic)
  I <- 100       # symptomatic
  A <- 0         # asymptomatic infected (post-symptomatic shedding, or reinfection of immune people)
  R <- 0.7 * N0  # immune to disease but not to infection
  c(S = S, E = E, I = I, A = A, R = R)
})

# ---- Model equations -------------------------------------------------
norovirus_simmons2013_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + E + I + A + R

    # Derived quantities
    Dinf <- eps_ea * inc_d + sym_d + eps_ea * shed_d  # effective infectious duration (2.55 days)
    b_nv <- R0_nv / Dinf                              # about 0.667/day
    lam <- b_nv * (I + eps_ea * (E + A)) / N

    # Rates of change
    dS <- mu_h * N - lam * S + 1 / theta_d * R
    dE <- lam * S - 1 / inc_d * E
    dI <- 1 / inc_d * E - 1 / sym_d * I
    dA <- 1 / sym_d * I - 1 / shed_d * A + lam * R
    dR <- 1 / shed_d * A - lam * R - 1 / theta_d * R

    list(c(dS, dE, dI, dA, dR))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = norovirus_simmons2013_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", E = "#E8C468", I = "#B71C1C", A = "#E57373", R = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Norovirus SEIAR with asymptomatic reinfection and waning disease immunity")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
