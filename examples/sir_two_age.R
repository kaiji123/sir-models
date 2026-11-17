## SIR with two age groups
## Children and adults with a who-acquires-infection-from-whom (WAIFW) matrix and ageing.
## Reference: Anderson & May (1991), Infectious Diseases of Humans, ch. 9; Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals, §3.1
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sc, Ic, Rc, Sa, Ia, Ra
## Flows (time unit: years):
##   (outside) -> Sc  : mu * N
##   Sc -> Ic         : lam_c * Sc
##   Ic -> Rc         : gamma * Ic
##   Sa -> Ia         : lam_a * Sa
##   Ia -> Ra         : gamma * Ia
##   Sc -> Sa         : ageing * Sc
##   Ic -> Ia         : ageing * Ic
##   Rc -> Ra         : ageing * Rc
##   Sa -> (outside)  : mu_a * Sa
##   Ia -> (outside)  : mu_a * Ia
##   Ra -> (outside)  : mu_a * Ra

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 50, by = 0.02)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  b_cc = 1000,      # transmission child -> child per year
  b_ca = 300,       # transmission adult -> child per year
  b_ac = 300,       # transmission child -> adult per year
  b_aa = 400,       # transmission adult -> adult per year
  gamma = 365 / 7,  # recovery rate per year
  ageing = 1 / 15,  # rate of moving from child to adult per year
  mu = 1 / 70       # birth rate per year (deaths in adults keep N constant)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sc = 2000,   # Susceptible children
  Ic = 10,     # Infectious children
  Rc = 17990,  # Recovered children
  Sa = 1000,   # Susceptible adults
  Ia = 0,      # Infectious adults
  Ra = 79000   # Recovered adults
)

# ---- Model equations -------------------------------------------------
sir_two_age_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sc + Ic + Rc + Sa + Ia + Ra

    # Derived quantities
    lam_c <- (b_cc * Ic + b_ca * Ia) / N  # force of infection on children
    lam_a <- (b_ac * Ic + b_aa * Ia) / N  # force of infection on adults
    mu_a <- mu * N / (Sa + Ia + Ra)       # adult death rate balancing births

    # Rates of change
    dSc <- mu * N - lam_c * Sc - ageing * Sc
    dIc <- lam_c * Sc - gamma * Ic - ageing * Ic
    dRc <- gamma * Ic - ageing * Rc
    dSa <- -lam_a * Sa + ageing * Sc - mu_a * Sa
    dIa <- lam_a * Sa - gamma * Ia + ageing * Ic - mu_a * Ia
    dRa <- gamma * Ia + ageing * Rc - mu_a * Ra

    list(c(dSc, dIc, dRc, dSa, dIa, dRa))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = sir_two_age_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sc = "#333333", Ic = "#B71C1C", Rc = "#4E8A3E", Sa = "#333333", Ia = "#B71C1C", Ra = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "SIR with two age groups")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
