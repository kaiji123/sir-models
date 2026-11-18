## Chronic wasting disease in deer with environmental prion accumulation
## Chronic wasting disease (prion), mule deer / white-tailed deer.
## Only 3 compartments, which makes this the cleanest demonstration of an environmental-reservoir model: with beta_env = 0 the epidemic is density-dependent and self-limiting, while with slow prion decay the environmental term creates a ratche
## Parameter values: 0 of 10 checked against the source; the rest are recalled or illustrative.
## Reference: Miller MW, Hobbs NT, Tavener SJ (2006). Dynamics of prion disease transmission in mule deer. Ecological Applications 16(6):2208–2214
## Source: https://doi.org/10.1890/1051-0761(2006)016[2208:DOPDTI]2.0.CO;2
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, Ii, Ev
## Flows (time unit: years):
##   (outside) -> S   : births
##   S -> Ii          : lambda * S
##   S -> (outside)   : (mu + harvest) * S
##   Ii -> (outside)  : (mu + alpha + harvest) * Ii
##   (outside) -> Ev  : (shed * Ii + carc * alpha * Ii)
##   Ev -> (outside)  : decay_e * Ev

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 100, by = 0.1)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta_dir = 0.0008,  # direct deer-to-deer transmission coefficient (density-dependent) [illustrative] (per deer per year)
  beta_env = 0.0002,  # environment-to-deer transmission coefficient [illustrative] (per prion unit per year)
  alpha = 0.55,       # CWD-induced mortality rate (mean time from infection to death ≈2 years) [unverified] (per year)
  mu = 0.15,          # natural (non-CWD) mortality rate [illustrative] (per year)
  bdeer = 0.35,       # maximum per-capita birth rate [illustrative] (per year)
  Kd = 1500,          # carrying capacity [illustrative] (deer)
  shed = 1.0,         # shedding from live infected animals (saliva, urine, faeces) [illustrative] (prion units per infected deer per year)
  carc = 10.0,        # prion deposited by a carcass [illustrative] (prion units per CWD death)
  decay_e = 0.05,     # loss of environmental infectivity (prions persist for years; this is the parameter the literature identifies as pivotal) [illustrative] (per year)
  harvest = 0         # hunter harvest rate (control lever; can be applied preferentially to Ii) [illustrative] (per year)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 1000,  # susceptible deer
  Ii = 1,    # infected deer (long subclinical + clinical infection, assumed infectious throughout)
  Ev = 0     # environmental prion load (infectious units; E avoided for clarity)
)

# ---- Model equations -------------------------------------------------
cwd_deer_environmental_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N_total <- S + Ii

    # Derived quantities
    N <- S + Ii                                # deer population
    lambda <- beta_dir * Ii + beta_env * Ev    # force of infection (direct + indirect, additive)
    births <- bdeer * N * pmax(1 - N / Kd, 0)  # logistic recruitment
    prev <- Ii / pmax(N, 1e-9)                 # prevalence

    # Rates of change
    dS <- births - lambda * S - (mu + harvest) * S
    dIi <- lambda * S - (mu + alpha + harvest) * Ii
    dEv <- (shed * Ii + carc * alpha * Ii) - decay_e * Ev

    list(c(dS, dIi, dEv))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = cwd_deer_environmental_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", Ii = "#B71C1C", Ev = "#E8C468")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Chronic wasting disease in deer with environmental prion accumulation")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
