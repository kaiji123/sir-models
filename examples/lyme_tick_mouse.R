## Lyme borreliosis: stage-structured Ixodes ticks with a mouse reservoir and deer as tick-only hosts [ILLUSTRATIVE-COMPOSITE]
## Lyme disease (Borrelia burgdorferi s.l.), temperate woodland in NE USA or Europe.
## Parameter values: 0 of 15 checked against the source; the rest are recalled or illustrative.
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: L, Ns, Ni, As, Ai, Ms, Mi
## Flows (time unit: days):
##   (outside) -> L   : eps_egg * a_A * Deer * (As + Ai)
##   L -> Ns          : a_L * s_L * (Hm - p_ML * Mi) * L
##   L -> Ni          : a_L * s_L * p_ML * Mi * L
##   L -> (outside)   : (a_L * Hm * (1 - s_L) + mu_L + mu_L2 * L) * L
##   Ns -> As         : a_N * s_N * (Hm - p_ML * Mi) * Ns
##   Ns -> Ai         : a_N * s_N * p_ML * Mi * Ns
##   Ns -> (outside)  : (a_N * Hm * (1 - s_N) + mu_N) * Ns
##   Ni -> Ai         : a_N * s_N * Hm * Ni
##   Ni -> (outside)  : (a_N * Hm * (1 - s_N) + mu_N) * Ni
##   As -> (outside)  : (a_A * Deer + mu_A) * As
##   Ai -> (outside)  : (a_A * Deer + mu_A) * Ai
##   (outside) -> Ms  : mu_m * K_m
##   Ms -> Mi         : a_N * p_NM * Ni * Ms
##   Ms -> (outside)  : mu_m * Ms
##   Mi -> (outside)  : mu_m * Mi

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 7300, by = 10)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  a_L = 1e-4,      # larval host-finding rate per mouse
  a_N = 1e-4,      # nymphal host-finding rate per mouse
  a_A = 0.005,     # adult host-finding rate per deer
  Deer = 5,        # deer density, constant; not competent for Borrelia (deer)
  s_L = 0.3,       # survival from fed larva to nymph
  s_N = 0.3,       # survival from fed nymph to adult
  eps_egg = 1000,  # effective larvae per fed adult (larvae/fed adult)
  mu_L = 0.01,     # larval death (1/day)
  mu_L2 = 4.8e-5,  # density-dependent larval death
  mu_N = 0.005,    # nymphal death (1/day)
  mu_A = 0.005,    # adult death (1/day)
  p_NM = 0.8,      # infected nymph → mouse transmission probability
  p_ML = 0.5,      # infected mouse → feeding larva transmission probability
  mu_m = 0.005,    # mouse death (1/day)
  K_m = 100        # mouse population; recruitment = mu_m·K_m (mice)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  L = 10000,  # questing larvae (born uninfected; no transovarial transmission)
  Ns = 1000,  # uninfected questing nymphs
  Ni = 0,     # infected questing nymphs
  As = 100,   # uninfected questing adults
  Ai = 0,     # infected adults
  Ms = 95,    # susceptible mice (reservoir)
  Mi = 5      # infected mice (lifelong)
)

# ---- Model equations -------------------------------------------------
lyme_tick_mouse_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Ms + Mi

    # Derived quantities
    Hm <- Ms + Mi          # mice
    Nticks <- L + Ns + Ni + As + Ai
    prev_mouse <- Mi / Hm  # mouse infection prevalence
    nip <- Ni / (Ns + Ni)  # nymphal infection prevalence (the human-risk index)

    # Rates of change
    dL <- eps_egg * a_A * Deer * (As + Ai) - a_L * s_L * (Hm - p_ML * Mi) * L - a_L * s_L * p_ML * Mi * L - (a_L * Hm * (1 - s_L) + mu_L + mu_L2 * L) * L
    dNs <- a_L * s_L * (Hm - p_ML * Mi) * L - a_N * s_N * (Hm - p_ML * Mi) * Ns - a_N * s_N * p_ML * Mi * Ns - (a_N * Hm * (1 - s_N) + mu_N) * Ns
    dNi <- a_L * s_L * p_ML * Mi * L - a_N * s_N * Hm * Ni - (a_N * Hm * (1 - s_N) + mu_N) * Ni
    dAs <- a_N * s_N * (Hm - p_ML * Mi) * Ns - (a_A * Deer + mu_A) * As
    dAi <- a_N * s_N * p_ML * Mi * Ns + a_N * s_N * Hm * Ni - (a_A * Deer + mu_A) * Ai
    dMs <- mu_m * K_m - a_N * p_NM * Ni * Ms - mu_m * Ms
    dMi <- a_N * p_NM * Ni * Ms - mu_m * Mi

    list(c(dL, dNs, dNi, dAs, dAi, dMs, dMi))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = lyme_tick_mouse_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(L = "#E8C468", Ns = "#333333", Ni = "#333333", As = "#E57373", Ai = "#E57373", Ms = "#F3C6D3", Mi = "#F3C6D3")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Lyme borreliosis: stage-structured Ixodes ticks with a mouse reservoir and deer as tick-only hosts [ILLUSTRATIVE-COMPOSITE]")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
