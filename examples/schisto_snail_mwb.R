## Schistosomiasis: human mean worm burden with SEI snails (Macdonald / Anderson–May hybrid)
## Schistosomiasis (S. haematobium / S. mansoni).
## Parameter values: 0 of 12 checked against the source; the rest are recalled or illustrative.
## Reference: Structure follows the Macdonald (1965) / Anderson & May (1985, 1991) mean-worm-burden + snail-infection framework, as implemented in, for example, Sokolow SH et al. PNAS 2015;112(31):9650–9655
## Source: https://doi.org/10.1073/pnas.1502651112
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Wm, Sn, En, In
## Flows (time unit: days):
##   (outside) -> Sn  : r_s * (1 - Nsn / K_s) * (Sn + En)
##   Sn -> En         : beta_sn * miracidia * Sn
##   En -> In         : sig_n * En
##   Sn -> (outside)  : mu_n * Sn
##   En -> (outside)  : mu_n * En
##   In -> (outside)  : mu_i * In
##   (outside) -> Wm  : lam_h * In
##   Wm -> (outside)  : (mu_wm + mu_h + tau_t) * Wm

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 30, by = 0.05)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  H = 1000,                 # human population, held constant [I] (persons)
  r_s = 0.1,                # snail intrinsic birth rate [M/I] (day^-1)
  K_s = 10000,              # snail carrying capacity [I] (snails)
  mu_n = 1 / 60,            # snail death rate [M] (day^-1)
  sig_n = 1 / 40,           # snail prepatent period (about 4–6 weeks) [M] (day^-1)
  mu_i = 1 / 10,            # infected snail death rate [M] (day^-1)
  beta_sn = 7.5e-7,         # miracidial infection of snails. [I] chosen to give about 5% snail prevalence at Wm = 50. (day^-1 per mated pair)
  lam_h = 1e-4,             # cercarial establishment [I] (worms day^-1 per infected snail per person)
  mu_wm = 1 / (3.3 * 365),  # adult worm death (about 3–4 year life span) [M] (day^-1)
  mu_h = 1 / (60 * 365),    # human death rate [I] (day^-1)
  k_agg = 0.2,              # negative binomial aggregation [M/I]
  tau_t = 0                 # continuous praziquantel treatment (MDA) rate. Pulsed MDA alternative: ifelse((t - 365floor(t/365)) < 7, -log(1 - 0.750.86)/7, 0) [I; 86% worm kill is an M value
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Wm = 50,    # mean adult worm burden per human (worms/person)
  Sn = 7000,  # susceptible snails
  En = 300,   # exposed (prepatent) snails
  In = 400    # infected (cercariae-shedding) snails
)

# ---- Model equations -------------------------------------------------
schisto_snail_mwb_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Derived quantities
    phi_m <- Wm / (Wm + 2 * k_agg + 1)  # crude mating-probability approximation. The exact Macdonald mating function is an integral; this is [I].
    Nsn <- Sn + En + In
    miracidia <- 0.5 * Wm * H * phi_m   # mated pairs, a proxy for egg output

    # Rates of change
    dWm <- lam_h * In - (mu_wm + mu_h + tau_t) * Wm
    dSn <- r_s * (1 - Nsn / K_s) * (Sn + En) - beta_sn * miracidia * Sn - mu_n * Sn
    dEn <- beta_sn * miracidia * Sn - sig_n * En - mu_n * En
    dIn <- sig_n * En - mu_i * In

    list(c(dWm, dSn, dEn, dIn))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = schisto_snail_mwb_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Wm = "#A5D6A7", Sn = "#333333", En = "#E8C468", In = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Schistosomiasis: human mean worm burden with SEI snails (Macdonald / Anderson–May hybrid)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
