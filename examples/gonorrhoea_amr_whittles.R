## Gonorrhoea in MSM with cefixime-sensitive and cefixime-resistant strains
## Gonorrhoea, antimicrobial resistance.
## Parameter values: 8 of 9 checked against the source; the rest are recalled or illustrative.
## Reference: Whittles LK, White PJ, Didelot X. Estimating the fitness cost and benefit of cefixime resistance in Neisseria gonorrhoeae to inform prescription policy: a modelling study. PLoS Med 2017;14(10):e1002416
## Source: https://doi.org/10.1371/journal.pmed.1002416
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, Us, Es, As, Tsc, Tso, Ur, Er, Ar, Trc, Tro
## Flows (time unit: years):
##   S -> Us    : foi_s * S
##   S -> Ur    : foi_r * S
##   Us -> Es   : psi * sig_inc * Us
##   Us -> As   : (1 - psi) * sig_inc * Us
##   Ur -> Er   : psi * sig_inc * Ur
##   Ur -> Ar   : (1 - psi) * sig_inc * Ur
##   As -> S    : nu * As
##   Ar -> S    : alpha_r * nu * Ar
##   Es -> Tsc  : pi_cef * mu_t * Es
##   Es -> Tso  : (1 - pi_cef) * mu_t * Es
##   Er -> Trc  : pi_cef * mu_t * Er
##   Er -> Tro  : (1 - pi_cef) * mu_t * Er
##   Tsc -> S   : rho * Tsc
##   Tso -> S   : rho * Tso
##   Tro -> S   : rho * Tro
##   Trc -> S   : (1 - phi_f) * rho * Trc
##   Trc -> Ar  : phi_f * rho * Trc

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 10, by = 0.01)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  theta = 5.2,    # transmission rate (/yr) [verified]
  psi = 0.63,     # proportion of infections becoming symptomatic [verified]
  sig_inc = 77,   # rate of leaving incubation (about 4.7 days) (/yr) [verified]
  nu = 1.8,       # natural recovery from asymptomatic infection, sensitive strain (/yr) [verified]
  alpha_r = 1.8,  # relative natural recovery rate of the resistant strain (fitness cost) [verified]
  mu_t = 136,     # treatment-seeking rate when symptomatic (/yr) [verified]
  rho = 53,       # recovery rate after treatment (/yr) [verified]
  phi_f = 0.83,   # cefixime treatment failure for resistant infection [verified]
  pi_cef = 0.3    # proportion treated with cefixime [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  S = 600000,  # susceptible MSM
  Us = 0,      # incubating, sensitive strain
  Es = 0,      # symptomatic, sensitive
  As = 618,    # asymptomatic, sensitive
  Tsc = 0,     # treated with cefixime, sensitive
  Tso = 0,     # treated with another antibiotic, sensitive
  Ur = 0,      # incubating, resistant strain
  Er = 0,      # symptomatic, resistant
  Ar = 49,     # asymptomatic, resistant
  Trc = 0,     # treated with cefixime, resistant
  Tro = 0      # treated with another antibiotic, resistant
)

# ---- Model equations -------------------------------------------------
gonorrhoea_amr_whittles_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + Us + Es + As + Tsc + Tso + Ur + Er + Ar + Trc + Tro

    # Derived quantities
    Inf_s <- Us + Es + As + Tsc + Tso  # infected with sensitive strain
    Inf_r <- Ur + Er + Ar + Trc + Tro  # infected with resistant strain
    foi_s <- theta * Inf_s / N         # force of infection, sensitive
    foi_r <- theta * Inf_r / N         # force of infection, resistant

    # Rates of change
    dS <- -foi_s * S - foi_r * S + nu * As + alpha_r * nu * Ar + rho * Tsc + rho * Tso + rho * Tro + (1 - phi_f) * rho * Trc
    dUs <- foi_s * S - psi * sig_inc * Us - (1 - psi) * sig_inc * Us
    dEs <- psi * sig_inc * Us - pi_cef * mu_t * Es - (1 - pi_cef) * mu_t * Es
    dAs <- (1 - psi) * sig_inc * Us - nu * As
    dTsc <- pi_cef * mu_t * Es - rho * Tsc
    dTso <- (1 - pi_cef) * mu_t * Es - rho * Tso
    dUr <- foi_r * S - psi * sig_inc * Ur - (1 - psi) * sig_inc * Ur
    dEr <- psi * sig_inc * Ur - pi_cef * mu_t * Er - (1 - pi_cef) * mu_t * Er
    dAr <- (1 - psi) * sig_inc * Ur - alpha_r * nu * Ar + phi_f * rho * Trc
    dTrc <- pi_cef * mu_t * Er - (1 - phi_f) * rho * Trc - phi_f * rho * Trc
    dTro <- (1 - pi_cef) * mu_t * Er - rho * Tro

    list(c(dS, dUs, dEs, dAs, dTsc, dTso, dUr, dEr, dAr, dTrc, dTro))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = gonorrhoea_amr_whittles_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", Us = "#333333", Es = "#E8C468", As = "#E57373", Tsc = "#333333", Tso = "#333333", Ur = "#333333", Er = "#E8C468", Ar = "#E57373", Trc = "#333333", Tro = "#333333")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Gonorrhoea in MSM with cefixime-sensitive and cefixime-resistant strains")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
