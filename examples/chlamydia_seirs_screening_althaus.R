## Chlamydia SEIRS, two sexes, symptomatic/asymptomatic, with screening
## Chlamydia trachomatis.
## Parameter values: 1 of 12 checked against the source; the rest are recalled or illustrative.
## Reference: Althaus CL, Heijne JCM, Roellin A, Low N. Transmission dynamics of Chlamydia trachomatis affect the impact of screening programmes. Epidemics 2010;2(3):123–131
## Source: https://doi.org/10.1016/j.epidem.2010.04.002
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sf, Ef, Iaf, Isf, Rf, Sm, Em, Iam, Ism, Rm
## Flows (time unit: days):
##   Sf -> Ef   : foi_f * Sf
##   Ef -> Isf  : p_sym_f * sig_inc * Ef
##   Ef -> Iaf  : (1 - p_sym_f) * sig_inc * Ef
##   Iaf -> Rf  : 1 / d_af * Iaf
##   Iaf -> Sf  : scr_f * Iaf
##   Isf -> Sf  : 1 / d_s * Isf
##   Rf -> Sf   : w_imm * Rf
##   Sm -> Em   : foi_m * Sm
##   Em -> Ism  : p_sym_m * sig_inc * Em
##   Em -> Iam  : (1 - p_sym_m) * sig_inc * Em
##   Iam -> Rm  : 1 / d_am * Iam
##   Iam -> Sm  : scr_m * Iam
##   Ism -> Sm  : 1 / d_s * Ism
##   Rm -> Sm   : w_imm * Rm

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 3650, by = 5)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  cp = 1 / 365,       # partner change rate (1 per year, equal by sex to keep partnerships balanced) (partners/day) [illustrative]
  beta_mf = 0.3,      # male-to-female transmission (per partnership) [illustrative]
  beta_fm = 0.3,      # female-to-male transmission (per partnership) [illustrative]
  sig_inc = 1 / 14,   # exit from latent period (/day) [illustrative]
  p_sym_f = 0.1,      # symptomatic fraction, women [illustrative]
  p_sym_m = 0.2,      # symptomatic fraction, men [illustrative]
  d_af = 433,         # asymptomatic duration, women (days) [verified]
  d_am = 433,         # asymptomatic duration, men (days) [illustrative]
  d_s = 30,           # time to treatment when symptomatic (days) [illustrative]
  w_imm = 1 / 365,    # waning of immunity (R back to S) (/day) [illustrative]
  scr_f = 0.2 / 365,  # screening plus treatment rate, women (/day) [illustrative]
  scr_m = 0.05 / 365  # screening plus treatment rate, men (/day) [illustrative]
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sf = 0.95,   # susceptible women
  Ef = 0,      # exposed (latent) women
  Iaf = 0.05,  # asymptomatic infected women
  Isf = 0,     # symptomatic infected women
  Rf = 0,      # immune women
  Sm = 0.95,   # susceptible men
  Em = 0,      # exposed men
  Iam = 0.05,  # asymptomatic infected men
  Ism = 0,     # symptomatic infected men
  Rm = 0       # immune men
)

# ---- Model equations -------------------------------------------------
chlamydia_seirs_screening_althaus_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sf + Ef + Iaf + Isf + Rf + Sm + Em + Iam + Ism + Rm

    # Derived quantities
    Nf <- Sf + Ef + Iaf + Isf + Rf            # women total
    Nm <- Sm + Em + Iam + Ism + Rm            # men total
    foi_f <- cp * beta_mf * (Iam + Ism) / Nm  # force of infection on women
    foi_m <- cp * beta_fm * (Iaf + Isf) / Nf  # force of infection on men

    # Rates of change
    dSf <- -foi_f * Sf + scr_f * Iaf + 1 / d_s * Isf + w_imm * Rf
    dEf <- foi_f * Sf - p_sym_f * sig_inc * Ef - (1 - p_sym_f) * sig_inc * Ef
    dIaf <- (1 - p_sym_f) * sig_inc * Ef - 1 / d_af * Iaf - scr_f * Iaf
    dIsf <- p_sym_f * sig_inc * Ef - 1 / d_s * Isf
    dRf <- 1 / d_af * Iaf - w_imm * Rf
    dSm <- -foi_m * Sm + scr_m * Iam + 1 / d_s * Ism + w_imm * Rm
    dEm <- foi_m * Sm - p_sym_m * sig_inc * Em - (1 - p_sym_m) * sig_inc * Em
    dIam <- (1 - p_sym_m) * sig_inc * Em - 1 / d_am * Iam - scr_m * Iam
    dIsm <- p_sym_m * sig_inc * Em - 1 / d_s * Ism
    dRm <- 1 / d_am * Iam - w_imm * Rm

    list(c(dSf, dEf, dIaf, dIsf, dRf, dSm, dEm, dIam, dIsm, dRm))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = chlamydia_seirs_screening_althaus_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sf = "#333333", Ef = "#E8C468", Iaf = "#B71C1C", Isf = "#B71C1C", Rf = "#4E8A3E", Sm = "#333333", Em = "#E8C468", Iam = "#B71C1C", Ism = "#B71C1C", Rm = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Chlamydia SEIRS, two sexes, symptomatic/asymptomatic, with screening")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
