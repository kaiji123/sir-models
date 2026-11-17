## Brucellosis in small ruminants with an environmental reservoir and human cases
## Brucellosis, sheep/goats → humans, with Brucella in the environment.
## The structural signature of this family is the additive direct + indirect force of infection with an explicit environmental bacterial compartment W, and the result that vaccination and disinfection of adults alone cannot eradicate because y
## Parameter values: 0 of 17 checked against the source; the rest are recalled or illustrative.
## Reference: Hou Q, Sun X, Zhang J, Liu Y, Wang Y, Jin Z (2013). Modeling the transmission dynamics of sheep brucellosis in Inner Mongolia Autonomous Region, China. Mathematical Biosciences 242(1):51–58
## Source: https://doi.org/10.1016/j.mbs.2012.11.012
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Ss, Es, Is, Vs, W, Sh, Eh, Iah, Ich
## Flows (time unit: years):
##   (outside) -> Ss  : As
##   Ss -> Es         : lambda_s * Ss
##   Ss -> Vs         : nu_s * Ss
##   Vs -> Ss         : wane_s * Vs
##   Es -> Is         : sigma_s * Es
##   Ss -> (outside)  : mu_s * Ss
##   (outside) -> W   : shed * Is
##   W -> (outside)   : (decay + disinf) * W
##   Sh -> Eh         : lambda_h * Sh
##   Eh -> Iah        : sigma_h * Eh
##   Iah -> Ich       : chron * Iah
##   Iah -> Sh        : gamma_a * Iah
##   Ich -> Sh        : gamma_c * Ich
##   Sh -> (outside)  : mu_h * Sh

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 30, by = 0.05)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  As = 2.5e5,      # sheep recruitment (lambs + purchases) [illustrative] (sheep per year)
  mu_s = 0.4,      # sheep removal rate (slaughter + natural death) [illustrative] (per year)
  beta_ss = 2e-7,  # direct sheep-to-sheep transmission coefficient [illustrative] (per sheep per year)
  beta_sw = 1e-6,  # indirect (environment-to-sheep) transmission coefficient [illustrative] (per unit load per year)
  sigma_s = 2.0,   # rate exposed sheep become infectious (≈6-month latency) [illustrative] (per year)
  nu_s = 0,        # sheep vaccination rate (control lever; S19/Rev-1) [illustrative] (per year)
  wane_s = 0.5,    # waning of sheep vaccine protection [illustrative] (per year)
  shed = 15,       # environmental shedding rate [illustrative] (load units per infectious sheep per year)
  decay = 3.0,     # natural decay of Brucella in the environment [illustrative] (per year)
  disinf = 0,      # additional removal of environmental load by disinfection (control lever) [illustrative] (per year)
  beta_hs = 1e-8,  # sheep-to-human transmission coefficient (occupational contact) [illustrative] (per sheep per year)
  beta_hw = 1e-8,  # environment-to-human transmission coefficient [illustrative] (per unit load per year)
  sigma_h = 4.0,   # rate exposed → acute human case (≈0.25 year onset) [verified as a published value for this quantity in the related literature] (per year)
  chron = 1.0,     # rate acute → chronic human infection [illustrative] (per year)
  gamma_a = 2.0,   # recovery rate from acute infection [illustrative] (per year)
  gamma_c = 0.4,   # recovery rate from chronic infection [illustrative] (per year)
  mu_h = 0.006     # human natural death rate [verified as a published value in the related literature] (per year)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Ss = 1e6,  # susceptible sheep
  Es = 1e4,  # exposed/latently infected sheep
  Is = 1e3,  # infectious sheep (shedding; abortion/parturition products)
  Vs = 0,    # vaccinated sheep
  W = 1e3,   # Brucella in the environment (bacterial load, arbitrary units)
  Sh = 2e6,  # susceptible humans
  Eh = 0,    # exposed humans
  Iah = 0,   # acutely infected humans
  Ich = 0    # chronically infected humans
)

# ---- Model equations -------------------------------------------------
brucellosis_sheep_human_environment_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N_total <- Sh + Eh + Iah + Ich

    # Derived quantities
    Ns <- Ss + Es + Is + Vs                                                                                                             # sheep population
    N <- Sh + Eh + Iah + Ich                                                                                                            # human population
    lambda_s <- beta_ss * Is + beta_sw * W                                                                                              # sheep force of infection
    lambda_h <- beta_hs * Is + beta_hw * W                                                                                              # human force of infection
    R0s <- (beta_ss + beta_sw * shed / (decay + disinf)) * As * sigma_s / (mu_s * (mu_s + sigma_s + nu_s) * (mu_s + shed * 0 + 1e-12))  # sheep reproduction number (write out explicitly in the generator; the point is that R0s has a direct and an indirect additive component, whi
    CHh_rate <- (sigma_h) * Eh                                                                                                          # rate of new cumulative reported human cases (integrate with cumsum for the cumulative count)

    # Rates of change
    dSs <- As - lambda_s * Ss - nu_s * Ss + wane_s * Vs - mu_s * Ss
    dEs <- lambda_s * Ss - sigma_s * Es
    dIs <- sigma_s * Es
    dVs <- nu_s * Ss - wane_s * Vs
    dW <- shed * Is - (decay + disinf) * W
    dSh <- -lambda_h * Sh + gamma_a * Iah + gamma_c * Ich - mu_h * Sh
    dEh <- lambda_h * Sh - sigma_h * Eh
    dIah <- sigma_h * Eh - chron * Iah - gamma_a * Iah
    dIch <- chron * Iah - gamma_c * Ich

    list(c(dSs, dEs, dIs, dVs, dW, dSh, dEh, dIah, dIch), c(CHh_rate = CHh_rate))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = brucellosis_sheep_human_environment_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Ss = "#333333", Es = "#E8C468", Is = "#B71C1C", Vs = "#2457C5", W = "#A5D6A7", Sh = "#333333", Eh = "#E8C468", Iah = "#B71C1C", Ich = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Brucellosis in small ruminants with an environmental reservoir and human cases")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")

output_names <- c("CHh_rate")
matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,
        col = seq_along(output_names), xlab = "Time (years)", ylab = "Value",
        main = "Brucellosis in small ruminants with an environmental reservoir and human cases: outputs")
legend("topright", legend = output_names, col = seq_along(output_names),
       lty = 1, lwd = 2, bty = "n")
