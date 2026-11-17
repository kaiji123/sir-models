## African trypanosomiasis: humans + cattle reservoir + tsetse with two-stage incubation (after Rogers 1988) [ILLUSTRATIVE-COMPOSITE]
## Rhodesian sleeping sickness / animal trypanosomiasis (T. brucei rhodesiense, Glossina spp.), East Africa.
## Parameter values: 0 of 17 checked against the source; the rest are recalled or illustrative.
## Reference: Rogers DJ (1988). A general model for the African trypanosomiases. Parasitology 97(1):193–212. https://doi.org/10.1017/S0031182000066853 . It is a multi-host, single-vector Ross–Macdonald generalisation with host-specific feeding fractions and…
## Source: https://doi.org/10.1017/S0031182000066853
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sh, Eh, Ih, Rh, Sc, Ec, Ic, Rc, Sf, Ef1, Ef2, If
## Flows (time unit: days):
##   (outside) -> Sh   : mu_h * Nh
##   Sh -> Eh          : lam_h * Sh
##   Eh -> Ih          : nu_h * Eh
##   Ih -> Rh          : gam_h * Ih
##   Ih -> (outside)   : d_h * Ih
##   Rh -> Sh          : rho_h * Rh
##   Sh -> (outside)   : mu_h * Sh
##   Eh -> (outside)   : mu_h * Eh
##   Ih -> (outside)   : mu_h * Ih
##   Rh -> (outside)   : mu_h * Rh
##   (outside) -> Sc   : mu_c * Nc
##   Sc -> Ec          : lam_c * Sc
##   Ec -> Ic          : nu_c * Ec
##   Ic -> Rc          : gam_c * Ic
##   Rc -> Sc          : rho_c * Rc
##   Sc -> (outside)   : mu_c * Sc
##   Ec -> (outside)   : mu_c * Ec
##   Ic -> (outside)   : mu_c * Ic
##   Rc -> (outside)   : mu_c * Rc
##   (outside) -> Sf   : mu_f * Nf
##   Sf -> Ef1         : lam_f * Sf
##   Ef1 -> Ef2        : 2 / eip * Ef1
##   Ef2 -> If         : 2 / eip * Ef2
##   Sf -> (outside)   : mu_f * Sf
##   Ef1 -> (outside)  : mu_f * Ef1
##   Ef2 -> (outside)  : mu_f * Ef2
##   If -> (outside)   : mu_f * If

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 7300, by = 10)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  a_f = 0.25,      # feeding rate (every 4 d) (1/day)
  f_hum = 0.1,     # fraction of meals on humans; the rest on cattle
  b_fh = 0.6,      # fly → human per-bite probability
  b_fc = 0.6,      # fly → cattle per-bite probability
  p_hf = 0.05,     # human → fly per-meal probability (low because mainly teneral flies are susceptible)
  p_cf = 0.05,     # cattle → fly per-meal probability
  mu_f = 0.03,     # tsetse mortality (1/day)
  eip = 25,        # EIP of T. brucei (two Erlang stages at rate 2/eip each) (days)
  nu_h = 0.0714,   # 1/14 d human incubation (1/day)
  gam_h = 0.0083,  # 1/120 d to treatment/recovery (1/day)
  d_h = 0.0017,    # disease death (1/day)
  rho_h = 0.0167,  # loss of protection after treatment (1/day)
  mu_h = 4.57e-5,  # human birth/death (1/day)
  nu_c = 0.1,      # cattle incubation (1/day)
  gam_c = 0.01,    # cattle recovery (100 d infectious) (1/day)
  rho_c = 0.02,    # cattle waning (1/day)
  mu_c = 3.65e-4   # cattle birth/death (7.5 yr) (1/day)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sh = 995,   # susceptible humans
  Eh = 0,     # incubating humans
  Ih = 5,     # infectious humans
  Rh = 0,     # treated/recovered (temporary)
  Sc = 490,   # susceptible cattle
  Ec = 0,     # incubating cattle
  Ic = 10,    # infectious cattle
  Rc = 0,     # recovered cattle
  Sf = 5000,  # susceptible tsetse
  Ef1 = 0,    # tsetse in EIP stage 1
  Ef2 = 0,    # tsetse in EIP stage 2
  If = 0      # infectious tsetse
)

# ---- Model equations -------------------------------------------------
trypanosomiasis_rogers1988_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sh + Eh + Ih + Rh

    # Derived quantities
    Nh <- Sh + Eh + Ih + Rh
    Nc <- Sc + Ec + Ic + Rc
    Nf <- Sf + Ef1 + Ef2 + If
    lam_h <- a_f * f_hum * b_fh * If / Nh                                   # human force of infection
    lam_c <- a_f * (1 - f_hum) * b_fc * If / Nc                             # cattle force of infection
    lam_f <- a_f * (f_hum * p_hf * Ih / Nh + (1 - f_hum) * p_cf * Ic / Nc)  # tsetse force of infection

    # Rates of change
    dSh <- mu_h * Nh - lam_h * Sh + rho_h * Rh - mu_h * Sh
    dEh <- lam_h * Sh - nu_h * Eh - mu_h * Eh
    dIh <- nu_h * Eh - gam_h * Ih - d_h * Ih - mu_h * Ih
    dRh <- gam_h * Ih - rho_h * Rh - mu_h * Rh
    dSc <- mu_c * Nc - lam_c * Sc + rho_c * Rc - mu_c * Sc
    dEc <- lam_c * Sc - nu_c * Ec - mu_c * Ec
    dIc <- nu_c * Ec - gam_c * Ic - mu_c * Ic
    dRc <- gam_c * Ic - rho_c * Rc - mu_c * Rc
    dSf <- mu_f * Nf - lam_f * Sf - mu_f * Sf
    dEf1 <- lam_f * Sf - 2 / eip * Ef1 - mu_f * Ef1
    dEf2 <- 2 / eip * Ef1 - 2 / eip * Ef2 - mu_f * Ef2
    dIf <- 2 / eip * Ef2 - mu_f * If

    list(c(dSh, dEh, dIh, dRh, dSc, dEc, dIc, dRc, dSf, dEf1, dEf2, dIf))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = trypanosomiasis_rogers1988_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sh = "#333333", Eh = "#E8C468", Ih = "#B71C1C", Rh = "#4E8A3E", Sc = "#333333", Ec = "#E8C468", Ic = "#B71C1C", Rc = "#4E8A3E", Sf = "#333333", Ef1 = "#E8C468", Ef2 = "#E8C468", If = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "African trypanosomiasis: humans + cattle reservoir + tsetse with two-stage incubation (after Rogers 1988) [ILLUSTRATIVE-COMPOSITE]")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
