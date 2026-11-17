## Visceral leishmaniasis (anthroponotic, Indian subcontinent) with asymptomatic, VL, treatment and PKDL; sandfly SEI (simplified after Stauch et al. 2011) [ILLUSTRATIVE-COMPOSITE]
## Visceral leishmaniasis (L. donovani, Phlebotomus argentipes), Bihar, India / Nepal.
## Parameter values: 0 of 17 checked against the source; the rest are recalled or illustrative.
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sh, Ah, Kh, Tx, Ph, Rh, Sf, Ef, If
## Flows (time unit: days):
##   (outside) -> Sh  : mu_h * Nh
##   Sh -> Ah         : lambda_h * Sh
##   Ah -> Kh         : p_K * w_A * Ah
##   Ah -> Rh         : (1 - p_K) * w_A * Ah
##   Kh -> Tx         : r_tx * Kh
##   Kh -> (outside)  : d_K * Kh
##   Tx -> Ph         : f_P * w_Tx * Tx
##   Tx -> Rh         : (1 - f_P) * w_Tx * Tx
##   Ph -> Rh         : w_P * Ph
##   Rh -> Sh         : rho_h * Rh
##   Sh -> (outside)  : mu_h * Sh
##   Ah -> (outside)  : mu_h * Ah
##   Kh -> (outside)  : mu_h * Kh
##   Tx -> (outside)  : mu_h * Tx
##   Ph -> (outside)  : mu_h * Ph
##   Rh -> (outside)  : mu_h * Rh
##   (outside) -> Sf  : mu_f * Nf
##   Sf -> Ef         : lambda_f * Sf
##   Ef -> If         : nu_f * Ef
##   Sf -> (outside)  : mu_f * Sf
##   Ef -> (outside)  : mu_f * Ef
##   If -> (outside)  : mu_f * If

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 7300, by = 10)  # days

# ---- Parameters ------------------------------------------------------
parameters <- c(
  a_sf = 0.25,       # sandfly biting rate (bite every 4 days) (1/day)
  p_fh = 0.1,        # sandfly → human per-bite infection probability
  k_A = 0.025,       # per-bite infectivity of asymptomatic humans to sandflies
  k_K = 0.5,         # per-bite infectivity of VL patients
  k_Tx = 0.25,       # per-bite infectivity of treated patients
  k_P = 0.3,         # per-bite infectivity of PKDL patients
  w_A = 0.00667,     # exit from asymptomatic (150 d) (1/day)
  p_K = 0.1,         # fraction of infections progressing to VL
  r_tx = 0.0167,     # VL onset → treatment (60 d) (1/day)
  d_K = 0.00167,     # VL death before treatment (1/day)
  w_Tx = 0.0357,     # treatment duration (28 d) (1/day)
  f_P = 0.025,       # fraction developing PKDL after treatment
  w_P = 0.00091,     # PKDL duration (3 yr) (1/day)
  rho_h = 0.000548,  # waning immunity (5 yr) (1/day)
  mu_h = 4.57e-5,    # human birth/death (60-yr life) (1/day)
  mu_f = 0.0714,     # sandfly death (14-d life) (1/day)
  nu_f = 0.2         # 1/EIP (5 d) (1/day)
)

# ---- Initial conditions ----------------------------------------------
state <- c(
  Sh = 9800,   # susceptible humans
  Ah = 150,    # asymptomatic infected (weakly infectious)
  Kh = 10,     # clinical VL (kala-azar, infectious)
  Tx = 0,      # under treatment
  Ph = 0,      # post-kala-azar dermal leishmaniasis (infectious)
  Rh = 40,     # recovered/immune
  Sf = 49000,  # susceptible sandflies
  Ef = 500,    # sandflies in extrinsic incubation
  If = 500     # infectious sandflies
)

# ---- Model equations -------------------------------------------------
vl_leishmaniasis_isc_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sh + Ah + Kh + Tx + Ph + Rh

    # Derived quantities
    Nh <- Sh + Ah + Kh + Tx + Ph + Rh
    Nf <- Sf + Ef + If
    lambda_h <- a_sf * p_fh * If / Nh                                     # human force of infection
    lambda_f <- a_sf * (k_A * Ah + k_K * Kh + k_Tx * Tx + k_P * Ph) / Nh  # sandfly force of infection

    # Rates of change
    dSh <- mu_h * Nh - lambda_h * Sh + rho_h * Rh - mu_h * Sh
    dAh <- lambda_h * Sh - p_K * w_A * Ah - (1 - p_K) * w_A * Ah - mu_h * Ah
    dKh <- p_K * w_A * Ah - r_tx * Kh - d_K * Kh - mu_h * Kh
    dTx <- r_tx * Kh - f_P * w_Tx * Tx - (1 - f_P) * w_Tx * Tx - mu_h * Tx
    dPh <- f_P * w_Tx * Tx - w_P * Ph - mu_h * Ph
    dRh <- (1 - p_K) * w_A * Ah + (1 - f_P) * w_Tx * Tx + w_P * Ph - rho_h * Rh - mu_h * Rh
    dSf <- mu_f * Nf - lambda_f * Sf - mu_f * Sf
    dEf <- lambda_f * Sf - nu_f * Ef - mu_f * Ef
    dIf <- nu_f * Ef - mu_f * If

    list(c(dSh, dAh, dKh, dTx, dPh, dRh, dSf, dEf, dIf))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = vl_leishmaniasis_isc_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sh = "#333333", Ah = "#E57373", Kh = "#333333", Tx = "#333333", Ph = "#F0A35E", Rh = "#4E8A3E", Sf = "#333333", Ef = "#E8C468", If = "#B71C1C")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (days)", ylab = "Number of individuals",
        main = "Visceral leishmaniasis (anthroponotic, Indian subcontinent) with asymptomatic, VL, treatment and PKDL; sandfly SEI (simplified after Stauch et al. 2011) [ILLUSTRATIVE-COMPOSITE]")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
