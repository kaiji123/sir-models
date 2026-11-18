## TB with treatment, acquired and transmitted MDR, and a fitness cost (after Dye & Williams 2000 and Cohen & Murray 2004)
## Tuberculosis.
## Parameter values: 0 of 14 checked against the source; the rest are recalled or illustrative.
## Reference: Dye C, Williams BG. "Criteria for the control of drug-resistant tuberculosis." PNAS 2000;97:8180–8185
## Source: https://doi.org/10.1073/pnas.140102797
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: S, Ls, Lr, Is, Ir, Rt
## Flows (time unit: years):
##   (outside) -> S   : (mu * N + mu_tb * (Is + Ir))
##   S -> Is          : p_fast * foi_s * S
##   S -> Ls          : (1 - p_fast) * foi_s * S
##   S -> Ir          : p_fast * foi_r * S
##   S -> Lr          : (1 - p_fast) * foi_r * S
##   Ls -> Is         : (v_react + p_fast * (1 - x_prot) * foi_s) * Ls
##   Ls -> Ir         : p_fast * (1 - x_prot) * foi_r * Ls
##   Lr -> Ir         : (v_react + p_fast * (1 - x_prot) * foi_r) * Lr
##   Lr -> Is         : p_fast * (1 - x_prot) * foi_s * Lr
##   Rt -> Is         : (relapse + p_fast * (1 - x_prot) * foi_s) * Rt
##   Rt -> Ir         : p_fast * (1 - x_prot) * foi_r * Rt
##   Is -> Rt         : (cure_nat + det * succ_s) * Is
##   Is -> Ir         : det * a_res * Is
##   Ir -> Rt         : (cure_nat + det * succ_r) * Ir
##   Is -> (outside)  : (mu + mu_tb) * Is

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 200, by = 0.2)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  beta_s = 10,       # effective contacts per infectious case per year (1/year) [illustrative]
  fit_r = 0.8,       # relative transmissibility of MDR
  p_fast = 0.14,     # fraction who progress quickly after infection [illustrative]
  v_react = 0.0005,  # reactivation from latency (1/year) [illustrative]
  x_prot = 0.65,     # protection of latent and recovered people against reinfection-driven fast progression [illustrative]
  cure_nat = 0.2,    # self-cure (1/year) [illustrative]
  mu_tb = 0.2,       # TB death, untreated (1/year) [illustrative]
  det = 0.7,         # detection and treatment rate (1/year) [illustrative]
  succ_s = 0.85,     # treatment success for DS-TB [illustrative]
  a_res = 0.03,      # fraction of DS treatments that produce MDR [illustrative]
  succ_r = 0.5,      # first-line treatment success for MDR [illustrative]
  relapse = 0.002,   # relapse from Rt to DS active disease (1/year) [illustrative]
  mu = 0.02,         # birth = non-TB death rate (1/year)
  N0 = 1e6           # value taken from the notes text [recalled]
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  S <- 0.7 * N0     # uninfected
  Ls <- 0.29 * N0   # latent, drug-sensitive
  Lr <- 0           # latent, MDR
  Is <- 0.002 * N0  # active TB, drug-sensitive
  Ir <- 10          # active TB, MDR
  Rt <- 0           # recovered/treated, latent infection cleared or controlled
  c(S = S, Ls = Ls, Lr = Lr, Is = Is, Ir = Ir, Rt = Rt)
})

# ---- Model equations -------------------------------------------------
tb_mdr_treatment_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- S + Ls + Lr + Is + Ir + Rt

    # Derived quantities
    foi_s <- beta_s * Is / N
    foi_r <- beta_s * fit_r * Ir / N

    # Rates of change
    dS <- (mu * N + mu_tb * (Is + Ir)) - p_fast * foi_s * S - (1 - p_fast) * foi_s * S - p_fast * foi_r * S - (1 - p_fast) * foi_r * S
    dLs <- (1 - p_fast) * foi_s * S - (v_react + p_fast * (1 - x_prot) * foi_s) * Ls - p_fast * (1 - x_prot) * foi_r * Ls
    dLr <- (1 - p_fast) * foi_r * S - (v_react + p_fast * (1 - x_prot) * foi_r) * Lr - p_fast * (1 - x_prot) * foi_s * Lr
    dIs <- p_fast * foi_s * S + (v_react + p_fast * (1 - x_prot) * foi_s) * Ls + p_fast * (1 - x_prot) * foi_s * Lr + (relapse + p_fast * (1 - x_prot) * foi_s) * Rt - (cure_nat + det * succ_s) * Is - det * a_res * Is - (mu + mu_tb) * Is
    dIr <- p_fast * foi_r * S + p_fast * (1 - x_prot) * foi_r * Ls + (v_react + p_fast * (1 - x_prot) * foi_r) * Lr + p_fast * (1 - x_prot) * foi_r * Rt + det * a_res * Is - (cure_nat + det * succ_r) * Ir
    dRt <- -(relapse + p_fast * (1 - x_prot) * foi_s) * Rt - p_fast * (1 - x_prot) * foi_r * Rt + (cure_nat + det * succ_s) * Is + (cure_nat + det * succ_r) * Ir

    list(c(dS, dLs, dLr, dIs, dIr, dRt))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = tb_mdr_treatment_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(S = "#333333", Ls = "#E8C468", Lr = "#E8C468", Is = "#B71C1C", Ir = "#B71C1C", Rt = "#4E8A3E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "TB with treatment, acquired and transmitted MDR, and a fitness cost (after Dye & Williams 2000 and Cohen & Murray 2004)")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
