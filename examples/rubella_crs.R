## Rubella with two age groups, routine infant vaccination, and congenital rubella syndrome (CRS) as an output
## Rubella / CRS.
## Parameter values: 0 of 9 checked against the source; the rest are recalled or illustrative.
##
## Generated from a compartmental model diagram (model format v1).
## Regenerating from the same diagram always gives this exact file.
##
## Compartments: Sch, Ich, Rch, Sad, Iad, Rad, V, CRS
## Flows (time unit: years):
##   (outside) -> Sch  : (1 - p_vac) * births
##   (outside) -> V    : p_vac * births
##   Sch -> Ich        : foi * Sch
##   Ich -> Rch        : g_rec * Ich
##   Sad -> Iad        : foi * Sad
##   Iad -> Rad        : g_rec * Iad
##   Sch -> Sad        : age_out * Sch
##   Ich -> Iad        : age_out * Ich
##   Rch -> Rad        : age_out * Rch
##   Sad -> (outside)  : mu_ad * Sad
##   V -> (outside)    : 1 / 70 * V
##   (outside) -> CRS  : crs_rate

library(deSolve)

# ---- Time ------------------------------------------------------------
times <- seq(0, 80, by = 0.1)  # years

# ---- Parameters ------------------------------------------------------
parameters <- c(
  R0_r = 7,           # illustrative; commonly cited range about 6–7
  g_rec = 365 / 11,   # 11-day infectious period (1/year) [illustrative]
  age_out = 1 / 15,   # ageing from child to adult (1/year)
  mu_ad = 1 / 55,     # adult death rate, so life expectancy is 70 (1/year) [illustrative]
  p_vac = 0.8,        # fraction of births vaccinated [illustrative]
  preg_frac = 0.015,  # fraction of adults who are women in the first 16 weeks of pregnancy at any time [illustrative]
  q_crs = 0.6,        # probability of CRS given infection in that window, averaged over weeks 0–16 [illustrative]
  births = 1e6 / 70,  # constant births, roughly balancing deaths (persons/year) [illustrative]
  N0 = 1e6            # value taken from the notes text [recalled]
)

# ---- Initial conditions ----------------------------------------------
state <- with(as.list(parameters), {
  Sch <- 0.1 * N0   # susceptible child, 0–14 y
  Ich <- 50         # infectious child
  Rch <- 0.1 * N0   # immune child
  Sad <- 0.05 * N0  # susceptible adult, 15+ y
  Iad <- 50         # infectious adult
  Rad <- 0.6 * N0   # immune adult
  V <- 0.15 * N0    # vaccinated, all ages
  CRS <- 0          # cumulative CRS births
  c(Sch = Sch, Ich = Ich, Rch = Rch, Sad = Sad, Iad = Iad, Rad = Rad, V = V, CRS = CRS)
})

# ---- Model equations -------------------------------------------------
rubella_crs_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Total population
    N <- Sch + Ich + Rch + Sad + Iad + Rad + V

    # Derived quantities
    foi <- R0_r * g_rec * (Ich + Iad) / N
    crs_rate <- q_crs * preg_frac * foi * Sad

    # Rates of change
    dSch <- (1 - p_vac) * births - foi * Sch - age_out * Sch
    dIch <- foi * Sch - g_rec * Ich - age_out * Ich
    dRch <- g_rec * Ich - age_out * Rch
    dSad <- -foi * Sad + age_out * Sch - mu_ad * Sad
    dIad <- foi * Sad - g_rec * Iad + age_out * Ich
    dRad <- g_rec * Iad + age_out * Rch
    dV <- p_vac * births - 1 / 70 * V
    dCRS <- crs_rate

    list(c(dSch, dIch, dRch, dSad, dIad, dRad, dV, dCRS))
  })
}

# ---- Solve -----------------------------------------------------------
out <- ode(y = state, times = times, func = rubella_crs_model, parms = parameters)
out <- as.data.frame(out)

# ---- Plot ------------------------------------------------------------
compartment_cols <- c(Sch = "#333333", Ich = "#B71C1C", Rch = "#4E8A3E", Sad = "#333333", Iad = "#B71C1C", Rad = "#4E8A3E", V = "#2457C5", CRS = "#8E0E0E")
matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,
        col = compartment_cols, xlab = "Time (years)", ylab = "Number of individuals",
        main = "Rubella with two age groups, routine infant vaccination, and congenital rubella syndrome (CRS) as an output")
legend("topright", legend = names(compartment_cols), col = compartment_cols,
       lty = 1, lwd = 2, bty = "n")
