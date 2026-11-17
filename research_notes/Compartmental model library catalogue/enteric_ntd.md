# Water-borne, enteric, environmental and NTD compartmental models (deterministic ODE specs)

Verification legend, used throughout:
- **[V]**: verified this session against the primary paper or a fetched secondary source (link given).
- **[M]**: recalled from the literature, not re-verified this session. Check against the source before treating it as the published value.
- **[I]**: illustrative. Chosen by the researcher so the model runs and gives plausible dynamics; it is not a published value.

Identifier conventions: the specs avoid R function and reserved names such as `beta`, `gamma`, `sigma`, `kappa`, `cov`, `scale`, `q`, `D`, `C`, `c`, `F`, `T` and `pi`, and write 2π as `6.2831853`. The compartment name `I` is kept for consistency with the existing SIWR model, but `I()` is a base R function. If the generator rejects it, rename it to `Inf_h`, not `Inf`, which is also reserved.

Not duplicated, as instructed: Tien & Earn SIWR cholera, the HEV SEIR-latrine model, and the generic SICR carrier model.

---

## For each model: compartments, flows with rate expressions, force of infection (including dose–response), and parameter values with units

### Takeaway
There are 13 specifications: cholera ×2, typhoid, polio (WPV+OPV/cVDPV), rotavirus, norovirus, hepatitis A, schistosomiasis, soil-transmitted helminths (STH), trachoma with MDA, cryptosporidiosis, leptospirosis and scabies with MDA. All fit the ≤14-compartment FROM > TO format. Only Codeço 2001 and parts of Simmons 2013 and Pitzer 2014 were checked number-by-number against the papers. Most other values are [M] or [I] and need checking before they are presented as published values.

### Cited Findings (model specifications)

### cholera_codeco2001 — Codeço endemic/epidemic cholera with aquatic reservoir and dose–response
- Disease/area: Cholera (Vibrio cholerae) | Category: Water, enteric & NTDs
- Source: Codeço CT. Endemic and epidemic dynamics of cholera: the role of the aquatic reservoir. *BMC Infect Dis* 2001;1:1. doi:10.1186/1471-2334-1-1 (PMID 11208258, PMC29087). Parameter values come from the paper's Table 2, "Community 2" [V] — [Europe PMC full text](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC29087/fullTextXML)
- Time unit: days; simulate 365–730 days for an epidemic, or 20+ years with seasonality in `a` for endemic dynamics.
- Compartments:
  - `S` — susceptible — 9999
  - `I` — infected — 1
  - `R` — recovered/removed. Added for bookkeeping: the original has only S and I with constant total H, and recovered people are implicitly not re-susceptible within the timescale — 0
  - `B` — toxigenic V. cholerae concentration in water (cells/ml) — 0 [initial value I] — not in N
- Parameters:
  - `mu_h` = 0.0001 (day^-1) — human birth = death rate (the paper's n) — Table 2 [V]. Community 3 uses 0.001.
  - `a` = 1 (day^-1) — rate of exposure to contaminated water — Table 2 Community 2 [V]. Community 1 uses 0.5.
  - `K` = 1e6 (cells/ml) — concentration giving a 50% chance of infection — [V]
  - `r` = 0.2 (day^-1) — recovery rate (5-day infection) — [V]
  - `nb_mb` = -0.33 (day^-1) — net growth rate of vibrios, nb − mb. It is negative, so the population cannot persist without human input — [V]
  - `e` = 10 (cells ml^-1 day^-1 person^-1) — contribution of each infected person to the aquatic population — [V]
- Derived variables:
  - `lam_B` = B/(K + B) — probability of infection on exposure (saturating dose–response)
  - `S_crit` = K*r*(-nb_mb)/(e*a) — critical susceptible number for invasion: 6600 for Community 2 and 13,200 for Community 1. These match the paper's Table 2 S_C values, which confirms the B/(K+B) form and the parameter reading [V by calculation].
- Flows:
  - ` > S : mu_h*N` (absolute; births)
  - `S > : mu_h`
  - `S > I : a*B/(K + B)`
  - `I > R : r`
  - `I > : mu_h`
  - `R > : mu_h`
  - ` > B : e*I` (absolute; shedding into water)
  - `B > : -nb_mb` (per capita; net loss 0.33/day)
- Notes:
  - The paper's Table 1 describes infection probability in words as a "logistic dose response". The algebraic form B/(K+B), with K the 50% infection dose, reproduces the published S_C values.
  - With Community 1 (a=0.5), N=10,000 is below S_C, so no epidemic occurs. Use a=1 to get an outbreak.
  - For seasonality, replace `a` with `a*(1 + amp*cos(6.2831853*t/365))`. This is an [I] extension.

### cholera_andrews_basu2011 — Haiti cholera with water reservoir, vaccination and water-sanitation intervention
- Disease/area: Cholera, Haiti 2010–11 | Category: Water, enteric & NTDs
- Source: Andrews JR, Basu S. Transmission dynamics and control of cholera in Haiti: an epidemic model. *Lancet* 2011;377(9773):1248–1255. doi:10.1016/S0140-6736(11)60273-0. PMC3172163 — [PMC record](https://pmc.ncbi.nlm.nih.gov/articles/PMC3172163)
  - The model was adapted from earlier cholera models. Fitted parameters were the contaminated-water consumption rate, the symptomatic mortality rate and the initial environmental vibrio concentration, using MCMC. Pathogenesis parameters came from the literature. Main results: 779,000 cases and 11,100 deaths projected for March–November 2011, and a 1%/week reduction in contaminated-water consumption averting 105,000 cases [V, abstract-level].
- Time unit: days; simulate 300 days.
- Compartments (structure reconstructed in the spirit of the paper [M/I]; the full text could not be fetched this session):
  - `S` — susceptible — 9,999,000 [I]
  - `I` — symptomatic infected — 1000 [I]
  - `A` — asymptomatic infected — 0
  - `R` — recovered/immune — 0
  - `V` — vaccinated (OCV) — 0
  - `W` — V. cholerae concentration in the water reservoir (cells/ml) — 1e3 [I] — not in N
- Parameters (all [I] or [M]; calibrate):
  - `beta_w` = 1.0 (day^-1) — contaminated-water exposure rate (fitted in the paper) [I]
  - `beta_h` = 0.1 (day^-1) — direct (short-cycle) transmission [I]
  - `kap_w` = 1e6 (cells/ml) — half-saturation dose, as in Codeço [M]
  - `p_sym` = 0.2 (–) — fraction of infections that are symptomatic [I]
  - `rel_a` = 0.1 (–) — relative shedding/infectiousness of asymptomatic infections [I]
  - `gam_c` = 0.2 (day^-1) — recovery rate (5 days) [M]
  - `cfr_rate` = 0.002 (day^-1) — symptomatic death rate (fitted in the paper) [I]
  - `xi_w` = 10 (cells ml^-1 day^-1 per symptomatic person) — shedding [I]
  - `del_w` = 0.033 (day^-1) — vibrio decay rate (30-day survival) [M, as used in Hartley/Morris/Smith 2006-style models]
  - `omega_r` = 1/1095 (day^-1) — waning of natural immunity (about 3 years) [M]
  - `nu_v` = 0 (day^-1) — vaccination rate; try 0.001 [I]
  - `eff_v` = 0.65 (–) — vaccine efficacy [I]
  - `omega_v` = 1/1095 (day^-1) — vaccine waning [I]
  - `red_w` = 0.001436 (day^-1) — equals −log(0.99)/7, a 1%/week reduction in contaminated-water consumption [V for the 1%/week scenario]
  - `t_int` = 120 (day) — start of the water intervention [I]
- Derived variables:
  - `water_mult` = exp(-red_w*max(t - t_int, 0))
  - `lam` = beta_w*water_mult*W/(kap_w + W) + beta_h*(I + rel_a*A)/N
- Flows:
  - `S > I : p_sym*lam`
  - `S > A : (1 - p_sym)*lam`
  - `I > R : gam_c`
  - `I > : cfr_rate`
  - `A > R : gam_c`
  - `R > S : omega_r`
  - `S > V : nu_v*eff_v`
  - `V > S : omega_v`
  - ` > W : xi_w*(I + rel_a*A)` (absolute)
  - `W > : del_w`
- Notes:
  - The published model was fitted per department, in 10 sub-populations. This is a single-patch simplification.
  - Treat this as a "Haiti-style" model. The exact published compartment list, priors and posteriors (paper appendix) were not verified.

### typhoid_pitzer2014 — Typhoid with chronic carriers, short- and long-cycle (water) transmission
- Disease/area: Typhoid fever (S. Typhi), South Asia (Kolkata, Delhi, Vellore) | Category: Water, enteric & NTDs
- Source: Pitzer VE, Bowles CC, Baker S, Kang G, Balaji V, Farrar JJ, Grenfell BT. Predicting the impact of vaccination on the transmission dynamics of typhoid in South Asia: a mathematical modeling study. *PLoS Negl Trop Dis* 2014;8(1):e2642. doi:10.1371/journal.pntd.0002642 — [PLoS NTD](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0002642). The equations are in Text S1 and Table 1 (an image), which were not readable this session.
- Time unit: weeks; simulate 52×50 weeks to reach equilibrium.
- Compartments:
  - `S1` — fully susceptible — 0.3*N0 [I]
  - `I1` — primary (clinical-capable) infection — 100 [I]
  - `Car` — chronic carrier. `C` is avoided because it is an R function — 0.005*N0 [I]
  - `R` — immune — 0.4*N0 [I]
  - `S2` — partially susceptible after waning — 0.3*N0 [I]
  - `I2` — subclinical reinfection — 0
  - `W` — water-supply contamination, scaled to prevalence units — 0 — not in N
- Parameters:
  - `R0_ty` = 2.8 (–) — basic reproduction number, Kolkata fit — [V]
  - `frac_short` = 0.9 (–) — share of R0 from short-cycle transmission — [V]
  - `r_car` = 0.01 (–) — relative infectiousness of carriers (best fit). The model also applies it to I2 [V for carriers; M/I for applying it to I2].
  - `theta_c` = 0.03 (–) — probability that a primary infection becomes a lifelong carrier. The paper uses age-specific values of 0.3–10.1%, so 0.03 is an [I] average [V range].
  - `p_clin` = 0.1 (–) — fraction of primary infections that are clinical typhoid. Used only for reporting [V].
  - `rep_f` = 0.0052 (–) — reporting fraction — [V]
  - `delta_i` = 0.25 (week^-1) — recovery rate (4-week infectiousness) [M]
  - `omega_r` = 0.0096 (week^-1) — waning R → S2 (about 104 weeks) [M]
  - `xi_w` = 0.333 (week^-1) — decay of water contamination (about 3 weeks) [M]
  - `mu_h` = 0.000385 (week^-1) — birth = death rate (50-year life) [I]. The paper used observed demography.
  - `amp_w` = 0.3 (–) — seasonal amplitude of long-cycle transmission [I]. The paper fitted seasonality.
- Derived variables:
  - `D_eff` = 1/(delta_i + mu_h) + theta_c*delta_i/(delta_i + mu_h)*r_car/mu_h — effective infectious duration per primary case (about 4.78 weeks)
  - `beta_s` = frac_short*R0_ty/D_eff — short-cycle rate (week^-1)
  - `beta_l` = (1 - frac_short)*R0_ty/D_eff — long-cycle rate. W is scaled so that its equilibrium equals the weighted prevalence.
  - `prev_w` = (I1 + r_car*Car + r_car*I2)/N
  - `lam` = beta_s*prev_w + beta_l*(1 + amp_w*cos(6.2831853*t/52))*W
- Flows:
  - ` > S1 : mu_h*N` (absolute)
  - `S1 > I1 : lam`
  - `I1 > Car : theta_c*delta_i`
  - `I1 > R : (1 - theta_c)*delta_i`
  - `R > S2 : omega_r`
  - `S2 > I2 : lam`
  - `I2 > R : delta_i`
  - `S1 > : mu_h`
  - `I1 > : mu_h`
  - `Car > : mu_h`
  - `R > : mu_h`
  - `S2 > : mu_h`
  - `I2 > : mu_h`
  - ` > W : xi_w*prev_w` (absolute; W relaxes to prev_w at rate xi_w)
  - `W > : xi_w`
- Notes:
  - The paper estimated immunity to clinical disease after recovery at ε ≈ 0 [V]. A primary-infection class can therefore be followed by subclinical reinfections.
  - Vaccines in the paper [V]: Ty21a efficacy 48%; ViPS 80% efficacy with 3-year duration; Vi-conjugate (ViCV) 95.6% efficacy with 19.2-year duration; 80% coverage. A vaccine extension adds `V` with `S1 > V : nu_v*eff_v` and `V > S1 : 1/duration`.
  - Sanitation was modelled by reducing R0_w linearly to 0 over years 5–35 [V]. Implement this by multiplying `beta_l` by `max(0, 1 - max(t - 260, 0)/1560)`.
  - Predicted incidence was 146 per 100,000 person-years [V].
  - The age structure is dropped in this spec.

### polio_wpv_opv — Wild poliovirus with OPV vaccine virus and circulating vaccine-derived (cVDPV) lineage
- Disease/area: Poliomyelitis (WPV, Sabin OPV, cVDPV) | Category: Water, enteric & NTDs
- Source: Structure inspired by Fine PEM, Carneiro IAM. Transmissibility and persistence of oral polio vaccine viruses: implications for the global poliomyelitis eradication initiative. *Am J Epidemiol* 1999;150(10):1001–1021. doi:10.1093/oxfordjournals.aje.a009924 [M citation]. The simplified OPV/cVDPV form follows Pons-Salort M, Grassly NC et al. (e.g. *PLoS Pathog* 2016 / *Science* 2020 work on Sabin2 cessation) [M]. Exact equations and values were not verified this session; all values are [I] unless stated.
- Time unit: days; simulate 20 years (7300 days).
- Compartments:
  - `S` — susceptible to intestinal infection — 0.1*N0
  - `IW` — wild-type infection — 10
  - `IV` — OPV (Sabin-like) infection — 0
  - `IVD` — cVDPV (reverted) infection — 0
  - `R` — intestinally immune — 0.9*N0
- Parameters:
  - `mu_h` = 1/(50*365) (day^-1) — birth/death [I]
  - `gam_p` = 1/28 (day^-1) — recovery (about 4 weeks of excretion) [M]
  - `R0_w` = 6 (–) — WPV R0; reported range about 4–13 in high-transmission settings [M]
  - `rel_v` = 0.4 (–) — Sabin virus transmissibility relative to WPV [I]
  - `rel_vd` = 0.9 (–) — cVDPV transmissibility relative to WPV [I]
  - `p_rev` = 0.01 (–) — fraction of Sabin onward transmissions producing revertant (cVDPV) virus [I]
  - `p_opv` = 0.8 (–) — routine OPV coverage at birth (cohort) [I]
  - `take_v` = 0.5 (–) — per-dose take, intestinal infection [I]
  - `omega_p` = 1/(5*365) (day^-1) — waning of intestinal immunity [I]
  - `par_ratio` = 0.005 (–) — paralysis-to-infection ratio, about 1/200 for WPV1 [M, widely cited]
- Derived variables:
  - `bw` = R0_w*gam_p
  - `bv` = rel_v*R0_w*gam_p
  - `bvd` = rel_vd*R0_w*gam_p
  - `inc_par` = par_ratio*(bw*IW + bvd*IVD)/N*S — paralytic incidence
- Flows:
  - ` > S : mu_h*N*(1 - p_opv*take_v)` (absolute)
  - ` > IV : mu_h*N*p_opv*take_v` (absolute; vaccinated newborns infected with Sabin virus)
  - `S > IW : bw*IW/N`
  - `S > IV : bv*IV/N*(1 - p_rev)`
  - `S > IVD : bv*IV/N*p_rev + bvd*IVD/N`
  - `IW > R : gam_p`
  - `IV > R : gam_p`
  - `IVD > R : gam_p`
  - `R > S : omega_p`
  - Deaths: `S > : mu_h`, `IW > : mu_h`, `IV > : mu_h`, `IVD > : mu_h`, `R > : mu_h`
- Notes:
  - Cross-immunity among WPV, Sabin and cVDPV is assumed complete: one R class.
  - Reversion is treated as a per-transmission probability, which is a gross simplification of stepwise genetic reversion.
  - To model OPV cessation, set `p_opv` to 0 after a time using `ifelse(t > t_stop, 0, 0.8)`.
  - Treat the whole parameter set as illustrative.

### rotavirus_pitzer2009 — Rotavirus with maternal immunity and repeated infections of decreasing severity
- Disease/area: Rotavirus gastroenteritis | Category: Water, enteric & NTDs
- Source: Pitzer VE, Viboud C, Simonsen L, Steiner C, Panozzo CA, Alonso WJ, Miller MA, Glass RI, Glasser JW, Parashar UD, Grenfell BT. Demographic variability, vaccination, and the spatiotemporal dynamics of rotavirus epidemics. *Science* 2009;325(5938):290–294. doi:10.1126/science.1172330 — [Princeton record](https://collaborate.princeton.edu/en/publications/demographic-variability-vaccination-and-the-spatiotemporal-dynami/). Relative-susceptibility values come from Velázquez et al. (*NEJM* 1996), as reproduced in a later model's parameter table — [PMC9746532](https://pmc.ncbi.nlm.nih.gov/articles/PMC9746532)
- Time unit: weeks; simulate 52×30 weeks.
- Compartments (0, 1, 2, 3+ prior infections):
  - `M` — maternal immunity — 0.01*N0
  - `S0` — naive susceptible — 0.05*N0
  - `I1` — first infection — 10
  - `R1` — temporary immunity after the 1st infection — 0
  - `S1` — susceptible after the 1st infection — 0.1*N0
  - `I2` — second infection — 0
  - `R2` — temporary immunity after the 2nd infection — 0
  - `S2` — susceptible after the 2nd infection — 0.2*N0
  - `I3` — third or later infection — 0
  - `R3` — temporary immunity after the 3rd or later infection — 0
  - `S3` — susceptible after 3+ infections — 0.64*N0
- Parameters:
  - `mu_h` = 0.0003 (week^-1) — birth/death. The paper used state-specific US birth rates (about 14–15 per 1000 per year) [I]
  - `omega_m` = 1/13 (week^-1) — loss of maternal immunity (13 weeks) [V via PMC9746532; M that Pitzer 2009 used it]
  - `gam_1` = 1 (week^-1) — recovery from the first infection (1 week) [M]
  - `gam_2` = 2 (week^-1) — recovery from later infections (half as long) [M]
  - `rho_inf` = 0.5 (–) — relative infectiousness of 2nd+ infections [M]
  - `sus_1` = 0.62 (–) — relative susceptibility after 1 infection [V, Velázquez via PMC9746532]
  - `sus_2` = 0.37 (–) — relative susceptibility after 2+ infections [V, Velázquez via PMC9746532]. Pitzer may have used 0.35 [M].
  - `omega_r` = 1/52 (week^-1) — waning of temporary immunity (1 year) [M]
  - `R0_rv` = 15 (–) — R0; fitted values were roughly 10–25 [M]
  - `amp_s` = 0.15 (–) — seasonal amplitude [M/I]
  - `phi_s` = 0 (week) — seasonal phase [I]
  - `ps1`, `ps2`, `ps3` = 0.47, 0.25, 0.32 (–) — proportion symptomatic for the 1st, 2nd and 3rd+ infections. Velázquez 1996-type values; the 4th infection is about 0.2 [M]. Used for outputs only.
- Derived variables:
  - `b0` = R0_rv*(gam_1 + mu_h)
  - `bt` = b0*(1 + amp_s*cos(6.2831853*(t - phi_s)/52))
  - `lam` = bt*(I1 + rho_inf*(I2 + I3))/N
  - `inc_sev` = lam*(ps1*S0 + ps2*sus_1*S1 + ps3*sus_2*S2) — incidence of rotavirus gastroenteritis (RVGE)
- Flows:
  - ` > M : mu_h*N` (absolute)
  - `M > S0 : omega_m`
  - `S0 > I1 : lam`
  - `I1 > R1 : gam_1`
  - `R1 > S1 : omega_r`
  - `S1 > I2 : sus_1*lam`
  - `I2 > R2 : gam_2`
  - `R2 > S2 : omega_r`
  - `S2 > I3 : sus_2*lam`
  - `I3 > R3 : gam_2`
  - `R3 > S3 : omega_r`
  - `S3 > I3 : sus_2*lam`
  - Deaths at `mu_h` from every human compartment: `M > : mu_h`, `S0 > : mu_h`, and so on (11 flows)
- Notes:
  - The published model is age-structured and has infection stages up to the 4th. This spec merges stages 3+, giving 11 compartments.
  - Vaccine extension: route a share of births or `M` outflow to `S1`, as a vaccine acting like a first infection [M].
  - Check the durations against the Science supplementary material (SOM).

### norovirus_simmons2013 — Norovirus SEIAR with asymptomatic reinfection and waning disease immunity
- Disease/area: Norovirus gastroenteritis | Category: Water, enteric & NTDs
- Source: Simmons K, Gambhir M, Leon J, Lopman B. Duration of immunity to norovirus gastroenteritis. *Emerg Infect Dis* 2013;19(8):1260–1267. doi:10.3201/eid1908.130472 — [CDC EID](https://wwwnc.cdc.gov/eid/article/19/8/13-0472_article). The parameter Table 3 and the technical appendix were not readable this session.
- Time unit: days; simulate 30 years (10,950 days) to equilibrium.
- Compartments:
  - `S` — susceptible to infection and disease — 0.3*N0
  - `E` — exposed (presymptomatic) — 0
  - `I` — symptomatic — 100
  - `A` — asymptomatic infected (post-symptomatic shedding, or reinfection of immune people) — 0
  - `R` — immune to disease but not to infection — 0.7*N0
- Parameters:
  - `inc_d` = 1 (day) — incubation, E → I [V, "about 1 day"]
  - `sym_d` = 2 (day) — symptomatic duration, I → A [V, "1 day, half the symptomatic phase"]
  - `shed_d` = 10 (day) — asymptomatic shedding, A → R [V, Model B]
  - `eps_ea` = 0.05 (–) — relative infectiousness of E and A: 0.05 in Model B, 0.25 in Model C, 0 in Model A [V]
  - `theta_d` = 5.1*365 (day) — duration of disease immunity, R → S. Fitted values: 5.1 years (Model A/B), 8.7 (C), 4.1 (D, 20% genetically resistant) [V]
  - `R0_nv` = 1.7 (–) — R0. Reported values: 1.64–1.88 for most models, 7.16 for C [V]
  - `mu_h` = 1/(70*365) (day^-1) — birth = death [I]
- Derived variables:
  - `Dinf` = eps_ea*inc_d + sym_d + eps_ea*shed_d — effective infectious duration (2.55 days)
  - `b_nv` = R0_nv/Dinf — about 0.667/day
  - `lam` = b_nv*(I + eps_ea*(E + A))/N
- Flows:
  - ` > S : mu_h*N` (absolute)
  - `S > E : lam`
  - `E > I : 1/inc_d`
  - `I > A : 1/sym_d`
  - `A > R : 1/shed_d`
  - `R > A : lam` (asymptomatic reinfection with boosting)
  - `R > S : 1/theta_d`
  - Deaths at `mu_h` from S, E, I, A and R
- Notes:
  - The paper's force of infection used age-specific contact rates: q1 for children under 5 and q2 for others. This is a single-age simplification, so the effective R0 is approximate.
  - All first infections are symptomatic [V].
  - Model D adds `G` (genetically resistant, 20% of births, never infected) [V].
  - Seasonal forcing of 6% was tested but excluded from the final estimates [V].

### hepatitisA_seir_vacc — Hepatitis A with age-independent symptomatic fraction and vaccination (simplified Van Effelterre-type)
- Disease/area: Hepatitis A (faecal–oral) | Category: Water, enteric & NTDs
- Source: Van Effelterre TP, Zink TK, Hoet BJ, Hausdorff WP, Rosenthal P. A mathematical model of hepatitis A transmission in the United States indicates value of universal childhood immunization. *Clin Infect Dis* 2006;43(2):158–164. doi:10.1086/505115 [M; DOI seen in search results but not opened] — [search lead](https://zendy.io/pdf-viewer/10.1086/505115). The original is age-structured and was calibrated to seroprevalence data. This spec is a simplified, illustrative structure.
- Time unit: days; simulate 50 years.
- Compartments:
  - `S` — susceptible — 0.5*N0
  - `E` — incubating — 0
  - `IA` — asymptomatic/anicteric infection — 10
  - `IS` — symptomatic (icteric) — 10
  - `R` — immune (lifelong) — 0.5*N0
  - `V` — vaccinated — 0
- Parameters (all [I]/[M]; clinical durations from standard HAV natural history):
  - `inc_h` = 28 (day) — incubation, range 15–50 [M]
  - `dur_h` = 21 (day) — infectious period [M/I]
  - `p_sym_h` = 0.3 (–) — fraction symptomatic. About 30% in young children rising to above 70% in adults [M]; the single value is [I].
  - `rel_ia` = 1 (–) — relative infectiousness of asymptomatic infections. Children are key transmitters [I]
  - `R0_ha` = 2.5 (–) — R0 [I]
  - `mu_h` = 1/(75*365) (day^-1) — birth/death [I]
  - `p_vacc` = 0 (–) — fraction of newborns vaccinated; try 0.8 [I]
  - `eff_hv` = 0.95 (–) — two-dose vaccine efficacy [M]
- Derived variables:
  - `b_ha` = R0_ha/dur_h
  - `lam` = b_ha*(IS + rel_ia*IA)/N
- Flows:
  - ` > S : mu_h*N*(1 - p_vacc*eff_hv)` (absolute)
  - ` > V : mu_h*N*p_vacc*eff_hv` (absolute)
  - `S > E : lam`
  - `E > IA : (1 - p_sym_h)/inc_h`
  - `E > IS : p_sym_h/inc_h`
  - `IA > R : 1/dur_h`
  - `IS > R : 1/dur_h`
  - Deaths at `mu_h` from all six compartments
- Notes:
  - Optional environmental route: add `W` (not in N) with ` > W : sh*(IS + IA)` (absolute), `W > : del` and an extra term `bw*W` in `lam`, mirroring the existing HEV model.
  - The age-dependent symptomatic fraction is the main lost feature, and it matters for the paradox of rising disease with falling transmission.

### schisto_snail_mwb — Schistosomiasis: human mean worm burden with SEI snails (Macdonald / Anderson–May hybrid)
- Disease/area: Schistosomiasis (S. haematobium / S. mansoni) | Category: Water, enteric & NTDs
- Source: Structure follows the Macdonald (1965) / Anderson & May (1985, 1991) mean-worm-burden + snail-infection framework, as implemented in, for example, Sokolow SH et al. *PNAS* 2015;112(31):9650–9655, doi:10.1073/pnas.1502651112 [M]. Related human–snail ODE models: Feng Z, Li C-C, Milner FA. *Math Biosci* 2002;177–178:271–286 (snail density dependence and infection age), and Feng Z, Milner FA, Minchella D. *Appl Math Lett* 2004 (parameter estimation, e.g. A_h ≈ 8/yr, A_s ≈ 25/yr) — [Purdue PDF list](https://www.math.purdue.edu/~fengz/documents/publications/aml04.pdf). Note: the Feng "2004" schistosomiasis model named in the brief is the AML parameter-estimation paper; the model itself is in the 2002 Math Biosci paper [V via search snippet]. Snail and worm values below are [I]/[M].
- Time unit: days; simulate 30 years.
- Compartments:
  - `Wm` — mean adult worm burden per human (worms/person) — 50 — not in N
  - `Sn` — susceptible snails — 7000 — not in N
  - `En` — exposed (prepatent) snails — 300 — not in N
  - `In` — infected (cercariae-shedding) snails — 400 — not in N
- Parameters:
  - `H` = 1000 (persons) — human population, held constant [I]
  - `r_s` = 0.1 (day^-1) — snail intrinsic birth rate [M/I]
  - `K_s` = 10000 (snails) — snail carrying capacity [I]
  - `mu_n` = 1/60 (day^-1) — snail death rate [M]
  - `sig_n` = 1/40 (day^-1) — snail prepatent period (about 4–6 weeks) [M]
  - `mu_i` = 1/10 (day^-1) — infected snail death rate [M]
  - `beta_sn` = 7.5e-7 (day^-1 per mated pair) — miracidial infection of snails. [I] chosen to give about 5% snail prevalence at Wm = 50.
  - `lam_h` = 1e-4 (worms day^-1 per infected snail per person) — cercarial establishment [I]
  - `mu_wm` = 1/(3.3*365) (day^-1) — adult worm death (about 3–4 year life span) [M]
  - `mu_h` = 1/(60*365) (day^-1) — human death rate [I]
  - `k_agg` = 0.2 (–) — negative binomial aggregation [M/I]
  - `tau_t` = 0 (day^-1) — continuous praziquantel treatment (MDA) rate. Pulsed MDA alternative: `ifelse((t - 365*floor(t/365)) < 7, -log(1 - 0.75*0.86)/7, 0)` [I; 86% worm kill is an M value for praziquantel]
- Derived variables:
  - `phi_m` = Wm/(Wm + 2*k_agg + 1) — crude mating-probability approximation. The exact Macdonald mating function is an integral; this is [I].
  - `Nsn` = Sn + En + In
  - `miracidia` = 0.5*Wm*H*phi_m — mated pairs, a proxy for egg output
- Flows:
  - ` > Sn : r_s*(1 - Nsn/K_s)*(Sn + En)` (absolute; infected snails castrated)
  - `Sn > En : beta_sn*miracidia`
  - `En > In : sig_n`
  - `Sn > : mu_n`
  - `En > : mu_n`
  - `In > : mu_i`
  - ` > Wm : lam_h*In` (absolute)
  - `Wm > : mu_wm + mu_h + tau_t`
- Notes:
  - Wm is a per-capita mean, not a population. The Wm flows are bookkeeping only. N has no human compartments, so use the parameter H.
  - The breakpoint (Allee effect from mating) is reproduced only qualitatively.
  - Prawn predation (Sokolow) is an add-on term `Sn > : pred*P` [M].
  - Calibrate `lam_h` and `beta_sn` to the local prevalence or egg counts.

### sth_mean_worm — Soil-transmitted helminths: mean worm burden with environmental reservoir of infective stages (Anderson & May)
- Disease/area: Ascaris / Trichuris / hookworm | Category: Water, enteric & NTDs
- Source: Anderson RM, May RM. *Infectious Diseases of Humans: Dynamics and Control.* Oxford University Press, 1991, chapters 16–17 (mean-worm-burden ODE with free-living reservoir, density-dependent fecundity and negative-binomial aggregation) [M]. Implemented for MDA in Truscott JE, Hollingsworth TD, Brooker SJ, Anderson RM. *Parasites & Vectors* 2014;7:266, doi:10.1186/1756-3305-7-266 [M]. Values are [M]/[I].
- Time unit: years; simulate 30 years.
- Compartments:
  - `Mw` — mean worm burden per host — 10 — not in N
  - `Lr` — infective stages in the environment, normalised per host — 5 — not in N
- Parameters:
  - `R0_sth` = 3 (–) — basic reproduction number (Ascaris about 1–5) [M]
  - `mu_w` = 1 (yr^-1) — adult worm death (Ascaris about 1 year) [M]
  - `mu_h` = 0.02 (yr^-1) — host death [I]
  - `mu_l` = 6 (yr^-1) — infective-stage loss (about 2 months) [I]
  - `k_agg` = 0.5 (–) — aggregation [M]
  - `z_dd` = 0.93 (–) — per-worm density-dependent fecundity factor, exp(−γ) [M/I]
  - `lam_e` = 1 (stages per worm per yr, normalised) — egg output [I; normalisation]
  - `tau_m` = 0 (yr^-1) — continuous MDA kill rate. Pulsed alternative: `ifelse((t - floor(t)) < 0.02, -log(1 - 0.75*0.95)/0.02, 0)` (75% coverage, 95% efficacy of albendazole against Ascaris) [I/M]
- Derived variables:
  - `beta_l` = R0_sth*(mu_w + mu_h)*mu_l/lam_e — uptake rate of infective stages (18.36/yr)
  - `fdd` = (1 + Mw*(1 - z_dd)/k_agg)^(-(k_agg + 1)) — density-dependent fecundity under negative-binomial aggregation
  - `prev_sth` = 1 - (1 + Mw/k_agg)^(-k_agg) — prevalence of infection
- Flows:
  - ` > Mw : beta_l*Lr` (absolute)
  - `Mw > : mu_w + mu_h + tau_m`
  - ` > Lr : lam_e*Mw*fdd` (absolute)
  - `Lr > : mu_l`
- Notes:
  - Uptake by hosts does not deplete Lr; Anderson & May neglect this.
  - Mating is ignored, which suits hookworm and Ascaris at moderate burdens.
  - Hookworm: mu_w ≈ 0.5/yr and mu_l ≈ 26/yr [M]. Trichuris: mu_w ≈ 1/yr [M].

### trachoma_sis_mda — Trachoma SIS with periodic mass azithromycin (Lietman 1999)
- Disease/area: Ocular Chlamydia trachomatis | Category: Water, enteric & NTDs
- Source: Lietman T, Porco T, Dawson C, Blower S. Global elimination of trachoma: how frequently should we administer mass chemotherapy? *Nat Med* 1999;5(5):572–576. doi:10.1038/8451 [M DOI]. Its conclusions [V]: treat annually where prevalence in children is below about 35%, and biannually in hyperendemic areas (>50% in children). The original was a stochastic model; the 2011 follow-up fitted non-linear stochastic SIS models (Lietman et al., *Epidemics* 2011) — [search summary](https://www.curehunter.com/public/pubmed10229236.do). Rates are [I].
- Time unit: weeks; simulate 10 years (520 weeks).
- Compartments:
  - `S` — uninfected — 500
  - `I` — infected with ocular chlamydia — 500
- Parameters:
  - `dur_tr` = 26 (week) — mean duration of infection in children [I]
  - `R0_tr` = 2 (–) — gives about 50% endemic prevalence (1 − 1/R0) [I]
  - `cover_m` = 0.8 (–) — MDA coverage [I]
  - `eff_m` = 0.95 (–) — azithromycin efficacy [M]
  - `int_m` = 52 (week) — treatment interval: 52 annual, 26 biannual [V that these intervals were compared]
  - `t_start` = 52 (week) — first MDA round [I]
  - `win_m` = 1 (week) — length of the treatment pulse [I]
- Derived variables:
  - `b_tr` = R0_tr/dur_tr
  - `mda_on` = ifelse(t >= t_start & (t - t_start - int_m*floor((t - t_start)/int_m)) < win_m, 1, 0)
  - `mda_rate` = mda_on*(-log(1 - cover_m*eff_m)/win_m) — a pulse giving the cover×eff kill in each round
- Flows:
  - `S > I : b_tr*I/N`
  - `I > S : 1/dur_tr + mda_rate`
- Notes:
  - Deterministic SIS cannot reach exact zero, so measure "elimination" as prevalence below 1/N.
  - Use small `win_m` together with the `hmax` solver option so pulses are not stepped over.
  - Imported infection: add ` > I : imp` (absolute) and `I > : imp/N`, or a background term in the force of infection.
- Yaws variant: the same SIS + MDA structure applies with azithromycin total community treatment.
  - Example values: dur ≈ 26 weeks, latent class optional [I].
  - Published yaws transmission models: e.g. Marks M et al. / Fitzpatrick C et al. 2014 *PLoS NTD*, doi:10.1371/journal.pntd.0003237 [M, not verified].

### crypto_oocyst — Cryptosporidiosis with environmental oocysts and exponential dose–response
- Disease/area: Cryptosporidium (also adaptable to Giardia cysts) | Category: Water, enteric & NTDs
- Source: Author-constructed SIRS + environment structure [I]. The dose–response is Haas CN et al., exponential model with r ≈ 0.004 per oocyst (IOWA strain): *J Am Water Works Assoc* 1996;88(9):131–136 [M]. Giardia exponential r ≈ 0.0199 (Rose, Haas & Regli 1991, *Am J Public Health* 81:709–713) [M]. All transmission scaling is [I].
- Time unit: days; simulate 3 years.
- Compartments:
  - `S` — susceptible — 9990
  - `I` — infected/shedding — 10
  - `R` — immune (partial, waning) — 0
  - `Oc` — oocysts per litre in the community water source — 0 — not in N
- Parameters:
  - `exp_c` = 1 (day^-1) — water-consumption events per day [I]
  - `vol_d` = 1 (L) — volume ingested per event [I]
  - `r_dose` = 0.004 (per oocyst) — exponential dose–response parameter [M]
  - `gam_c` = 1/14 (day^-1) — recovery and end of shedding [M/I]
  - `sh_vol` = 2e-4 (oocysts L^-1 day^-1 per infected) — shedding divided by water volume. [I] chosen to give R0 ≈ 2.
  - `mu_o` = 0.05 (day^-1) — oocyst inactivation and flushing [I]
  - `omega_c` = 1/365 (day^-1) — waning immunity [I]
  - `mu_h` = 1/(70*365) (day^-1) — birth/death [I]
  - `rem_t` = 0 (–) — fractional log-removal by treatment, for example 0.99 [I]
- Derived variables:
  - `dose` = vol_d*Oc*(1 - rem_t)
  - `lam` = exp_c*(1 - exp(-r_dose*dose)) — daily infection probability as a hazard
  - `R0_cr` = exp_c*r_dose*vol_d*sh_vol*N/(mu_o*(gam_c + mu_h)) — linearised R0 (about 2)
- Flows:
  - ` > S : mu_h*N` (absolute)
  - `S > I : lam`
  - `I > R : gam_c`
  - `R > S : omega_c`
  - Deaths at `mu_h` from S, I and R
  - ` > Oc : sh_vol*I` (absolute)
  - `Oc > : mu_o`
- Notes:
  - Add a pulse ` > Oc : ifelse(t > 100 & t < 107, pulse_o, 0)` (absolute) to mimic a contamination event, as in the Milwaukee 1993 outbreak.
  - For Giardia, use `r_dose` = 0.0199 and longer cyst survival.

### lepto_rodent_env — Leptospirosis: rodent reservoir, environmental leptospires, humans as dead-end hosts
- Disease/area: Leptospirosis | Category: Water, enteric & NTDs
- Source: Structure after Holt J, Davis S, Leirs H. A model of leptospirosis infection in an African rodent to determine risk to humans: seasonal fluctuations and the impact of rodent control. *Acta Trop* 2006;99(2–3):218–225. doi:10.1016/j.actatropica.2006.08.003 [M]. Environmental compartment as in later environment-mediated models, e.g. Baca-Carrasco et al. 2015 [M]. Values are [I].
- Time unit: days; simulate 5 years.
- Compartments:
  - `Sh` — susceptible humans — 9990
  - `Ih` — infected (clinical) humans — 0
  - `Rh` — recovered humans — 0
  - `Sr` — susceptible rodents — 900 — not in N
  - `Ir` — infected (chronically shedding) rodents — 100 — not in N
  - `Bl` — leptospires in the environment (arbitrary units) — 0 — not in N
- Parameters:
  - `b_r` = 1/180 (day^-1) — rodent birth rate per rodent [I]
  - `K_r` = 1000 (rodents) — rodent carrying capacity [I]
  - `mu_r` = 1/365 (day^-1) — rodent death [I]
  - `bd_r` = 0.002 (day^-1 per infected rodent) — direct rodent transmission, scaled by `Ir/Nr` [I]
  - `be_r` = 0.05 (day^-1) — rodent infection from the environment [I]
  - `sh_r` = 1 (units day^-1 per infected rodent) — shedding [I]
  - `del_b` = 0.1 (day^-1) — leptospire decay [I]
  - `K_b` = 1000 (units) — half-saturation [I]
  - `be_h` = 1e-4 (day^-1) — human exposure rate. It scales with flooding: multiply by `(1 + amp_f*cos(6.2831853*t/365))` [I]
  - `gam_h` = 1/14 (day^-1) — human recovery [M/I]
  - `cfr_h` = 0.005 (day^-1) — disease death hazard (CFR about 5–10% in severe cases) [I]
  - `omega_h` = 1/(5*365) (day^-1) — waning serovar immunity [I]
  - `mu_h` = 1/(70*365) (day^-1) — human birth/death [I]
- Derived variables:
  - `Nr` = Sr + Ir
  - `envB` = Bl/(K_b + Bl)
- Flows:
  - ` > Sr : b_r*Nr*(1 - Nr/K_r)` (absolute; vertical transmission ignored) [I]
  - `Sr > Ir : bd_r*Ir/Nr + be_r*envB`
  - `Sr > : mu_r`
  - `Ir > : mu_r`
  - ` > Bl : sh_r*Ir` (absolute)
  - `Bl > : del_b`
  - ` > Sh : mu_h*N` (absolute)
  - `Sh > Ih : be_h*envB`
  - `Ih > Rh : gam_h`
  - `Ih > : cfr_h`
  - `Rh > Sh : omega_h`
  - Deaths at `mu_h` from Sh, Ih and Rh
- Notes:
  - Humans do not feed back into transmission.
  - Rodent control: add `Sr > : ctrl` and `Ir > : ctrl`.
  - The logistic rodent birth term can go negative if Nr > K_r. Wrap it in `max(0, ...)`.

### scabies_mda — Scabies SIS with egg-stage reinfection and two-dose ivermectin MDA
- Disease/area: Scabies (Sarcoptes scabiei) | Category: Water, enteric & NTDs (skin NTD)
- Source: Lydeamore MJ, Campbell PT, Regan DG, Tong SYC, Andrews RM, Steer AC, Romani L, Kaldor JM, McVernon J, McCaw JM. A biological model of scabies infection dynamics and treatment informs mass drug administration strategies to increase the likelihood of elimination. *Math Biosci* 2019;309:163–173. doi:10.1016/j.mbs.2018.08.007 [M]. Simplified here to two infected states; values are [I]/[M].
- Time unit: days; simulate 3 years.
- Compartments:
  - `S` — uninfested — 0.7*N0
  - `I` — infested (active mites) — 0.3*N0
  - `Ie` — post-treatment carriage of surviving eggs, which re-hatch unless a second dose is given — 0
- Parameters:
  - `R0_sc` = 1.5 (–) — R0 [I]
  - `dur_sc` = 365 (day) — untreated infestation duration, chronic [I]
  - `hatch` = 1/7 (day^-1) — egg hatch/maturation back to active infestation (ivermectin is not ovicidal) [M]
  - `cover_s` = 0.85 (–) — MDA coverage [I]
  - `eff_s` = 0.95 (–) — ivermectin kill of active mites [I]
  - `t1` = 100 (day) — first MDA dose [I]
  - `t2` = 110 (day) — second dose, 7–14 days later [M on interval]
  - `win_s` = 1 (day) — pulse width [I]
- Derived variables:
  - `b_sc` = R0_sc/dur_sc
  - `dose_on` = ifelse((t >= t1 & t < t1 + win_s) | (t >= t2 & t < t2 + win_s), 1, 0)
  - `mda_k` = dose_on*(-log(1 - cover_s*eff_s)/win_s)
- Flows:
  - `S > I : b_sc*I/N`
  - `I > S : 1/dur_sc`
  - `I > Ie : mda_k` (mites killed, eggs survive)
  - `Ie > I : hatch` (re-infestation from surviving eggs)
  - `Ie > S : mda_k` (second dose clears the hatched juveniles)
- Notes:
  - Ie is assumed non-infectious.
  - Permethrin, being ovicidal, can be modelled as `I > S : mda_k` instead.
  - The `|` operator inside `ifelse` is allowed if the generator accepts comparisons. Otherwise use `min(1, (t >= t1)*(t < t1 + win_s) + (t >= t2)*(t < t2 + win_s))`.

### Inferences
- Codeço 2001 is the only model whose parameter set could be fully checked. Its published S_C values (6600 and 13,200) follow from S_C = K·r·(mb−nb)/(e·a), which confirms that the dose–response is B/(K+B).
- Most published models in this area are age-structured: Pitzer typhoid and rotavirus, Simmons norovirus, Van Effelterre HAV, and Truscott STH. The single-age specs above therefore reproduce qualitative dynamics, not the published calibrations.
- Pulsed MDA (trachoma, scabies, STH, schistosomiasis) fits the ODE-only format with `ifelse`/`floor` windows and a log-kill rate −log(1 − coverage·efficacy)/window. The solver step must be limited (`hmax` ≤ window).

### Gaps
- Lymphatic filariasis and onchocerciasis were skipped. Published models (EPIFIL, LYMFASIM, ONCHOSIM, EPIONCHO) are vector-based, age-structured and/or individual-based, and do not reduce cleanly to ≤14 compartments without vectors.
- Guinea worm was skipped. No standard published compartmental ODE model was found or verified. An [I] host–copepod SEI-L structure is possible: about 1-year prepatency, ~3-week emergence, larvae in water.
- The numeric tables and supplements could not be retrieved for Andrews & Basu 2011 (PMC page blocked), Pitzer 2009 (Science SOM), Pitzer 2014 (Table 1 image / Text S1), Simmons 2013 (Table 3 / technical appendix), Feng et al. (PDF not text-extractable), and Lietman 1999 (paywalled).

---

## What is the primary source (authors, year, journal, DOI/URL), and where do the parameter values come from?

### Takeaway
The primary citations are listed per model above. Parameter values are traced to the primary paper only where marked [V]. Elsewhere they come from recalled secondary literature [M] or are illustrative [I]. The supplementary tables named below should be pulled before the library labels values as "published".

### Cited Findings
- Codeço 2001: *BMC Infect Dis* 1:1, doi:10.1186/1471-2334-1-1, PMC29087. Table 1 defines H, n, a, K, r, nb, mb and e. Table 2 gives the three hypothetical communities (n = 0.0001/0.0001/0.001; a = 0.5/1/1; K = 1e6; r = 0.2; nb−mb = −0.33; e = 10; S_C = 13,200/6600/6600) — [Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC29087/fullTextXML)
- Andrews & Basu 2011: *Lancet* 377:1248–1255. They calibrated the water consumption rate, symptomatic mortality and initial environmental concentration by MCMC (100,000 iterations), and projected 779,000 cases and 11,100 deaths — [PMC3172163](https://pmc.ncbi.nlm.nih.gov/articles/PMC3172163); [CIDRAP](https://www.cidrap.umn.edu/foodborne-disease/researchers-predict-high-2011-cholera-toll-haiti)
- Pitzer et al. 2014 typhoid: R0 = 2.8 (2.6–2.9), with 90% short-cycle; r = 0.01; ε ≈ 0; f = 0.0052; 10% of infections clinical; carriers 0.3–10.1% by age, lifelong; vaccine efficacies 48%, 80% and 95.6%; sanitation scenario over years 5–35 — [PLoS NTD](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0002642)
- Pitzer et al. 2009 rotavirus: *Science* 325:290–294, doi:10.1126/science.1172330. Birth-rate variation explains the spatiotemporal epidemic pattern; reduced and lagged post-vaccination epidemics are predicted — [Princeton](https://collaborate.princeton.edu/en/publications/demographic-variability-vaccination-and-the-spatiotemporal-dynami/)
- Rotavirus secondary parameters: maternal immunity 13 weeks (Linhares); relative risk of infection 0.62 after the 1st and 0.37 after the 2nd+ infection (Velázquez); latent period 1 day — [PMC9746532](https://pmc.ncbi.nlm.nih.gov/articles/PMC9746532)
- Simmons et al. 2013 norovirus: *EID* 19(8). SEIAR structure, with R → A reinfection and R → S waning at 1/θ. Duration of immunity 4.1–8.7 years. E and A relative infectiousness 5% (B) or 25% (C); shedding 10 days; R0 1.64–1.88 for most models — [CDC EID](https://wwwnc.cdc.gov/eid/article/19/8/13-0472_article)
- Discrepancy: the abstract gives the Model C immunity CI as 6.8–11.3 years, while the Results section text gives 4.0–11.3 — [CDC EID](https://wwwnc.cdc.gov/eid/article/19/8/13-0472_article)
- Schistosomiasis, Feng group: the model is in Feng, Li & Milner *Math Biosci* 2002 (density dependence and age of infection in snails). Parameter estimation is in Feng, Milner & Minchella, *Appl Math Lett* 2004 (A_h ≈ 8/yr, A_s ≈ 25/yr per snippet) — [Purdue publications](https://www.math.purdue.edu/~fengz/documents/publications/mathbio02.pdf)
- Lietman et al. 1999: *Nat Med* 5(5):572–576. Recommends annual MDA where prevalence in children is below 35%, and biannual where it is above 50% — [abstract](https://www.curehunter.com/public/pubmed10229236.do); [Burnet summary of later work](https://burnet.edu.au/knowledge-and-media/research-articles/the-epidemiological-dynamics-of-infectious-trachoma-may-facilitate-elimination)
- Van Effelterre hepatitis A: the 2006 paper is in *Clin Infect Dis* (US universal childhood immunization), not *Epidemiology and Infection*. Later Van Effelterre models cover Mexico and Brazil — [PMC3551882](https://pmc.ncbi.nlm.nih.gov/articles/PMC3551882); [PMC5557237](https://pmc.ncbi.nlm.nih.gov/articles/PMC5557237)

### Inferences
- For library labelling: mark Codeço as "published values". Mark Pitzer 2014 and Simmons 2013 as "partly published, structure simplified". Mark all others as "structure after source; values illustrative".

### Gaps
- Not verified this session: DOIs marked [M] (Lietman 1999 doi:10.1038/8451, Fine & Carneiro 1999, Van Effelterre CID 2006 doi:10.1086/505115, Holt 2006, Lydeamore 2019, Sokolow 2015, Truscott 2014, Haas 1996), and every [M] parameter value (rotavirus durations and symptomatic fractions, typhoid durations and water decay, polio R0 and paralysis ratio, STH life spans, Cryptosporidium r = 0.004).
- Sources still to pull: the Science SOM for Pitzer 2009, Text S1 for Pitzer 2014, the EID technical appendix for Simmons 2013, and the Lancet appendix for Andrews & Basu 2011. These are the most important to check.
