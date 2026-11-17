# Classic textbook compartmental models plus intervention, treatment, behaviour and AMR extensions

Scope: 17 deterministic ODE model specifications that do not duplicate the already-implemented list. Every one is written in the generator's flow format ("FROM > TO : rate", per-capita by default, "absolute" marked, blank side = outside).

**Verification status, read first.** I verified these citations and DOIs by web search in this session: Funk et al. 2009, Gumel et al. 2004, Bonhoeffer et al. 1997 (including its equations and its Figure 2 parameters, read from PMC), Wang & Ruan 2004, Zhonghua & Yaohong 2009, Capasso & Serio 1978, Gomes et al. 2004, Nowak & May 1994, Kribs-Zaleta & Velasco-Hernández 2000, Granich et al. 2009 and Spicknall et al. 2013. The other DOIs and ISBNs come from my prior knowledge and are labelled "(from memory, not re-verified)". The companion-code sites (modelinginfectiousdiseases.org and the Warwick mirror) and the full Gumel 2004 PDF could not be fetched (connection refused, or a captcha). Most numeric parameter values below are therefore **illustrative**, and each one is labelled as such. Values attributed to a source are flagged when they are recalled rather than read.

**Identifier caution.** `beta` and `gamma` are also base-R function names (`beta()`, `gamma()`). Inside a deSolve `with(as.list(c(state, parms)), ...)` block they shadow the functions harmlessly, and the existing library almost certainly uses them already. If the generator's validator rejects R function names strictly, rename them to `beta_` / `gamma_` or `b_inf` / `g_rec`. I avoid `F, T, c, t, pi, C, D` (except as compartment D where noted) and `I()`-style collisions beyond `I`, which the existing models already use. `C()` is a stats function, so no compartment is named `C`, and treated compartments are named `Tr`.

## What are the compartments, flows, rate expressions, forces of infection and typical parameter values for each model?

### Takeaway
All 17 models below fit the flow format with at most 9 compartments. Each spec lists its compartments, flows, force of infection and parameters. Structure is taken from the cited source wherever possible. Parameter values are mostly illustrative but plausible, and each spec says which values come from the source: Bonhoeffer 1997 is the only model whose values were read directly from the paper, and the Gumel 2004 values are recalled and need checking against its Table 1.

### Cited Findings

---

### sitr_treatment — SEITR influenza-type model with antiviral treatment
- Disease/area: pandemic/seasonal influenza, antiviral treatment | Category: treatment
- Source: Brauer F. (2008) "Compartmental models in epidemiology", in Brauer, van den Driessche & Wu (eds), *Mathematical Epidemiology*, Lecture Notes in Mathematics 1945, Springer, pp. 19–79, doi:10.1007/978-3-540-78911-6_2 (from memory, not re-verified); also Brauer & Castillo-Chavez (2012) *Mathematical Models in Population Biology and Epidemiology*, 2nd ed., Springer, doi:10.1007/978-1-4614-1686-9 (from memory). The structure follows Brauer's "treatment model" (treated people are less infectious and recover faster). The values are illustrative.
- Time unit: days ; simulate 200 days
- Compartments:
  - `S` — susceptible — 9990
  - `E` — latent — 0
  - `I` — infectious, untreated — 10
  - `Tr` — infectious, under treatment — 0
  - `R` — recovered — 0
- Parameters:
  - `beta` = 0.5 (1/day) — transmission rate — illustrative
  - `kappa` = 1/1.5 (1/day) — progression from latent — illustrative (influenza latent period about 1–2 d)
  - `tau` = 0.3 (1/day) — treatment uptake rate — illustrative
  - `alpha` = 1/4 (1/day) — untreated recovery rate — illustrative
  - `eta` = 1/2.5 (1/day) — recovery rate under treatment — illustrative
  - `delta` = 0.5 (–) — relative infectiousness of treated people — illustrative
- Derived variables:
  - `lambda` = beta*(I + delta*Tr)/N — force of infection
  - `R0` = beta*(1/(alpha+tau) + delta*tau/((alpha+tau)*eta)) — basic reproduction number (about 1.26 with these values; set beta = 0.8 for about 2.0)
- Flows:
  - `S > E : lambda`
  - `E > I : kappa`
  - `I > Tr : tau`
  - `I > R : alpha`
  - `Tr > R : eta`
- Notes: Closed population, so N is constant. Use it to show how treatment coverage lowers both R0 and the final size. Setting tau = 0 recovers SEIR.

---

### sir_saturated_treatment — SIR with saturated (capacity-limited) treatment
- Disease/area: generic, limited treatment capacity | Category: treatment / backward bifurcation
- Source: Wang W. & Ruan S. (2004) "Bifurcations in an epidemic model with constant removal rate of the infectives", *J. Math. Anal. Appl.* 291:775–793, doi:10.1016/j.jmaa.2003.11.043 ([MaRDI](https://portal.mardi4nfdi.de/wiki/Publication:1433379)). That paper uses a constant, piecewise treatment rate. The saturated form h(I) = rI/(1+αI) comes from Zhonghua Z. & Yaohong S. (2009) "Qualitative analysis of a SIR epidemic model with saturated treatment rate", *J. Appl. Math. Comput.*, doi:10.1007/s12190-009-0315-9 ([Springer](https://link.springer.com/doi/10.1007/s12190-009-0315-9)). See also Zhang X. & Liu X. (2008) *J. Math. Anal. Appl.* 348:433–443, doi:10.1016/j.jmaa.2008.07.042 (from memory). The values are illustrative.
- Time unit: days ; simulate 3650 days (10 y)
- Compartments:
  - `S` — susceptible — 9990
  - `I` — infectious — 10
  - `R` — removed — 0
- Parameters:
  - `Lambda` = 0.391 (persons/day) — recruitment, equal to mu*10000 — illustrative
  - `mu` = 1/(70*365) (1/day) — natural death — illustrative
  - `beta` = 3e-5 (1/(person·day)) — mass-action transmission coefficient — illustrative
  - `gamma` = 0.1 (1/day) — natural recovery — illustrative
  - `d_dis` = 0.001 (1/day) — disease death — illustrative
  - `r_max` = 0.1 (1/day) — maximal per-capita treatment rate at low prevalence — illustrative
  - `alpha_t` = 0.01 (1/person) — saturation constant (capacity effect) — illustrative
- Derived variables:
  - `treat_pc` = r_max/(1 + alpha_t*I) — per-capita treatment rate, so the total treatment is r_max*I/(1+alpha_t*I)
- Flows:
  - ` > S : Lambda` (absolute)
  - `S > I : beta*I`
  - `I > R : gamma`
  - `I > R : treat_pc`
  - `S > : mu`
  - `I > : mu + d_dis`
  - `R > : mu`
- Notes: Mass-action (density-dependent) incidence as in the source papers. When treatment saturates, a backward bifurcation can occur, so an endemic state can persist with R0 < 1. Here R0 = beta*(Lambda/mu)/(mu+d_dis+gamma+r_max), about 1.5 with these values.

---

### sir_saturated_incidence — SIR with saturated incidence (Capasso–Serio)
- Disease/area: cholera (Bari 1973 epidemic), generic | Category: nonlinear incidence / behaviour
- Source: Capasso V. & Serio G. (1978) "A generalization of the Kermack–McKendrick deterministic epidemic model", *Math. Biosci.* 42:43–61, doi:10.1016/0025-5564(78)90006-8 ([MaRDI](https://portal.mardi4nfdi.de/wiki/Publication:1254217)). Hethcote 2000 (SIAM Rev. 42:599–653, doi:10.1137/S0036144500371907, from memory) reviews nonlinear incidences. The values are illustrative.
- Time unit: days ; simulate 365 days
- Compartments:
  - `S` — susceptible — 9990
  - `I` — infectious — 10
  - `R` — recovered — 0
- Parameters:
  - `beta` = 5e-5 (1/(person·day)) — transmission coefficient — illustrative
  - `alpha_s` = 0.002 (1/person) — saturation ("psychological"/crowding) constant — illustrative
  - `gamma` = 0.2 (1/day) — recovery — illustrative
- Derived variables:
  - `lambda` = beta*I/(1 + alpha_s*I) — saturating force of infection
- Flows:
  - `S > I : lambda`
  - `I > R : gamma`
- Notes: The incidence g(I)S with g bounded represents protective behaviour or saturation at high prevalence. R0 = beta*S0/gamma (about 2.5) is unchanged, but the peak is flattened. Add births and deaths (mu) if an endemic version is wanted.

---

### sir_awareness — SIR with aware (behaviour-changing) susceptibles
- Disease/area: generic emerging infection; awareness/behaviour | Category: behaviour
- Source: The concept and qualitative result come from Funk S., Gilad E., Watkins C., Jansen V.A.A. (2009) "The spread of awareness and its impact on epidemic outbreaks", *PNAS* 106(16):6872–6877, doi:10.1073/pnas.0810762106 ([LSHTM](https://researchonline.lshtm.ac.uk/id/eprint/611254)). This mean-field structure is close to Kiss I.Z., Cassell J., Recker M., Simon P.L. (2010) "The impact of information transmission on epidemic outbreaks", *Math. Biosci.* 225:1–10, doi:10.1016/j.mbs.2009.11.009 (from memory). **The structure is a simplified synthesis and the values are illustrative.** Funk et al.'s own model is network-based, with graded awareness levels.
- Time unit: days ; simulate 200 days
- Compartments:
  - `S` — unaware susceptible — 9990
  - `Sa` — aware susceptible (reduced contacts) — 0
  - `I` — infectious — 10
  - `R` — recovered — 0
- Parameters:
  - `beta` = 0.4 (1/day) — transmission rate — illustrative
  - `gamma` = 0.2 (1/day) — recovery — illustrative
  - `phi_a` = 2.0 (1/day) — awareness generation from contact with or news of cases — illustrative
  - `omega_f` = 0.05 (1/day) — awareness fading — illustrative
  - `sigma_a` = 0.3 (–) — relative susceptibility of aware susceptibles — illustrative
- Derived variables:
  - `lambda` = beta*I/N — force of infection
- Flows:
  - `S > Sa : phi_a*I/N`
  - `Sa > S : omega_f`
  - `S > I : lambda`
  - `Sa > I : sigma_a*lambda`
  - `I > R : gamma`
- Notes: Awareness is triggered only by prevalence, so the disease-free state has Sa = 0 and the epidemic threshold R0 = beta/gamma is unchanged. Awareness reduces the final size and the peak, which reproduces Funk et al.'s well-mixed finding that awareness shrinks outbreaks without changing the threshold. Optionally add `S > Sa : phi_s*Sa/N` for word-of-mouth spread among aware people.

---

### svirs_waning_vaccine — SVIRS with imperfect (leaky) and waning vaccine and waning natural immunity
- Disease/area: e.g. influenza, pertussis-like | Category: vaccination
- Source: The structure follows Arino J., McCluskey C.C., van den Driessche P. (2003) "Global results for an epidemic model with vaccination that exhibits backward bifurcation", *SIAM J. Appl. Math.* 64:260–276, doi:10.1137/S0036139902413829 (from memory, not re-verified), which has an SIRS model with vaccination of newborns and susceptibles, vaccine waning and a leaky vaccine. The values are illustrative.
- Time unit: years ; simulate 100 years
- Compartments:
  - `S` — susceptible — 9000
  - `V` — vaccinated (partially protected) — 0
  - `I` — infectious — 10
  - `R` — recovered — 990
- Parameters:
  - `mu` = 1/70 (1/yr) — birth = death rate — illustrative
  - `p_nb` = 0.5 (–) — fraction of newborns vaccinated — illustrative
  - `nu_v` = 0.1 (1/yr) — vaccination rate of susceptibles — illustrative
  - `theta_v` = 0.2 (1/yr) — vaccine waning (mean 5 y) — illustrative
  - `sigma_v` = 0.2 (–) — relative susceptibility of vaccinees (vaccine efficacy 80%) — illustrative
  - `beta` = 100 (1/yr) — transmission rate — illustrative
  - `gamma` = 365/7 (1/yr) — recovery (7-day infectious period) — illustrative
  - `omega_r` = 0.1 (1/yr) — natural immunity waning — illustrative
- Derived variables:
  - `lambda` = beta*I/N — force of infection
- Flows:
  - ` > S : (1 - p_nb)*mu*N` (absolute)
  - ` > V : p_nb*mu*N` (absolute)
  - `S > V : nu_v`
  - `V > S : theta_v`
  - `S > I : lambda`
  - `V > I : sigma_v*lambda`
  - `I > R : gamma`
  - `R > S : omega_r`
  - `S > : mu`
  - `V > : mu`
  - `I > : mu`
  - `R > : mu`
- Notes: N is constant. This differs from the existing "SVIR leaky vaccine" model in having vaccine waning, newborn vaccination and waning natural immunity. Backward bifurcation is possible when sigma_v > 0 and vaccination is high.

---

### siv_backward_bifurcation — SIS with vaccination (Kribs-Zaleta & Velasco-Hernández)
- Disease/area: generic SIS infection with imperfect waning vaccine | Category: vaccination / multistability
- Source: Kribs-Zaleta C.M. & Velasco-Hernández J.X. (2000) "A simple vaccination model with multiple endemic states", *Math. Biosci.* 164:183–201, doi:10.1016/S0025-5564(00)00003-1 ([MaRDI](https://portal.mardi4nfdi.de/wiki/Publication:1572822)). The paper's equations were not read in full, so the flow structure below is the standard SIV form described in its abstract and later literature. The values are illustrative.
- Time unit: years ; simulate 50 years
- Compartments:
  - `S` — susceptible — 9900
  - `I` — infectious — 100
  - `V` — vaccinated — 0
- Parameters:
  - `mu` = 0.02 (1/yr) — birth/death (turnover) — illustrative
  - `beta` = 20 (1/yr) — transmission — illustrative
  - `gamma` = 10 (1/yr) — recovery, returning to susceptible — illustrative
  - `phi_v` = 1.0 (1/yr) — vaccination rate — illustrative
  - `theta_v` = 0.1 (1/yr) — vaccine waning — illustrative
  - `sigma_v` = 0.3 (–) — relative susceptibility of the vaccinated — illustrative
- Derived variables:
  - `lambda` = beta*I/N
- Flows:
  - ` > S : mu*N` (absolute)
  - `S > I : lambda`
  - `I > S : gamma`
  - `S > V : phi_v`
  - `V > S : theta_v`
  - `V > I : sigma_v*lambda`
  - `S > : mu`
  - `I > : mu`
  - `V > : mu`
- Notes: This is the canonical small model where an endemic equilibrium coexists with a stable disease-free state when the vaccine-reduced R_phi < 1. Run it from different initial I to show bistability (the parameters may need tuning to land in the bistable region).

---

### sirv_booster — vaccination with waning and booster doses
- Disease/area: pertussis, COVID-19, HPV-like | Category: vaccination / booster
- Source: **Illustrative structure.** The rationale for waning and boosting follows Hethcote 2000 (SIAM Rev., doi:10.1137/S0036144500371907, from memory) and the vaccination chapters of Keeling & Rohani (2008) *Modeling Infectious Diseases in Humans and Animals*, Princeton University Press, ISBN 978-0-691-11617-4 (from memory). No single primary source supplies these values.
- Time unit: years ; simulate 50 years
- Compartments:
  - `S` — never vaccinated, never infected — 9000
  - `V1` — primary-vaccinated (high protection) — 0
  - `W` — waned (partially susceptible) — 0
  - `V2` — boosted — 0
  - `I` — infectious — 10
  - `R` — recovered — 990
- Parameters:
  - `mu` = 1/70 (1/yr) — birth/death — illustrative
  - `p_nb` = 0.9 (–) — infant primary coverage — illustrative
  - `omega1` = 1/4 (1/yr) — waning of primary protection — illustrative
  - `omega2` = 1/10 (1/yr) — waning after booster — illustrative
  - `omega_r` = 1/20 (1/yr) — waning of infection-derived immunity, flowing into W — illustrative
  - `b_rate` = 0.2 (1/yr) — booster uptake among waned — illustrative
  - `eps1` = 0.1 (–) — relative susceptibility in V1 — illustrative
  - `eps2` = 0.1 (–) — relative susceptibility in V2 — illustrative
  - `eps_w` = 0.5 (–) — relative susceptibility of the waned — illustrative
  - `beta` = 300 (1/yr) — transmission — illustrative
  - `gamma` = 365/21 (1/yr) — recovery — illustrative
- Derived variables:
  - `lambda` = beta*I/N
- Flows:
  - ` > S : (1 - p_nb)*mu*N` (absolute)
  - ` > V1 : p_nb*mu*N` (absolute)
  - `V1 > W : omega1`
  - `V2 > W : omega2`
  - `R > W : omega_r`
  - `W > V2 : b_rate`
  - `S > I : lambda`
  - `V1 > I : eps1*lambda`
  - `V2 > I : eps2*lambda`
  - `W > I : eps_w*lambda`
  - `I > R : gamma`
  - `S > : mu`
  - `V1 > : mu`
  - `W > : mu`
  - `V2 > : mu`
  - `I > : mu`
  - `R > : mu`
- Notes: N is constant. Compare b_rate = 0 against b_rate > 0. This differs from the existing SIRWS model (which uses natural boosting through exposure) because boosting here is programmatic.

---

### seiqhr_sars — SEIQHR quarantine and isolation (Gumel et al. 2004, SARS)
- Disease/area: SARS 2003 (Toronto, Hong Kong, Singapore, Beijing) | Category: quarantine / isolation
- Source: Gumel A.B. et al. (2004) "Modelling strategies for controlling SARS outbreaks", *Proc. R. Soc. B* 271:2223–2232, doi:10.1098/rspb.2004.2800, PMC1691853 ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC1691853/)). The search confirmed six classes: susceptible, asymptomatic (exposed), quarantined, symptomatic, isolated and recovered, with N their sum. **The full text could not be retrieved (captcha), so the parameter values below are recalled approximations of the Toronto set and must be checked against Gumel et al.'s Table 1 before use.**
- Time unit: days ; simulate 200 days
- Compartments:
  - `S` — susceptible — 4000000 (illustrative city size)
  - `E` — exposed/asymptomatic, mildly infectious — 10
  - `Q` — quarantined exposed — 0
  - `I` — symptomatic infectious — 1
  - `H` — isolated/hospitalised (Gumel's J) — 0
  - `R` — recovered — 0
  - `Dd` — cumulative SARS deaths — 0 — not in N
- Parameters:
  - `Pi_b` = 136 (persons/day) — recruitment (recalled; or set to mu*N0) — verify
  - `mu` = 1/(70*365) (1/day) — natural death — illustrative
  - `beta` = 0.2 (1/day) — transmission rate — recalled, verify
  - `eps_E` = 0.3 (–) — relative infectiousness of E — recalled, verify
  - `eps_Q` = 0 (–) — relative infectiousness of Q — recalled
  - `eps_H` = 0.36 (–) — relative infectiousness of the isolated (Toronto hospital transmission) — recalled, verify
  - `gamma1` = 0.1 (1/day) — quarantine rate of exposed — recalled, verify
  - `gamma2` = 0.5 (1/day) — isolation rate of symptomatic — recalled, verify
  - `kappa1` = 0.1 (1/day) — progression E→I — recalled
  - `kappa2` = 0.125 (1/day) — progression Q→H — recalled, verify
  - `sigma1` = 0.0337 (1/day) — recovery of I — recalled, verify
  - `sigma2` = 0.0386 (1/day) — recovery of H — recalled, verify
  - `d1` = 0.0068 (1/day) — SARS death rate in I — recalled, verify
  - `d2` = 0.0045 (1/day) — SARS death rate in H — recalled, verify
- Derived variables:
  - `lambda` = beta*(I + eps_E*E + eps_Q*Q + eps_H*H)/N
- Flows:
  - ` > S : Pi_b` (absolute)
  - `S > E : lambda`
  - `E > Q : gamma1`
  - `E > I : kappa1`
  - `Q > H : kappa2`
  - `I > H : gamma2`
  - `I > R : sigma1`
  - `H > R : sigma2`
  - `I > Dd : d1`
  - `H > Dd : d2`
  - `S > : mu`
  - `E > : mu`
  - `Q > : mu`
  - `I > : mu`
  - `H > : mu`
  - `R > : mu`
- Notes: The headline result is that optimal isolation of symptomatic cases is more effective than sub-optimal isolation combined with quarantine ([search summary of the abstract](https://pmc.ncbi.nlm.nih.gov/articles/PMC1691853/)). This model differs from the existing SIQR and SEIHRD models because it has both quarantine of exposed people and isolation of cases, and transmission continues from isolation. Vital dynamics are negligible over 200 days and can be dropped.

---

### hiv_test_and_treat — HIV universal test-and-treat (simplified Granich et al. 2009)
- Disease/area: HIV, generalised heterosexual epidemic (South Africa) | Category: test-and-treat
- Source: Granich R.M., Gilks C.F., Dye C., De Cock K.M., Williams B.G. (2009) "Universal voluntary HIV testing with immediate antiretroviral therapy as a strategy for elimination of HIV transmission: a mathematical model", *Lancet* 373:48–57, doi:10.1016/S0140-6736(08)61697-9 ([LSHTM](https://researchonline.lshtm.ac.uk/id/eprint/312973); [summary](https://natap.org/2010/HIV/032010_01.htm)). The original model has four CD4 stages and behaviour change. This version collapses them into one untreated stage. The ART relative infectiousness of 0.01 matches the roughly 99% reduction Granich assumed, and the 1.5%/yr dropout is recalled; both should be verified. Everything else is illustrative.
- Time unit: years ; simulate 50 years
- Compartments:
  - `S` — susceptible adults — 8000000
  - `Iu` — HIV+ undiagnosed — 1500000
  - `Id` — HIV+ diagnosed, not on ART — 300000
  - `A` — on ART — 200000
- Parameters:
  - `Lambda` = 250000 (persons/yr) — entry into the adult population — illustrative
  - `mu` = 1/35 (1/yr) — exit/natural death of adults — illustrative
  - `beta` = 0.5 (1/yr) — transmission rate from untreated infected — illustrative (chosen for about 15–20% prevalence)
  - `eps_A` = 0.01 (–) — relative infectiousness on ART — Granich (about 99% reduction), verify
  - `eps_d` = 0.8 (–) — relative infectiousness of the diagnosed untreated (behaviour change) — illustrative
  - `alpha_h` = 1/10 (1/yr) — HIV death rate untreated (about 10-year survival) — illustrative / textbook
  - `alpha_a` = 0.01 (1/yr) — HIV death rate on ART — illustrative
  - `test_rate` = 1.0 (1/yr) — testing rate (annual testing in the Granich scenario) — Granich scenario
  - `link_art` = 0.9 (–) — fraction starting ART immediately on diagnosis — illustrative
  - `art_init` = 0.5 (1/yr) — later ART initiation among the diagnosed — illustrative
  - `dropout` = 0.015 (1/yr) — ART dropout — recalled from Granich, verify
- Derived variables:
  - `lambda` = beta*(Iu + eps_d*Id + eps_A*A)/N
- Flows:
  - ` > S : Lambda` (absolute)
  - `S > Iu : lambda`
  - `Iu > A : test_rate*link_art`
  - `Iu > Id : test_rate*(1 - link_art)`
  - `Id > A : art_init`
  - `A > Id : dropout`
  - `S > : mu`
  - `Iu > : mu + alpha_h`
  - `Id > : mu + alpha_h`
  - `A > : mu + alpha_a`
- Notes: Set test_rate = 0 for the baseline. Granich projected incidence below 1/1000/yr within about 10 years of full implementation ([LSHTM](https://researchonline.lshtm.ac.uk/id/eprint/312973)). Commentators flagged its optimistic adherence assumptions. This differs from the existing "HIV stages" model.

---

### amr_bonhoeffer — sensitive and resistant strains with treatment and de novo resistance (single drug)
- Disease/area: bacterial infection with antibiotic treatment | Category: antimicrobial resistance
- Source: Bonhoeffer S., Lipsitch M., Levin B.R. (1997) "Evaluating treatment protocols to prevent antibiotic resistance", *PNAS* 94(22):12106–12111, doi:10.1073/pnas.94.22.12106, PMC23718 ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC23718)). The reduced equations in Appendix A2 were read: dx/dt = λ − dx − bxy + (r_w + fh(1−s))(1−ρ)y + r_rρy. Parameters are from Figure 2: λ = 100, d = 1, c = 1.5, h = 1, r_w = 0, s = 10⁻³, r_r = 0 or 0.1. **The paper gives no numeric value for b, so b below is illustrative.** Review context: Spicknall I.H. et al. (2013) *Am. J. Epidemiol.* 178:508–520, doi:10.1093/aje/kwt017 ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3736756)).
- Time unit: arbitrary (dimensionless; treat as weeks) ; simulate 100 time units
- Compartments:
  - `X` — uninfected hosts — 99
  - `Yw` — infected with wild-type (sensitive) — 1
  - `Yr` — infected with resistant strain — 0
- Parameters:
  - `lambda_in` = 100 (hosts/time) — entry of uninfected hosts — Bonhoeffer Fig. 2
  - `d_x` = 1 (1/time) — death/removal of uninfected — Bonhoeffer Fig. 2
  - `c_inf` = 1.5 (1/time) — death/removal of infected hosts — Bonhoeffer Fig. 2
  - `b_tr` = 0.05 (1/(host·time)) — mass-action transmission — illustrative (R0 = b·λ/(d·c) ≈ 3.3 untreated)
  - `h_max` = 1 (1/time) — cure rate under treatment — Bonhoeffer Fig. 2
  - `f_tr` = 0.5 (–) — fraction of infected treated — illustrative (the paper varies 0–1)
  - `s_res` = 0.001 (–) — fraction of treated wild-type infections that acquire resistance — Bonhoeffer Fig. 2
  - `r_w` = 0 (1/time) — spontaneous recovery, wild-type — Bonhoeffer Fig. 2
  - `r_r` = 0.1 (1/time) — spontaneous recovery, resistant (cost of resistance Δr = r_r − r_w) — Bonhoeffer Fig. 2B
- Derived variables:
  - `res_frac` = Yr/(Yw + Yr) — proportion of infections that are resistant (ρ)
- Flows:
  - ` > X : lambda_in` (absolute)
  - `X > : d_x`
  - `X > Yw : b_tr*Yw`
  - `X > Yr : b_tr*Yr`
  - `Yw > X : r_w + f_tr*h_max*(1 - s_res)`
  - `Yw > Yr : f_tr*h_max*s_res`
  - `Yr > X : r_r`
  - `Yw > : c_inf`
  - `Yr > : c_inf`
- Notes: Resistance takes over when f·h > Δr, as the ρ equation shows. The source's multi-drug extension (Yw, Ya, Yb, Yab with f_a, f_b, f_ab and probability q of double resistance) has 5 compartments, also fits the format, and can be added as `amr_bonhoeffer_multidrug` using the Fig. 1B equations quoted in the PMC text (cycling, 50:50 mixing, combination therapy). The paper's conclusion is that combination therapy is usually best.

---

### amr_hospital_lipsitch — hospital colonisation with sensitive and resistant bacteria and patient turnover
- Disease/area: nosocomial antibiotic-resistant bacteria (e.g. VRE, MRSA) | Category: antimicrobial resistance
- Source: Lipsitch M., Bergstrom C.T., Levin B.R. (2000) "The epidemiology of antibiotic resistance in hospitals: paradoxes and prescriptions", *PNAS* 97(4):1938–1943, doi:10.1073/pnas.97.4.1938 (from memory, not re-verified this session). **The structure is reconstructed from memory of that paper: uninfected X, colonised with sensitive S, colonised with resistant R, admissions carrying each type, and drug 1 (to which R is resistant) plus drug 2. The values are illustrative.**
- Time unit: days ; simulate 365 days
- Compartments:
  - `X` — uncolonised patients — 90
  - `Sc` — colonised with sensitive strain — 9
  - `Rc` — colonised with resistant strain — 1
- Parameters:
  - `mu_d` = 0.1 (1/day) — discharge rate = admission rate (10-day stay) — illustrative
  - `m_s` = 0.10 (–) — fraction admitted colonised with sensitive — illustrative
  - `m_r` = 0.01 (–) — fraction admitted colonised with resistant — illustrative
  - `beta_s` = 1.0 (1/day) — transmission, sensitive — illustrative
  - `beta_r` = 1.0 (1/day) — transmission, resistant — illustrative
  - `tau1` = 0.1 (1/day) — treatment rate with drug 1 (clears sensitive only) — illustrative
  - `tau2` = 0.05 (1/day) — treatment rate with drug 2 (clears both) — illustrative
  - `gamma` = 0.03 (1/day) — spontaneous clearance — illustrative
- Derived variables:
  - `prop_res` = Rc/(Sc + Rc)
- Flows:
  - ` > X : (1 - m_s - m_r)*mu_d*N` (absolute)
  - ` > Sc : m_s*mu_d*N` (absolute)
  - ` > Rc : m_r*mu_d*N` (absolute)
  - `X > Sc : beta_s*Sc/N`
  - `X > Rc : beta_r*Rc/N`
  - `Sc > X : tau1 + tau2 + gamma`
  - `Rc > X : tau2 + gamma`
  - `X > : mu_d`
  - `Sc > : mu_d`
  - `Rc > : mu_d`
- Notes: Total admissions equal total discharges, so N is constant. The paper's key insight is that in hospitals, with rapid turnover, resistance prevalence responds within weeks to changes in antibiotic use. Infection control (lowering beta) and reducing tau1 are the levers.

---

### two_strain_cross_immunity — two-strain SIR with partial cross-immunity (reduced susceptibility)
- Disease/area: influenza strains, dengue serotypes (without ADE), RSV | Category: multi-strain
- Source: Keeling M.J. & Rohani P. (2008), Chapter 4 ("Multi-pathogen / multi-host"), Section 4.1 and companion Program 4.1 (from memory; the companion site was unreachable). Original framework: Castillo-Chavez C., Hethcote H.W., Andreasen V., Levin S.A., Liu W.M. (1989) "Epidemiological models with age structure, proportionate mixing, and cross-immunity", *J. Math. Biol.* 27:233–258, doi:10.1007/BF00275810 (from memory). The values are illustrative.
- Time unit: days ; simulate 365*50 days (long-term dynamics)
- Compartments:
  - `S` — naive to both — 0.9*N0 (N0 = 100000)
  - `I1` — first infection, strain 1 — 10
  - `I2` — first infection, strain 2 — 10
  - `R1` — recovered from 1, susceptible to 2 — 0.05*N0
  - `R2` — recovered from 2, susceptible to 1 — 0.05*N0
  - `J1` — infected with strain 1 after recovering from 2 — 0
  - `J2` — infected with strain 2 after recovering from 1 — 0
  - `R12` — immune to both — 0
- Parameters:
  - `mu` = 1/(50*365) (1/day) — birth/death — illustrative
  - `beta1` = 0.5 (1/day) — strain 1 transmission — illustrative
  - `beta2` = 0.45 (1/day) — strain 2 transmission — illustrative
  - `gamma1` = 1/7 (1/day) — recovery, strain 1 — illustrative
  - `gamma2` = 1/7 (1/day) — recovery, strain 2 — illustrative
  - `sigma12` = 0.5 (–) — susceptibility to strain 2 after strain 1 (1 = no cross-immunity, 0 = full) — illustrative
  - `sigma21` = 0.5 (–) — susceptibility to strain 1 after strain 2 — illustrative
- Derived variables:
  - `lambda1` = beta1*(I1 + J1)/N
  - `lambda2` = beta2*(I2 + J2)/N
- Flows:
  - ` > S : mu*N` (absolute)
  - `S > I1 : lambda1`
  - `S > I2 : lambda2`
  - `I1 > R1 : gamma1`
  - `I2 > R2 : gamma2`
  - `R1 > J2 : sigma12*lambda2`
  - `R2 > J1 : sigma21*lambda1`
  - `J1 > R12 : gamma1`
  - `J2 > R12 : gamma2`
  - `S > : mu`
  - `I1 > : mu`
  - `I2 > : mu`
  - `R1 > : mu`
  - `R2 > : mu`
  - `J1 > : mu`
  - `J2 > : mu`
  - `R12 > : mu`
- Notes: There are 8 compartments and no simultaneous coinfection. Adding seasonal forcing to beta (beta*(1 + amp*cos(2*3.14159*t/365))) gives the out-of-phase oscillations discussed by Keeling & Rohani. Spell out the constant, because `pi` is disallowed. A variant with reduced transmissibility instead of reduced susceptibility puts the factor in lambda rather than on the R1 > J2 rate.

---

### superinfection_nowak_may — two-strain superinfection (virulence evolution)
- Disease/area: evolution of virulence; generic parasite | Category: multi-strain / superinfection
- Source: Nowak M.A. & May R.M. (1994) "Superinfection and the evolution of parasite virulence", *Proc. R. Soc. Lond. B* 255:81–89, doi:10.1098/rspb.1994.0012 ([search record](https://pmc.ncbi.nlm.nih.gov/articles/PMC1691217)). Also Anderson R.M. & May R.M. (1991) *Infectious Diseases of Humans*, Oxford University Press, ISBN 978-0-19-854040-3 (from memory). This is the two-strain special case. The values are illustrative.
- Time unit: arbitrary (years) ; simulate 200 years
- Compartments:
  - `X` — uninfected hosts — 990
  - `Y1` — infected with the more virulent strain 1 — 5
  - `Y2` — infected with the less virulent strain 2 — 5
- Parameters:
  - `Lambda` = 100 (hosts/yr) — host recruitment — illustrative
  - `mu` = 0.1 (1/yr) — natural host death — illustrative
  - `b1` = 0.005 (1/(host·yr)) — transmission, strain 1 — illustrative
  - `b2` = 0.004 (1/(host·yr)) — transmission, strain 2 — illustrative
  - `v1` = 0.5 (1/yr) — virulence (extra death), strain 1 — illustrative
  - `v2` = 0.2 (1/yr) — virulence, strain 2 — illustrative
  - `s_sup` = 1.0 (–) — superinfection coefficient (strain 1 superinfecting Y2 hosts) — illustrative
- Derived variables: none
- Flows:
  - ` > X : Lambda` (absolute)
  - `X > Y1 : b1*Y1`
  - `X > Y2 : b2*Y2`
  - `Y2 > Y1 : s_sup*b1*Y1`
  - `X > : mu`
  - `Y1 > : mu + v1`
  - `Y2 > : mu + v2`
- Notes: Without superinfection (s_sup = 0) the strain with the highest R0 = b·(Λ/μ)/(μ+v) excludes the other. With superinfection, the more virulent strain can coexist or win even with a lower R0, which is the paper's point that superinfection selects for higher virulence and polymorphism.

---

### sir_reinfection_partial — SIR with reinfection of recovered (partial immunity; reinfection threshold)
- Disease/area: TB, pertussis, RSV, malaria-like partial immunity | Category: reinfection
- Source: Gomes M.G.M., White L.J., Medley G.F. (2004) "Infection, reinfection, and vaccination under suboptimal immune protection: epidemiological perspectives", *J. Theor. Biol.* 228(4):539–549, doi:10.1016/j.jtbi.2004.02.015 ([Gulbenkian repository](https://arca.igc.gulbenkian.pt/entities/publication/7e81bd2c-5bd7-47cf-8601-830a31e3552a/full)). Critique: Breban R. & Blower S. (2005) "The reinfection threshold does not exist", *J. Theor. Biol.* (same journal). The values are illustrative.
- Time unit: years ; simulate 200 years
- Compartments:
  - `S` — susceptible — 9990
  - `I` — infectious — 10
  - `R` — recovered, partially immune — 0
  - `V` — vaccinated at birth (same partial protection as R) — 0
- Parameters:
  - `mu` = 1/70 (1/yr) — birth/death — illustrative
  - `gamma` = 52 (1/yr) — recovery (1 week) — illustrative
  - `beta` = 260 (1/yr) — transmission (R0 ≈ 5) — illustrative
  - `sigma_r` = 0.25 (–) — relative susceptibility to reinfection (reinfection threshold at R0 = 1/sigma_r = 4) — illustrative
  - `p_nb` = 0 (–) — newborn vaccine coverage (set to 0.9 to show vaccine failure above the reinfection threshold) — illustrative
- Derived variables:
  - `lambda` = beta*I/N
- Flows:
  - ` > S : (1 - p_nb)*mu*N` (absolute)
  - ` > V : p_nb*mu*N` (absolute)
  - `S > I : lambda`
  - `R > I : sigma_r*lambda`
  - `V > I : sigma_r*lambda`
  - `I > R : gamma`
  - `S > : mu`
  - `I > : mu`
  - `R > : mu`
  - `V > : mu`
- Notes: Gomes et al. report that above the reinfection threshold reinfection dominates, prevalence is high and vaccination fails to protect ([abstract via repository](https://arca.igc.gulbenkian.pt/entities/publication/7e81bd2c-5bd7-47cf-8601-830a31e3552a/full)). Whether the threshold is a true bifurcation was disputed by Breban & Blower (2005). This differs from SIRS because immunity is partial rather than waning.

---

### sis_core_group_gonorrhea — two-risk-group SIS with screening and treatment (Hethcote & Yorke core group)
- Disease/area: gonorrhoea and other bacterial STIs | Category: risk structure / test-and-treat
- Source: Hethcote H.W. & Yorke J.A. (1984) *Gonorrhea Transmission Dynamics and Control*, Lecture Notes in Biomathematics 56, Springer, doi:10.1007/978-3-662-07544-9 (from memory, not re-verified). Also Anderson & May (1991), Chapter 11 (from memory). The core-group concept (a small, highly active group sustaining transmission) is from Hethcote & Yorke. **The values are illustrative.**
- Time unit: years ; simulate 20 years
- Compartments:
  - `Sh` — susceptible, core (high activity) — 1980
  - `Ih` — infected, core — 20
  - `Sl` — susceptible, non-core — 97900
  - `Il` — infected, non-core — 100
- Parameters:
  - `beta_hh` = 30 (1/yr) — core-to-core transmission — illustrative
  - `beta_hl` = 1.0 (1/yr) — infection of core by non-core — illustrative
  - `beta_lh` = 0.6 (1/yr) — infection of non-core by core — illustrative
  - `beta_ll` = 0.5 (1/yr) — non-core to non-core — illustrative
  - `gamma` = 4 (1/yr) — natural clearance or symptomatic treatment (3-month duration) — illustrative
  - `scr_h` = 2 (1/yr) — screening-and-treatment rate in the core — illustrative
  - `scr_l` = 0.2 (1/yr) — screening-and-treatment rate in the non-core — illustrative
- Derived variables:
  - `Nh` = Sh + Ih
  - `Nl` = Sl + Il
  - `lambda_h` = beta_hh*Ih/Nh + beta_hl*Il/Nl
  - `lambda_l` = beta_lh*Ih/Nh + beta_ll*Il/Nl
- Flows:
  - `Sh > Ih : lambda_h`
  - `Sl > Il : lambda_l`
  - `Ih > Sh : gamma + scr_h`
  - `Il > Sl : gamma + scr_l`
- Notes: Group sizes are constant (no demography). Check the next-generation R0, because the transmission values may need tuning for an endemic state. Use it to show that targeting screening to the core (raising scr_h) reduces prevalence far more than the same effort spread uniformly, which is Hethcote & Yorke's key policy message. This differs from the existing "SIR two age groups" model.

---

### sir_vertical_transmission — SIR with vertical (mother-to-child) transmission
- Disease/area: hepatitis B, rubella (congenital), HTLV, Chagas | Category: classic textbook
- Source: Busenberg S. & Cooke K. (1993) *Vertically Transmitted Diseases: Models and Dynamics*, Biomathematics 23, Springer, doi:10.1007/978-3-642-75301-5 (from memory). Also Brauer & Castillo-Chavez (2012) and Hethcote 2000 (from memory). The values are illustrative.
- Time unit: years ; simulate 200 years
- Compartments:
  - `S` — susceptible — 9900
  - `I` — infectious — 100
  - `R` — recovered/immune — 0
- Parameters:
  - `b_birth` = 1/60 (1/yr) — per-capita birth rate — illustrative
  - `mu` = 1/60 (1/yr) — natural death — illustrative
  - `q_v` = 0.3 (–) — probability a child of an infected mother is born infected — illustrative
  - `beta` = 20 (1/yr) — horizontal transmission — illustrative
  - `gamma` = 10 (1/yr) — recovery — illustrative (a chronic HBV-like case would use a small gamma)
- Derived variables:
  - `lambda` = beta*I/N
- Flows:
  - ` > S : b_birth*(S + R + (1 - q_v)*I)` (absolute)
  - ` > I : b_birth*q_v*I` (absolute)
  - `S > I : lambda`
  - `I > R : gamma`
  - `S > : mu`
  - `I > : mu`
  - `R > : mu`
- Notes: Vertical transmission adds q_v·b_birth to the infected's "reproductive output", so R0 = (beta + q_v*b_birth)/(gamma + mu). It matters most for chronic infections.

---

### sidr_diagnosis — acute infection with diagnosis, isolation-at-home and treatment (no carriers)
- Disease/area: generic acute bacterial or viral infection with case finding (e.g. diphtheria, scabies, early COVID testing) | Category: diagnosis / test-and-treat
- Source: **Illustrative structure** in the spirit of the treatment and control chapters of Vynnycky E. & White R.G. (2010) *An Introduction to Infectious Disease Modelling*, Oxford University Press, ISBN 978-0-19-856576-5, and Martcheva M. (2015) *An Introduction to Mathematical Epidemiology*, Springer, doi:10.1007/978-1-4899-7612-3 (both from memory). No primary-source values.
- Time unit: days ; simulate 200 days
- Compartments:
  - `S` — susceptible — 9990
  - `Iu` — infectious, undiagnosed — 10
  - `Id` — diagnosed and treated, reduced contacts — 0
  - `R` — recovered — 0
- Parameters:
  - `beta` = 0.4 (1/day) — transmission — illustrative
  - `delta_d` = 0.2 (1/day) — diagnosis rate (mean 5 days to diagnosis) — illustrative
  - `kappa_d` = 0.2 (–) — relative infectiousness after diagnosis (isolation plus treatment) — illustrative
  - `gamma_u` = 1/10 (1/day) — recovery, undiagnosed — illustrative
  - `gamma_d` = 1/3 (1/day) — recovery, diagnosed and treated — illustrative
- Derived variables:
  - `lambda` = beta*(Iu + kappa_d*Id)/N
  - `R0` = beta/(gamma_u + delta_d) + beta*delta_d*kappa_d/((gamma_u + delta_d)*gamma_d)
- Flows:
  - `S > Iu : lambda`
  - `Iu > Id : delta_d`
  - `Iu > R : gamma_u`
  - `Id > R : gamma_d`
- Notes: R0 is about 1.53 with these values, against 4 without diagnosis. Differs from the existing SIQR because diagnosis both shortens infectiousness (treatment) and only partially isolates. It is also a minimal template for "test-and-treat" scenarios: vary delta_d.

### Inferences
- All 17 models are expressible in the required format with at most 9 compartments. The largest are two_strain_cross_immunity (8) and seiqhr_sars (7, including the not-in-N death accumulator).
- Every model conserves N except through explicit births and deaths, admissions and discharges, or the not-in-N accumulator `Dd`. In closed models (sitr_treatment, sir_saturated_incidence, sir_awareness, sis_core_group_gonorrhea, sidr_diagnosis) every flow is internal. In models with demography, births balance deaths when the absolute birth term equals mu*N (svirs_waning_vaccine, siv_backward_bifurcation, sirv_booster, two_strain_cross_immunity, sir_reinfection_partial, amr_hospital_lipsitch). In sir_vertical_transmission births equal deaths when b_birth = mu.
- Models showing backward bifurcation or bistability (sir_saturated_treatment, svirs_waning_vaccine, siv_backward_bifurcation) are good teaching additions, but the illustrative values may need tuning to land in the bistable region. This should be checked numerically when implemented.
- The amr_bonhoeffer multi-drug version can be implemented directly from the PMC equations if a 13th–18th model is wanted.

## What is the primary source for each model, and where do the parameter values come from?

### Takeaway
Every model has a named primary source, given in its "Source" field above. Parameter values come directly from the paper only for amr_bonhoeffer (Figure 2, read from PMC). The Gumel 2004 and Granich 2009 values are recalled and need verifying against the original tables. All other values are illustrative and are labelled as such.

### Cited Findings
- Funk et al. 2009, PNAS 106(16):6872–6877, doi:10.1073/pnas.0810762106: in well-mixed populations, awareness reduces outbreak size but not the epidemic threshold — [LSHTM](https://researchonline.lshtm.ac.uk/id/eprint/611254), [Omics DI](https://omicsdi.org/dataset/biostudies-other/S-EPMC2672559)
- Gumel et al. 2004, Proc R Soc B 271:2223–2232, doi:10.1098/rspb.2004.2800: six-class model (S, asymptomatic, quarantined, symptomatic, isolated, recovered). Optimal isolation beats sub-optimal isolation plus quarantine — [PMC1691853](https://pmc.ncbi.nlm.nih.gov/articles/PMC1691853/)
- Bonhoeffer, Lipsitch & Levin 1997, PNAS 94:12106–12111, doi:10.1073/pnas.94.22.12106: equations and Figure 2 and 3 values (λ = 100, d = 1, c = 1.5, h = 1, r_w = 0, s = 10⁻³, r_r ∈ {0, 0.1}; b unspecified). Combination therapy is generally best — [PMC23718](https://pmc.ncbi.nlm.nih.gov/articles/PMC23718)
- Wang & Ruan 2004 (constant removal rate), JMAA 291:775–793, doi:10.1016/j.jmaa.2003.11.043 — [MaRDI](https://portal.mardi4nfdi.de/wiki/Publication:1433379). Saturated treatment h(I) = βI/(1+αI) — [Springer, doi:10.1007/s12190-009-0315-9](https://link.springer.com/doi/10.1007/s12190-009-0315-9)
- Capasso & Serio 1978, Math Biosci 42:43–61, doi:10.1016/0025-5564(78)90006-8 — [MaRDI](https://portal.mardi4nfdi.de/wiki/Publication:1254217)
- Gomes, White & Medley 2004, J Theor Biol 228:539–549, doi:10.1016/j.jtbi.2004.02.015. Reinfection threshold, with vaccination ineffective above it — [repository](https://arca.igc.gulbenkian.pt/entities/publication/7e81bd2c-5bd7-47cf-8601-830a31e3552a/full)
- Nowak & May 1994, Proc R Soc B 255:81–89, doi:10.1098/rspb.1994.0012. Superinfection selects for higher virulence and polymorphism — [PMC citing record](https://pmc.ncbi.nlm.nih.gov/articles/PMC1691217)
- Kribs-Zaleta & Velasco-Hernández 2000, Math Biosci 164:183–201, doi:10.1016/S0025-5564(00)00003-1. A two-dimensional SIS model with vaccination shows backward bifurcation — [MaRDI](https://portal.mardi4nfdi.de/wiki/Publication:1572822)
- Granich et al. 2009, Lancet 373:48–57, doi:10.1016/S0140-6736(08)61697-9. Annual universal testing with immediate ART was projected to bring incidence below 1/1000/yr within 10 years — [LSHTM](https://researchonline.lshtm.ac.uk/id/eprint/312973), [NATAP](https://natap.org/2010/HIV/032010_01.htm)
- Spicknall et al. 2013, Am J Epidemiol 178:508–520, doi:10.1093/aje/kwt017. All AMR transmission models share sensitive–resistant competition, and model structure is rarely tested — [PMC3736756](https://pmc.ncbi.nlm.nih.gov/articles/PMC3736756)

### Inferences
- Before release, pull the Gumel 2004 Table 1 (PMC PDF or the electronic appendix) to replace the recalled SARS values, and confirm Granich's ART dropout and infectiousness-reduction figures.
- The Lipsitch 2000 hospital model structure and the Arino 2003, Zhang & Liu 2008, Kiss 2010, Castillo-Chavez 1989, Hethcote & Yorke 1984 and Busenberg & Cooke 1993 citations come from memory. Verify their DOIs before publishing the catalogue.

### Gaps
- The Keeling & Rohani companion programs (modelinginfectiousdiseases.org and the Warwick mirror) were unreachable (connection refused or DNS failure). Program-specific default parameters, for example for Program 4.1 (two strains), could not be confirmed.
- The full text of Gumel et al. 2004 was blocked by a captcha, so its parameter table was not verified.
- The transmission coefficient b in Bonhoeffer 1997 has no numeric value in the text, so it is illustrative here.
- No primary source was retrieved for the booster model (sirv_booster) or the diagnosis model (sidr_diagnosis). Both are explicitly illustrative constructions.
- The exact equations of the Kribs-Zaleta & Velasco-Hernández and Lipsitch 2000 models were not read in full this session. Their flow structures are standard reconstructions.
