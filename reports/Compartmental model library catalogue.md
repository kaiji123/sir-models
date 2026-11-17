# Eighty-nine compartmental models, ready for code generation

The research produced **89 deterministic ODE model specifications across six domains**. Each one can be written in a generator's compartments + flows + parameters + derived-variables grammar (`FROM > TO : rate`, per-capita unless marked absolute) with no more than 14 compartments. After merging three near-duplicates, about 86 are distinct. The structures are sound and mostly follow named primary sources. The parameter values are the weak point. **Only about a dozen models have parameter sets read directly from the primary paper this session**: Bonhoeffer 1997 (antimicrobial resistance, AMR), Arino 2008 (influenza antivirals), SIDARTHE (COVID-19, arXiv v1), Chowell 2004 (SARS), Chitnis 2008 (malaria), Bowman 2005 (West Nile), Gao 2016 (Zika), Yakob & Clements 2013 (chikungunya), Whittles 2017 (gonorrhoea AMR), Martin 2013 (hepatitis C in people who inject drugs, PWID), Codeço 2001 (cholera) and the two transmission rates of Zinsstag 2009 (dog rabies). A second tier has some values verified: Pitzer 2014 (typhoid), Simmons 2013 (norovirus), syphilis durations from Garnett 1997 via Feldman & Mishra 2019, the chlamydia asymptomatic duration from Althaus 2010, and the RSV rates traced to Weber 2001. Every other model uses recalled or illustrative numbers. The library should ship them with a provenance flag rather than present them as published estimates. The catalogue also turned up a handful of specification bugs and source conflicts that need fixing before code generation (section 8).

## 1. How the provenance labels work, and why they matter more than the structures

Each researcher tagged every parameter with one of three levels, which this report standardises as **V** (verified against the primary text or a fetched secondary source this session), **R** (recalled from the literature, not re-checked) and **I** (illustrative, chosen so the model runs with plausible dynamics). A model is labelled by its weakest load-bearing value. Model structure is generally more reliable than the numbers. Where a published model was age-structured, stochastic, spatial or individual-based, the specification is a simplified deterministic, single-age analogue. It reproduces the source's qualitative mechanism, not its calibration. Pitzer's typhoid and rotavirus models, Hethcote's pertussis model, Van Effelterre's hepatitis A model, Whittles's gonorrhoea model and Kaplan's smallpox model are all in this position.

Some conventions are shared across the specifications, and the generator should enforce them:

- **Identifier names.** Identifiers avoid R base and stats function names and reserved words (`beta`, `gamma`, `sigma`, `kappa`, `q`, `c`, `C`, `D`, `I`, `T`, `F`, `t`, `pi`, `df`, `dt`). Some notes still use bare `beta`, `gamma` or `I`. These work inside a deSolve `with()` block but should be renamed if the validator is strict, for example to `beta_`, `g_rec` or `Inf_h`. The constant 2π is spelled out as `6.2831853`.
- **Accumulators.** Cumulative counters (deaths, case counts, congenital rubella syndrome (CRS) births, vaccine doses) are marked "not in N" and must be excluded from force-of-infection denominators.
- **Non-human compartments.** Vectors, reservoirs, environmental pathogen loads and within-host cell or virion pools are also marked "not in N".
- **Time-varying inputs.** Interventions, seasonality and pulsed mass drug administration (MDA) use `ifelse`, `floor` and products of comparisons. A pulse is applied as the log-kill rate −log(1 − coverage·efficacy)/window. When pulses are used, the solver step must be limited (`hmax` ≤ window).
- **Fixed delays.** Fixed delays become Erlang stages. Examples are the extrinsic incubation period (EIP) in vectors and incubation in COVID-19.

The shorthand below writes death flows as `X,Y,Z > : mu`, meaning one per-capita exit flow from each listed compartment, and marks absolute flows with **[abs]**.

## 2. Classic, intervention, behaviour and AMR models: 17 specifications, one with source-read values

This group extends the textbook SIR/SEIR family with treatment, saturation, awareness, imperfect vaccines, quarantine, test-and-treat, resistance, multi-strain competition and vertical transmission. Only **Bonhoeffer, Lipsitch & Levin 1997** was checked number-by-number. Its reduced equations and Figure 2 values were read from PMC (λ = 100, d = 1, c = 1.5, h = 1, r_w = 0, s = 10⁻³, r_r = 0 or 0.1), but **the paper gives no numeric transmission coefficient b** ([PMC23718](https://pmc.ncbi.nlm.nih.gov/articles/PMC23718)). The Gumel 2004 SARS values and the Granich 2009 ART figures are recalled. The full Gumel text was blocked by a captcha, although its six-class structure (S, asymptomatic, quarantined, symptomatic, isolated, recovered) was confirmed ([PMC1691853](https://pmc.ncbi.nlm.nih.gov/articles/PMC1691853/)).

| ID | Disease / area | Compartments | Force of infection | Typical parameters (units) | Primary source | Prov. |
|---|---|---|---|---|---|---|
| sitr_treatment | Influenza with antivirals | S,E,I,Tr,R | beta·(I+delta·Tr)/N | beta 0.5/d; kappa 1/1.5/d; tau 0.3/d; alpha 0.25/d; eta 0.4/d; delta 0.5; R0≈1.25 | Brauer 2008, *Math Epidemiology* LNM 1945 (R) | I |
| sir_saturated_treatment | Generic, capacity-limited treatment | S,I,R | beta·I (mass action) | Lambda 0.391/d [abs]; mu 1/(70·365)/d; beta 3e-5/(person·d); gamma 0.1/d; d_dis 0.001/d; r_max 0.1/d; alpha_t 0.01/person; R0≈1.5 | Wang & Ruan 2004 ([MaRDI](https://portal.mardi4nfdi.de/wiki/Publication:1433379)); Zhonghua & Yaohong 2009 ([Springer](https://link.springer.com/doi/10.1007/s12190-009-0315-9)) | I |
| sir_saturated_incidence | Cholera (Bari 1973) / generic | S,I,R | beta·I/(1+alpha_s·I) | beta 5e-5/(person·d); alpha_s 0.002/person; gamma 0.2/d | Capasso & Serio 1978 ([MaRDI](https://portal.mardi4nfdi.de/wiki/Publication:1254217)) | I |
| sir_awareness | Emerging infection, behaviour | S,Sa,I,R | beta·I/N (aware: ×sigma_a) | beta 0.4/d; gamma 0.2/d; phi_a 2/d; omega_f 0.05/d; sigma_a 0.3 | Funk et al. 2009 PNAS ([LSHTM](https://researchonline.lshtm.ac.uk/id/eprint/611254)); simplified mean-field version | I |
| svirs_waning_vaccine | Influenza/pertussis-like | S,V,I,R | beta·I/N | mu 1/70/yr; p_nb 0.5; nu_v 0.1/yr; theta_v 0.2/yr; sigma_v 0.2; beta 100/yr; gamma 365/7/yr; omega_r 0.1/yr | Arino, McCluskey & van den Driessche 2003 SIAM J Appl Math (R) | I |
| siv_backward_bifurcation | SIS with imperfect, waning vaccine | S,I,V | beta·I/N | mu 0.02/yr; beta 20/yr; gamma 10/yr; phi_v 1/yr; theta_v 0.1/yr; sigma_v 0.3 | Kribs-Zaleta & Velasco-Hernández 2000 ([MaRDI](https://portal.mardi4nfdi.de/wiki/Publication:1572822)) | I |
| sirv_booster | Pertussis/COVID/HPV-like boosters | S,V1,W,V2,I,R | beta·I/N | p_nb 0.9; omega1 1/4/yr; omega2 1/10/yr; omega_r 1/20/yr; b_rate 0.2/yr; eps1=eps2 0.1; eps_w 0.5; beta 300/yr; gamma 365/21/yr | Illustrative (Hethcote 2000; Keeling & Rohani 2008) | I |
| seiqhr_sars | SARS 2003, quarantine + isolation | S,E,Q,I,H,R,Dd* | beta·(I+eps_E·E+eps_Q·Q+eps_H·H)/N | beta 0.2/d; eps_E 0.3; eps_H 0.36; gamma1 0.1; gamma2 0.5; kappa1 0.1; kappa2 0.125; sigma1 0.0337; sigma2 0.0386; d1 0.0068; d2 0.0045 (all /d) | Gumel et al. 2004 Proc R Soc B ([PMC1691853](https://pmc.ncbi.nlm.nih.gov/articles/PMC1691853/)) | R |
| hiv_test_and_treat | HIV, South Africa, universal test-and-treat (UTT) | S,Iu,Id,A | beta·(Iu+eps_d·Id+eps_A·A)/N | Lambda 250000/yr [abs]; mu 1/35/yr; beta 0.5/yr; eps_A 0.01; eps_d 0.8; alpha_h 0.1/yr; test_rate 1/yr; link_art 0.9; art_init 0.5/yr; dropout 0.015/yr | Granich et al. 2009 Lancet ([LSHTM](https://researchonline.lshtm.ac.uk/id/eprint/312973)) | R/I |
| amr_bonhoeffer | Bacterial infection, single-drug resistance | X,Yw,Yr | b_tr·Y (mass action, per strain) | lambda_in 100 [abs]; d_x 1; c_inf 1.5; h_max 1; s_res 1e-3; r_w 0; r_r 0.1 (V); b_tr 0.05 (I); f_tr 0.5 (I) | Bonhoeffer et al. 1997 PNAS ([PMC23718](https://pmc.ncbi.nlm.nih.gov/articles/PMC23718)) | **V** (except b, f) |
| amr_hospital_lipsitch | Nosocomial VRE/MRSA | X,Sc,Rc | beta_s·Sc/N; beta_r·Rc/N | mu_d 0.1/d; m_s 0.10; m_r 0.01; beta_s=beta_r 1/d; tau1 0.1/d; tau2 0.05/d; gamma 0.03/d | Lipsitch, Bergstrom & Levin 2000 PNAS (structure R) | I |
| two_strain_cross_immunity | Influenza strains, dengue without ADE, RSV | S,I1,I2,R1,R2,J1,J2,R12 | beta_k·(I_k+J_k)/N | beta1 0.5/d; beta2 0.45/d; gamma 1/7/d; sigma12=sigma21 0.5; mu 1/(50·365)/d | Keeling & Rohani 2008 ch.4; Castillo-Chavez et al. 1989 (R) | I |
| superinfection_nowak_may | Virulence evolution | X,Y1,Y2 | b_k·Y_k (mass action) | Lambda 100/yr; mu 0.1/yr; b1 0.005, b2 0.004 /(host·yr); v1 0.5, v2 0.2 /yr; s_sup 1 | Nowak & May 1994 Proc R Soc B ([record](https://pmc.ncbi.nlm.nih.gov/articles/PMC1691217)) | I |
| sir_reinfection_partial | TB/pertussis/RSV partial immunity | S,I,R,V | beta·I/N (R,V: ×sigma_r) | mu 1/70/yr; gamma 52/yr; beta 260/yr (R0≈5); sigma_r 0.25; p_nb 0–0.9 | Gomes, White & Medley 2004 JTB ([repository](https://arca.igc.gulbenkian.pt/entities/publication/7e81bd2c-5bd7-47cf-8601-830a31e3552a/full)) | I |
| sis_core_group_gonorrhea | Gonorrhoea, core group + screening | Sh,Ih,Sl,Il | beta_xy·I_y/N_y summed | beta_hh 30, beta_hl 1, beta_lh 0.6, beta_ll 0.5 /yr; gamma 4/yr; scr_h 2/yr; scr_l 0.2/yr | Hethcote & Yorke 1984 LNBiomath 56 (R) | I |
| sir_vertical_transmission | HBV, rubella, HTLV, Chagas | S,I,R | beta·I/N | b_birth=mu 1/60/yr; q_v 0.3; beta 20/yr; gamma 10/yr; R0=(beta+q_v·b_birth)/(gamma+mu) | Busenberg & Cooke 1993 (R) | I |
| sidr_diagnosis | Generic case finding / test-and-treat | S,Iu,Id,R | beta·(Iu+kappa_d·Id)/N | beta 0.4/d; delta_d 0.2/d; kappa_d 0.2; gamma_u 0.1/d; gamma_d 1/3/d; R0≈1.5 (4 without diagnosis) | Illustrative (Vynnycky & White 2010; Martcheva 2015) | I |

\* `Dd` is an accumulator and is not in N.

```text
sitr_treatment:        S>E:lambda; E>I:kappa; I>Tr:tau; I>R:alpha; Tr>R:eta
sir_saturated_treatment: >S:Lambda[abs]; S>I:beta*I; I>R:gamma; I>R:r_max/(1+alpha_t*I); S,R>:mu; I>:mu+d_dis
sir_saturated_incidence: S>I:beta*I/(1+alpha_s*I); I>R:gamma
sir_awareness:         S>Sa:phi_a*I/N; Sa>S:omega_f; S>I:lambda; Sa>I:sigma_a*lambda; I>R:gamma
svirs_waning_vaccine:  >S:(1-p_nb)*mu*N[abs]; >V:p_nb*mu*N[abs]; S>V:nu_v; V>S:theta_v; S>I:lambda; V>I:sigma_v*lambda; I>R:gamma; R>S:omega_r; S,V,I,R>:mu
siv_backward_bifurcation: >S:mu*N[abs]; S>I:lambda; I>S:gamma; S>V:phi_v; V>S:theta_v; V>I:sigma_v*lambda; S,I,V>:mu
sirv_booster:          >S:(1-p_nb)*mu*N[abs]; >V1:p_nb*mu*N[abs]; V1>W:omega1; V2>W:omega2; R>W:omega_r; W>V2:b_rate; S>I:lambda; V1>I:eps1*lambda; V2>I:eps2*lambda; W>I:eps_w*lambda; I>R:gamma; all>:mu
seiqhr_sars:           >S:Pi_b[abs]; S>E:lambda; E>Q:gamma1; E>I:kappa1; Q>H:kappa2; I>H:gamma2; I>R:sigma1; H>R:sigma2; I>Dd:d1; H>Dd:d2; S,E,Q,I,H,R>:mu
hiv_test_and_treat:    >S:Lambda[abs]; S>Iu:lambda; Iu>A:test_rate*link_art; Iu>Id:test_rate*(1-link_art); Id>A:art_init; A>Id:dropout; S>:mu; Iu,Id>:mu+alpha_h; A>:mu+alpha_a
amr_bonhoeffer:        >X:lambda_in[abs]; X>:d_x; X>Yw:b_tr*Yw; X>Yr:b_tr*Yr; Yw>X:r_w+f_tr*h_max*(1-s_res); Yw>Yr:f_tr*h_max*s_res; Yr>X:r_r; Yw,Yr>:c_inf
amr_hospital_lipsitch: >X:(1-m_s-m_r)*mu_d*N[abs]; >Sc:m_s*mu_d*N[abs]; >Rc:m_r*mu_d*N[abs]; X>Sc:beta_s*Sc/N; X>Rc:beta_r*Rc/N; Sc>X:tau1+tau2+gamma; Rc>X:tau2+gamma; X,Sc,Rc>:mu_d
two_strain_cross_immunity: >S:mu*N[abs]; S>I1:lambda1; S>I2:lambda2; I1>R1:gamma1; I2>R2:gamma2; R1>J2:sigma12*lambda2; R2>J1:sigma21*lambda1; J1>R12:gamma1; J2>R12:gamma2; all>:mu
superinfection_nowak_may: >X:Lambda[abs]; X>Y1:b1*Y1; X>Y2:b2*Y2; Y2>Y1:s_sup*b1*Y1; X>:mu; Y1>:mu+v1; Y2>:mu+v2
sir_reinfection_partial: >S:(1-p_nb)*mu*N[abs]; >V:p_nb*mu*N[abs]; S>I:lambda; R>I:sigma_r*lambda; V>I:sigma_r*lambda; I>R:gamma; S,I,R,V>:mu
sis_core_group_gonorrhea: Sh>Ih:lambda_h; Sl>Il:lambda_l; Ih>Sh:gamma+scr_h; Il>Sl:gamma+scr_l
sir_vertical_transmission: >S:b_birth*(S+R+(1-q_v)*I)[abs]; >I:b_birth*q_v*I[abs]; S>I:lambda; I>R:gamma; S,I,R>:mu
sidr_diagnosis:        S>Iu:lambda; Iu>Id:delta_d; Iu>R:gamma_u; Id>R:gamma_d
```

Several of these models have a teaching result to reproduce.

- **Awareness (Funk et al. 2009).** In a well-mixed population, awareness shrinks the outbreak without moving the epidemic threshold. This is the paper's finding ([LSHTM](https://researchonline.lshtm.ac.uk/id/eprint/611254)).
- **Reinfection threshold (Gomes et al. 2004).** Above the threshold, vaccination stops protecting the population ([repository](https://arca.igc.gulbenkian.pt/entities/publication/7e81bd2c-5bd7-47cf-8601-830a31e3552a/full)). Breban & Blower (2005) dispute whether this is a true bifurcation.
- **Universal test-and-treat (Granich et al. 2009).** Annual testing with immediate ART was projected to cut incidence below 1/1000/yr within about 10 years ([LSHTM](https://researchonline.lshtm.ac.uk/id/eprint/312973)).
- **Backward bifurcation.** The saturated-treatment, SVIRS and SIV models can all show a backward bifurcation. Their illustrative values may need tuning before they land in the bistable region.

## 3. Respiratory and childhood infections: 16 specifications, four with source-read parameters

The respiratory set is the strongest on provenance.

- **Arino et al. 2008** take their pandemic-influenza values from Longini et al. 2004 and convert them to rates. Treatment of 80% within a day becomes φ = −ln(0.2) = 1.61/day, and b0 = 0.581 gives R0 = 1.35 under standard incidence with 40% of symptomatic cases withdrawing from mixing ([Arino et al. 2008 PDF](https://julien-arino.github.io/assets/pdf/papers/2008_ArinoBrauerPvdDWatmoughWu-JTB253.pdf)). One caveat: text extraction misaligned Table 1. The relapse rate θ = 0.0179 and the relative factors σ_S = 0.7 and σ = 0.2 were assigned by row order and checked against the Table 2–3 captions.
- **SIDARTHE.** The day-1 Italian values give R0 = 2.38, and the day-4 measures cut it to 1.66 ([Giordano et al., arXiv:2003.09861](https://arxiv.org/abs/2003.09861)). These are the arXiv v1 values. The *Nature Medicine* version is recalled as differing (μ = 0.017, τ = 0.01, unchecked). The Greek labels on the day-22 scenario values were lost in extraction.
- **Chowell et al. EID 2004.** The SARS SEIJR values were read from the article ([Chowell et al. 2004 EID](https://wwwnc.cdc.gov/eid/article/10/7/03-0647_article)). Its reported median R0 values (1.10 Hong Kong, 1.17 Singapore, 0.58 Toronto) look low next to the usual R0 ≈ 3 for SARS. They probably describe reproduction under control, so check the definition before using them as validation targets.
- **Weber RSV rates.** Latency of about 4 days, infectiousness of about 10 days and immunity of about 200 days trace to Weber 2001. Later work calls the 200-day immunity figure weakly supported ([ANU thesis](https://openresearch-repository.anu.edu.au/server/api/core/bitstreams/e3e1d619-dc5a-4d65-9747-7e31ca0ab4ad/content)).

| ID | Disease / area | Compartments | Force of infection | Typical parameters (units) | Primary source | Prov. |
|---|---|---|---|---|---|---|
| flu_arino2008_antiviral | Pandemic influenza, antivirals/prophylaxis | S,ST,L,LT,Ip,IpT,A,AT,R,Dead* | b0·Lam/N_M, where N_M = N−q_w(Ip+IpT) | k_lat 0.526; alpha_r=eta_r 0.244; alpha_T=eta_T 0.323 (/d); p_sym 0.667; tau_av 0.4; f_rec 0.98; delta_a 0.5; phi 1.61/d; theta 0.0179/d; sig_S 0.7; sig 0.2; q_w 0.4; b0 0.581/d | Arino et al. 2008 JTB ([PDF](https://julien-arino.github.io/assets/pdf/papers/2008_ArinoBrauerPvdDWatmoughWu-JTB253.pdf)) | **V** |
| flu_two_strain_cross | Co-circulating influenza strains | S,I1,I2,R1,R2,J1,J2,R12 | b_k·(I_k+J_k)/N | b1 0.5, b2 0.45, g_rec 0.25 (/d); cx 0.3; mu 3.9e-5/d; omega 1/730/d | Castillo-Chavez 1989; Andreasen, Lin & Levin 1997 (R) | I |
| flu_seasonal_vax | Seasonal influenza SIRS with vaccination | S,Ip,R,V | R0_base·g_rec·(1+amp·cos(2π(t−t_peak)/365))·Ip/N | R0_base 1.3; g_rec 0.25/d; amp 0.2; omega_R 1/730/d; omega_V 1/240/d; cov 0.45; ve 0.5; campaign day 274 for 60 d | Keeling & Rohani 2008 ch.5; Biggerstaff 2014 R≈1.28 (R) | I |
| covid_sidarthe | COVID-19, Italy 2020 | S,Iu,Dg,A,Rd,Th,H,Ex* | (alpha·Iu+beta·Dg+gamma·A+delta·Rd)/Npop | alpha 0.570→0.422 (t≥3); beta=delta 0.011→0.0057; gamma 0.456→0.285; epsilon 0.171→0.143 (t≥11); theta 0.371; zeta=eta 0.125; mu 0.012; nu 0.027; tau 0.003; lambda=rho 0.034; kappa=xi=sigma 0.017 (/d) | Giordano et al. 2020 Nat Med; [arXiv:2003.09861](https://arxiv.org/abs/2003.09861) | **V** (arXiv v1) |
| covid_kucharski2020 | Early Wuhan SARS-CoV-2 | S,E1,E2,I1,I2,Rm,Qr*,Cr* | Rt/inf_per·(I1+I2)/N | inc_per 5.2 d; inf_per 2.9 d; rep_delay 6.1 d; R 2.35→1.05 at day 23 | Kucharski et al. 2020 Lancet Infect Dis (R) | R |
| covid_vax_waning | COVID-19 with waning, boosting, hospitalisation | S,V,Wv,E,Ip,Hc,R,Dead* | R0_eff·g_rec·Ip/N | R0_eff 3; sig_lat 1/3; g_rec 1/5; h_rec 1/8 (/d); p_hosp 0.02; p_die 0.15; ve_V 0.8; ve_W 0.3; wane_V 1/120/d; wane_R 1/180/d | Illustrative; waning pattern after Andrews et al. 2022 NEJM (R) | I |
| rsv_weber2001_graded | RSV, seasonal, graded reinfection | M,S1,E1,I1,R1,S2,E2,I2,R2 | b0·(1+b1·cos(2πt+phase))·(I1+I2)/N | b0 60/yr (I); b1 0.3; sig_lat 91.25/yr; nu 36.5/yr; gam_w 1.825/yr (V via secondary); xi_m 2/yr; f2 0.5; mu 0.0133/yr | Weber, Weber & Milligan 2001 Math Biosci; MEmilio MSEIRS4 ([docs](https://memilio.readthedocs.io/en/stable/_sources/cpp/models/omseirs4.rst.txt)) | partial V |
| pertussis_hethcote_mild | Pertussis with mild reinfection in the waned | S,Ip,R,V,W,Im | R0_p·(g_rec+mu)·(Ip+rel_m·Im)/N | R0_p 15; g_rec 17.4/yr; rel_m 0.5; w_nat 1/30/yr; w_vac 1/10/yr; p_vac 0.9 | Hethcote 1997 Math Biosci (age removed) | I |
| measles_mseir_2dose | Measles, two-dose MSEIR | M,S0,Sf,V,E,Ip,R | R0_m·(g_rec+mu)·Ip/N | R0_m 15; sig_lat 365/8/yr; g_rec 365/5/yr; xi_m 2/yr; nu1 1.5/yr; nu2 0.3/yr; e1 0.93; e2 0.97 | Hethcote 2000 SIAM Rev; Guerra et al. 2017 R0 review (R) | I |
| tb_mdr_treatment | TB with acquired and transmitted MDR | S,Ls,Lr,Is,Ir,Rt | beta_s·Is/N; beta_s·fit_r·Ir/N | beta_s 10/yr; fit_r 0.8; p_fast 0.14; v_react 5e-4/yr; x_prot 0.65; det 0.7/yr; succ_s 0.85; a_res 0.03; succ_r 0.5; mu 0.02/yr | Dye & Williams 2000 PNAS; Cohen & Murray 2004 Nat Med ([PMC2652755](https://pmc.ncbi.nlm.nih.gov/articles/PMC2652755)) | I (mu, fit_r V via snippet) |
| pneumo_vt_nvt_lipsitch | Pneumococcal serotype replacement | X,Y1,Y2,Y12 + vaccinated Xv,Y1v,Y2v,Y12v | b·(carriers of type)/N | b_vt 0.06, b_nvt 0.055, clr 1/42 (/d); comp 0.3; ve 0.9; p_vac 0.9; mu 1/(5·365)/d | Lipsitch 1997 PNAS ([PMC21091](https://pmc.ncbi.nlm.nih.gov/articles/PMC21091)) | I |
| varicella_zoster_brisson | VZV with exogenous boosting and zoster | M,S,E,Ip,Lb,Zs,Zo,Rz,V | R0_v·(g_rec+mu)·(Ip+rel_z·Zo)/N | R0_v 8; w_boost 1/20/yr (V, Brisson 2002); react 0.004/yr; rel_z 0.2; boost_f 1 (≈0.9 defensible) | Brisson et al. 2000 Epidemiol Infect; 2002 Vaccine ([LSHTM](https://researchonline.lshtm.ac.uk/id/eprint/6784)) | I (+1 V) |
| rubella_crs | Rubella, two age groups, CRS output | Sch,Ich,Rch,Sad,Iad,Rad,V,CRS* | R0_r·g_rec·(Ich+Iad)/N | R0_r 7; g_rec 365/11/yr; age_out 1/15/yr; mu_ad 1/55/yr; p_vac 0.8; preg_frac 0.015; q_crs 0.6 | Anderson & May 1983 J Hyg; Miller et al. 1982 Lancet (R) | I |
| mumps_waning_lewnard | Mumps with waning vaccine immunity | S,E,Ip,R,V | R0_mu·(g_rec+mu)·Ip/N | R0_mu 7; sig_lat 365/17/yr; g_rec 365/8/yr; p_vac 0.92; w_vac 1/27/yr (R) | Lewnard & Grad 2018 Sci Transl Med | R/I |
| meningo_carriage_trotter | MenC carriage, conjugate vaccine (Hib by relabelling) | S,Ca,Rim,V,Cases* | R0_c·(clr+mu)·Ca/N | R0_c 1.36 (V, abstract); clr 2/yr; w_nat 0.5/yr; ve_c 0.75; w_vac 0.1/yr; ccr 1e-4 | Trotter, Gay & Edmunds 2005 AJE ([LSHTM](https://researchonline.lshtm.ac.uk/6806/)) | partial V |
| sars_chowell_seijr | SARS, diagnosis/isolation (MERS template) | S1,S2,E,Ip,J,R,Dead* | beta_s·(Ip+q_e·E+l_iso·J)/N | beta_s 0.25/d; q_e 0.1; l_iso 0.43 HK / 0.49 SG / 0.10 Toronto; p_low 1/3; rho_hr 0.77; k_inc 0.15707; alpha_d 0.2061; g1 0.035285; g2 0.0426; delta_d 0.0279 (/d) | Chowell et al. 2003 JTB; 2004 EID ([CDC](https://wwwnc.cdc.gov/eid/article/10/7/03-0647_article)) | **V** |

```text
flu_arino2008_antiviral: S>L:foi; ST>LT:sig_S*foi; L>Ip:p_sym*k_lat; L>A:(1-p_sym)*k_lat; L>LT:phi_L; LT>L:theta_L; LT>IpT:p_sym*tau_av*k_latT; LT>AT:(1-p_sym*tau_av)*k_latT; Ip>R:f_rec*alpha_r; Ip>Dead:(1-f_rec)*alpha_r; Ip>IpT:phi_I; IpT>Ip:theta_I; IpT>R:fT_rec*alpha_T; IpT>Dead:(1-fT_rec)*alpha_T; A>R:eta_r; A>AT:phi_A; AT>A:theta_A; AT>R:eta_T
   (Lam = eps_L*L + eps_L*sig_L*LT + (1-q_w)*Ip + (1-q_w)*sig_I*IpT + delta_a*A + delta_a*sig_A*AT; eps_L = 0)
flu_two_strain_cross:  as two_strain_cross_immunity, plus R1>S:omega; R2>S:omega; R12>R1:omega/2; R12>R2:omega/2
flu_seasonal_vax:      >S:mu*N[abs]; S>Ip:foi; Ip>R:g_rec; R>S:omega_R; S>V:ve*vrate; V>S:omega_V; S,Ip,R,V>:mu
   (vrate = in_camp*(-log(1-cov)/camp_len); in_camp = (tyear>=camp_start)*(tyear<camp_start+camp_len))
covid_sidarthe:        S>Iu:foi; Iu>Dg:eps_v; Iu>A:zeta; Iu>H:lambda_h; Dg>Rd:eta; Dg>H:rho_h; A>Rd:theta; A>Th:mu_t; A>H:kappa_h; Rd>Th:nu_t; Rd>H:xi_h; Th>Ex:tau_t; Th>H:sigma_h
covid_kucharski2020:   S>E1:foi; E1>E2:2/inc_per; E2>I1:2/inc_per; I1>I2:2/inf_per; I2>Rm:2/inf_per; >Qr:p_rep*(2/inc_per)*E2[abs]; Qr>Cr:1/rep_delay
covid_vax_waning:      S>E:foi; V>E:(1-ve_V)*foi; Wv>E:(1-ve_W)*foi; S>V:vac_rate; V>Wv:wane_V; Wv>V:boost_rate; R>Wv:wane_R; E>Ip:sig_lat; Ip>Hc:p_hosp*g_rec; Ip>R:(1-p_hosp)*g_rec; Hc>R:(1-p_die)*h_rec; Hc>Dead:p_die*h_rec
rsv_weber2001_graded:  >M:mu*N[abs]; M>S1:xi_m; S1>E1:foi1; E1>I1:sig_lat; I1>R1:nu; R1>S2:gam_w; S2>E2:f2*foi1; E2>I2:sig_lat; I2>R2:nu; R2>S2:gam_w; all>:mu
pertussis_hethcote_mild: >S:(1-p_vac)*mu*N[abs]; >V:p_vac*mu*N[abs]; S>Ip:foi; Ip>R:g_rec; R>W:w_nat; V>W:w_vac; W>Im:foi; Im>R:g_rec; all>:mu
measles_mseir_2dose:   >M:mu*N[abs]; M>S0:xi_m; S0>V:e1*nu1; S0>Sf:(1-e1)*nu1; Sf>V:e2*nu2; S0>E:foi; Sf>E:foi; E>Ip:sig_lat; Ip>R:g_rec; all>:mu
tb_mdr_treatment:      >S:mu*N+mu_tb*(Is+Ir)[abs]; S>Is:p_fast*foi_s; S>Ls:(1-p_fast)*foi_s; S>Ir:p_fast*foi_r; S>Lr:(1-p_fast)*foi_r; Ls>Is:v_react+p_fast*(1-x_prot)*foi_s; Ls>Ir:p_fast*(1-x_prot)*foi_r; Lr>Ir:v_react+p_fast*(1-x_prot)*foi_r; Lr>Is:p_fast*(1-x_prot)*foi_s; Rt>Is:relapse+p_fast*(1-x_prot)*foi_s; Rt>Ir:p_fast*(1-x_prot)*foi_r; Is>Rt:cure_nat+det*succ_s; Is>Ir:det*a_res; Ir>Rt:cure_nat+det*succ_r; Is,Ir>:mu+mu_tb; S,Ls,Lr,Rt>:mu
pneumo_vt_nvt_lipsitch: >X:(1-p_vac)*mu*N[abs]; >Xv:p_vac*mu*N[abs]; X>Y1:foi1; X>Y2:foi2; Y1>Y12:comp*foi2; Y2>Y12:comp*foi1; Y1,Y2>X:clr; Y12>Y2:clr; Y12>Y1:clr; Xv>Y1v:(1-ve)*foi1; Xv>Y2v:foi2; Y1v>Y12v:comp*foi2; Y2v>Y12v:comp*(1-ve)*foi1; Y1v,Y2v>Xv:clr; Y12v>Y2v:clr; Y12v>Y1v:clr; all>:mu
varicella_zoster_brisson: >M:(1-p_vac)*mu*N[abs]; >V:p_vac*mu*N[abs]; M>S:xi_m; S>E:foi; E>Ip:sig_lat; Ip>Lb:g_rec; Lb>Zs:w_boost; Zs>Lb:boost_f*foi; Zs>Zo:react; Zo>Rz:z_rec; all>:mu
rubella_crs:           >Sch:(1-p_vac)*births[abs]; >V:p_vac*births[abs]; Sch>Ich:foi; Ich>Rch:g_rec; Sad>Iad:foi; Iad>Rad:g_rec; Sch>Sad,Ich>Iad,Rch>Rad:age_out; Sad,Iad,Rad>:mu_ad; V>:1/70; >CRS:q_crs*preg_frac*foi*Sad[abs]
mumps_waning_lewnard:  >S:(1-p_vac)*mu*N[abs]; >V:p_vac*mu*N[abs]; S>E:foi; E>Ip:sig_lat; Ip>R:g_rec; V>S:w_vac; S>V:boost_v; all>:mu
meningo_carriage_trotter: >S:(1-p_vac)*mu*N[abs]; >V:p_vac*mu*N[abs]; S>Ca:foi; V>Ca:(1-ve_c)*foi; Ca>Rim:clr; Rim>S:w_nat; V>S:w_vac; S>V:vac_rate; Rim>V:vac_rate; >Cases:ccr*acq[abs]; S,Ca,Rim,V>:mu
sars_chowell_seijr:    S1>E:foi; S2>E:p_low*foi; E>Ip:k_inc; Ip>J:alpha_d; Ip>R:g1; Ip>Dead:delta_d; J>R:g2; J>Dead:delta_d
```

Some secondary sources support these models.

- **Lipsitch 1997.** With two serotypes, the rise in non-vaccine-type carriage cannot exceed the fall in vaccine-type carriage ([PMC21091](https://pmc.ncbi.nlm.nih.gov/articles/PMC21091)). Whether the paper allows dual carriage still needs checking; setting `comp = 0` gives complete exclusion.
- **Cohen & Murray 2004.** The unfit MDR strain has relative fitness 0.3 and the fit strain 0.8–1.2. Mutation to the unfit phenotype is 100 times as frequent as to the fit one ([PMC2652755](https://pmc.ncbi.nlm.nih.gov/articles/PMC2652755)).
- **Brisson 2002.** Exposure to varicella boosts zoster immunity for about 20 years (95% CI 7–41) ([LSHTM](https://researchonline.lshtm.ac.uk/id/eprint/6784)).
- **Trotter 2005.** Results were sensitive to the duration of vaccine protection and to vaccine efficacy against carriage ([LSHTM](https://researchonline.lshtm.ac.uk/6806/)).

## 4. Vector-borne diseases: 14 specifications, with Chitnis 2008 and Bowman 2005 fully verified

Two vector models were read in full.

- **Chitnis, Hyman & Cushing 2008.** Equations 1a–g and the Table 3 baseline sets were read from the author preprint. The high-transmission set gives R0 = 4.4 and the low-transmission set R0 = 1.1. The biting term moves smoothly from frequency dependence to mass action ([Chitnis et al. 2008 preprint](https://archive.math.arizona.edu/cushing/DeterminingImportantParametersInTheSpreadofMalaria%20BMB2008.pdf)). This entry also covers the Ngwa & Shu SEIRS–SEI lineage.
- **Bowman et al. 2005.** Table 1 was read from the YorkSpace copy. Mosquito recruitment Π_M and death μ_M are listed as "Variable", so those two values are illustrative ([Bowman et al. 2005](https://yorkspace.library.yorku.ca/bitstream/handle/10315/908/2005_bmb_wnv.pdf;sequence=1)).

Two more have verified parameter tables but equations rebuilt from text, because the equations appear only as images.

- **Gao et al. 2016 (Zika).** R0 = 2.055, with sexual transmission contributing about 3.04% ([PMC4911567](https://pmc.ncbi.nlm.nih.gov/articles/PMC4911567)).
- **Yakob & Clements 2013 (chikungunya).** β1 = 0.14 and β2 = 0.40 were fitted to Réunion data ([PLoS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0057448)).

| ID | Disease / area | Compartments (vectors not in N) | Force of infection | Typical parameters (units) | Primary source | Prov. |
|---|---|---|---|---|---|---|
| malaria_chitnis2008_seirs_sei | P. falciparum malaria | Sh,Eh,Ih,Rh; Sv,Ev,Iv | bites = sig_v·Nv·sig_h·Nh/(sig_v·Nv+sig_h·Nh); lambda_h = (bites/Nh)·bet_hv·Iv/Nv; lambda_v = (bites/Nv)·(bet_vh·Ih+bet_vh_r·Rh)/Nh | High [low]: Lambda_h 0.033 [0.041]/d; psi_h 1.1e-4 [5.5e-5]; psi_v 0.13; sig_v 0.5 [0.33]; sig_h 19 [4.3]; bet_hv 0.022; bet_vh 0.48 [0.24]; bet_vh_r 0.048 [0.024]; nu_h 0.1; nu_v 0.091 [0.083]; gam_h 0.0035; delta_h 9e-5 [1.8e-5]; rho_h 5.5e-4 [2.7e-3]; mu1h 1.6e-5 [8.8e-6]; mu2h 3e-7 [2e-7]; mu1v 0.033; mu2v 2e-5 [4e-5] | Chitnis, Hyman & Cushing 2008 Bull Math Biol ([preprint](https://archive.math.arizona.edu/cushing/DeterminingImportantParametersInTheSpreadofMalaria%20BMB2008.pdf)) | **V** |
| malaria_superinfection_boosting | Stable endemic malaria | Sh,Ih,Rh; Sv,Ev,Iv | h = a·p_mh·Iv/Nh; r_eff = h/(exp(h/r)−1); w_loss = h/(exp(h·tau_imm)−1) | a_bite 0.3/d; p_mh 0.1; p_hm 0.5; r_clear 0.005/d; tau_imm 1825 d; mu_v 0.1/d; nu_v 0.1/d; m 5 | Macdonald/Dietz superinfection; Aron 1988; Smith et al. 2012 PLoS Pathog (R) | I |
| dengue_two_serotype | DENV, two serotypes | Sh,I1,I2,R1,R2,I12,I21,Rh; Sv,Iv1,Iv2 | lam_k = a·p_vh·Iv_k/Nh; lamv_k = a·p_hv·(I_k+I_secondary,k)/Nh | a 0.5/d; p_vh=p_hv 0.5; gam_h 1/7/d; mu_h 4.2e-5/d; mu_v 1/14/d; A_v 1428/d; x_sec 1 (>1 = ADE) | Esteva & Vargas 2003 J Math Biol (R) | I |
| zika_gao2016_sexual | Zika, Latin America 2015–16 | Sh,Eh,I1,I2,Ah,Rh; Sv,Ev,Iv | lambda_h = a·b·Iv/Nh + bet_sex·(kappa·Eh+I1+tau_c·I2)/Nh; lambda_v = a·c·(eta·Eh+I1)/Nh | a 0.5; b 0.4; c 0.5; eta 0.1; bet_sex 0.05/d; kappa 0.6; tau_c 0.3; theta 0.18; m 5; 1/nu_h 5 d; 1/nu_v 10 d; 1/gam_h1 5 d; 1/gam_h2 20 d; 1/gam_ha 7 d; 1/mu_v 14 d | Gao et al. 2016 Sci Rep ([PMC4911567](https://pmc.ncbi.nlm.nih.gov/articles/PMC4911567)) | **V** (params) |
| chikungunya_yakob2013_reunion | Chikungunya, Réunion 2005–06 | Sh,Eh,Ish,Iah,Rh; Xv,Yv,Zv | beta1·Zv/Nv; beta2·(Ish+Iah)/Nh | beta1 0.14/d; beta2 0.40/d; lat_h=lat_v 0.5/d; omega 0.25/d; phi_s 0.97; gam_h 0.25/d; mu_v 0.05/d; R0≈4.1 (calculated) | Yakob & Clements 2013 PLoS ONE ([link](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0057448)) | **V** (params) |
| wnv_bowman2005 | West Nile virus, North America | Mu,Mi,Bu,Bi (not in N); Sh,Eh,Ih,Hh,Rh | lam_M = b1·bet1·Bi/NB; lam_B = b1·bet2·Mi/NB; lam_H = b2·bet3·Mi/NH | Pi_B 1000/d; mu_B 0.001/d; Pi_H 30/d; mu_H 3.914e-5/d; b1=b2 0.09/d; bet1 0.16; bet2=bet3 0.88; d_B 5e-5; d_H 5e-7; alpha_h 1/14; delta_h 1; tau_h 1/14; Pi_M 20000, mu_M 0.05 (I) | Bowman et al. 2005 Bull Math Biol ([YorkSpace](https://yorkspace.library.yorku.ca/bitstream/handle/10315/908/2005_bmb_wnv.pdf;sequence=1)) | **V** |
| lyme_tick_mouse | Lyme borreliosis | L,Ns,Ni,As,Ai (ticks); Ms,Mi | Ms→Mi at a_N·p_NM·Ni; larval infection at a_L·s_L·p_ML·Mi | a_L=a_N 1e-4/(mouse·d); a_A 0.005; Deer 5; s_L=s_N 0.3; eps_egg 1000; p_NM 0.8; p_ML 0.5; mu_m 0.005/d; K_m 100 | Rosà et al. 2003 JTB; Rosà & Pugliese 2007 (R) | I |
| vl_leishmaniasis_isc | Visceral leishmaniasis, India/Nepal | Sh,Ah,Kh,Tx,Ph,Rh; Sf,Ef,If | a_sf·p_fh·If/Nh; a_sf·(k_A·Ah+k_K·Kh+k_Tx·Tx+k_P·Ph)/Nh | a_sf 0.25/d; p_fh 0.1; k_A 0.025; k_K 0.5; k_Tx 0.25; k_P 0.3; w_A 1/150/d; p_K 0.1; f_P 0.025; mu_f 1/14/d; nu_f 0.2/d | Stauch et al. 2011 PLoS NTD ([DOAJ](https://doaj.org/article/090474c305d14de69ead51f19b49887d)) | I |
| chagas_household_cohen2001 | Chagas, rural Argentina | Hs,Hi,Gs,Gi; Vs,Vi | bite share f_x = w_x·N_x/wsum (including chickens); lam_hum = a_b·f_hum·p_vh·Vi/Nh | a_b 1/7/d; w_dog 2.5; N_chk 5; p_vh 5.8e-4 (R, Nouvellet 2013); p_gv 0.49; p_hv 0.05; mu_dog 9.13e-4/d; K_bug 1000 | Cohen & Gürtler 2001 Science (R) | I |
| yellow_fever_urban_vacc | Urban yellow fever | Sh,Eh,Ih,Sev,Rh,Vh; Sv,Ev,Iv | a·p_vh·Iv/Nh; a·p_hv·Ih/Nh | nu_h 1/4.3/d; nu_v 1/10/d (Johansson 2010, R); p_sev 0.12; cfr_sev 0.47 (Johansson 2014, R); v_eff 0.99 (WHO, R); v_rate 0.01/d from day 60 | Johansson et al. 2010, 2014; WHO 2013 (R) | R/I |
| rvf_gaff2007 | Rift Valley fever | Pa,Qa,Sa,Ea,Ia,Pc,Sc,Ec,Ic (vectors); Sl,El,Il,Rl | (a_a·p_al·Ia + a_c·p_cl·Ic)/Nl | q_tov 0.05; th_a 0.05/d; mu_ea 0.001/d; eps 1/7/d; a 0.3/d; p 0.3; gam_l 0.25/d; d_l 0.025/d; flood window 60 d (all I) | Gaff, Hartley & Leahy 2007 EJDE ([abstract](https://ejde.math.txstate.edu/Volumes/2007/115/abstr.html)) | I |
| plague_keeling_gilligan2000 | Bubonic plague | Srat,Irat,Rrat,Nf,Ff (not in N); Shum,Ihum,Rhum | lam_rat = bet_R·Ff·(1−exp(−a_s·Nrat))/Nrat; human: bet_H·Ff·exp(−a_s·Nrat)/Nhum | r_R 5/yr; p_inh 0.975; K_R 2500; d_R 0.2/yr; bet_R 4.7/yr; a_s 0.004; m_R 20/yr; g_R 0.02; r_F 20/yr; K_F 6.57; d_F 10/yr | Keeling & Gilligan 2000 Nature ([link](https://www.nature.com/articles/35038073)); Proc R Soc B ([PMC1690796](https://pmc.ncbi.nlm.nih.gov/articles/PMC1690796)) | R (human block I) |
| trypanosomiasis_rogers1988 | Rhodesian HAT / animal trypanosomiasis | Sh,Eh,Ih,Rh; Sc,Ec,Ic,Rc; Sf,Ef1,Ef2,If | lam_h = a_f·f_hum·b_fh·If/Nh; lam_f summed over hosts | a_f 0.25/d; f_hum 0.1; b 0.6; p 0.05; mu_f 0.03/d; eip 25 d (2 Erlang stages) | Rogers 1988 Parasitology (R) | I |
| bluetongue_gubbins2008 | BTV-8, northern Europe | Scow,Icow,Rcow,Ssh,Ish,Rsh; Sv,Ev1–Ev3,Iv | temperature-driven a(T), mu(T), nu(T); host-preference weighting phi_c | a(T) = 0.0002·T·(T−3.7)·(41.9−T)^(1/2.7); mu(T) = 0.009·e^(0.16T); nu(T) = 0.018·(T−13.4); b_vh 0.9; p_hv 0.05; r_cow 1/20.6/d; r_sh 1/16/d | Gubbins et al. 2008 J R Soc Interface (R) | R |

```text
malaria_chitnis2008:  >Sh:Lambda_h[abs]; >Sh:psi_h*Nh[abs]; Sh>Eh:lambda_h; Eh>Ih:nu_h; Ih>Rh:gam_h; Rh>Sh:rho_h; Ih>:delta_h; Sh,Eh,Ih,Rh>:fh(=mu1h+mu2h*Nh); >Sv:psi_v*Nv[abs]; Sv>Ev:lambda_v; Ev>Iv:nu_v; Sv,Ev,Iv>:fv(=mu1v+mu2v*Nv)
malaria_superinfection_boosting: >Sh:mu_h*Nh[abs]; Sh>Ih:h_foi; Ih>Rh:r_eff; Rh>Sh:w_loss; Sh,Ih,Rh>:mu_h; >Sv:mu_v*Nv[abs]; Sv>Ev:lambda_v; Ev>Iv:nu_v; Sv,Ev,Iv>:mu_v   (h_foi floored at 1e-10)
dengue_two_serotype:  >Sh:mu_h*Nh[abs]; Sh>I1:lam1; Sh>I2:lam2; I1>R1,I2>R2:gam_h; R1>I12:x_sec*lam2; R2>I21:x_sec*lam1; I12,I21>Rh:gam_h; humans>:mu_h; >Sv:A_v[abs]; Sv>Iv1:lamv1; Sv>Iv2:lamv2; Sv,Iv1,Iv2>:mu_v
zika_gao2016:         Sh>Eh:lambda_h; Eh>I1:theta*nu_h; Eh>Ah:(1-theta)*nu_h; I1>I2:gam_h1; I2>Rh:gam_h2; Ah>Rh:gam_ha; >Sv:mu_v*Nv[abs]; Sv>Ev:lambda_v; Ev>Iv:nu_v; Sv,Ev,Iv>:mu_v
chikungunya_yakob2013: Sh>Eh:lambda_h; Eh>Ish:phi_s*lat_h; Eh>Iah:(1-phi_s)*lat_h; Ish,Iah>Rh:gam_h; >Xv:mu_v*Nv[abs]; Xv>Yv:lambda_v; Yv>Zv:lat_v; Xv,Yv,Zv>:mu_v
wnv_bowman2005:       >Mu:Pi_M[abs]; Mu>Mi:lam_M; Mu,Mi>:mu_M; >Bu:Pi_B[abs]; Bu>Bi:lam_B; Bu,Bi>:mu_B; Bi>:d_B; >Sh:Pi_H[abs]; Sh>Eh:lam_H; Eh>Ih:alpha_h; Ih>Hh:delta_h; Hh>Rh:tau_h; Hh>:d_H; humans>:mu_H
lyme_tick_mouse:      >L:eps_egg*a_A*Deer*(As+Ai)[abs]; L>Ns:a_L*s_L*(Hm-p_ML*Mi); L>Ni:a_L*s_L*p_ML*Mi; L>:a_L*Hm*(1-s_L)+mu_L+mu_L2*L; Ns>As:a_N*s_N*(Hm-p_ML*Mi); Ns>Ai:a_N*s_N*p_ML*Mi; Ns>:a_N*Hm*(1-s_N)+mu_N; Ni>Ai:a_N*s_N*Hm; Ni>:a_N*Hm*(1-s_N)+mu_N; As,Ai>:a_A*Deer+mu_A; >Ms:mu_m*K_m[abs]; Ms>Mi:a_N*p_NM*Ni; Ms,Mi>:mu_m
vl_leishmaniasis_isc: >Sh:mu_h*Nh[abs]; Sh>Ah:lambda_h; Ah>Kh:p_K*w_A; Ah>Rh:(1-p_K)*w_A; Kh>Tx:r_tx; Kh>:d_K; Tx>Ph:f_P*w_Tx; Tx>Rh:(1-f_P)*w_Tx; Ph>Rh:w_P; Rh>Sh:rho_h; humans>:mu_h; >Sf:mu_f*Nf[abs]; Sf>Ef:lambda_f; Ef>If:nu_f; Sf,Ef,If>:mu_f
chagas_household:     >Hs:mu_hum*Nh[abs]; Hs>Hi:lam_hum; Hs,Hi>:mu_hum; >Gs:mu_dog*Ng[abs]; Gs>Gi:lam_dog; Gs,Gi>:mu_dog; >Vs:r_bug*Nv[abs]; Vs>Vi:lam_bug; Vs,Vi>:mu_bug+(r_bug-mu_bug)*Nv/K_bug
yellow_fever_urban:   Sh>Eh:lambda_h; Sh>Vh:vacc; Eh>Ih:nu_h; Ih>Sev:p_sev*gam_h; Ih>Rh:(1-p_sev)*gam_h; Sev>:cfr_sev*w_sev; Sev>Rh:(1-cfr_sev)*w_sev; >Sv:mu_v*Nv[abs]; Sv>Ev:lambda_v; Ev>Iv:nu_v; Sv,Ev,Iv>:mu_v
rvf_gaff2007:         >Pa:b_a*(Sa+Ea)+b_a*(1-q_tov)*Ia[abs]; >Qa:b_a*q_tov*Ia[abs]; Pa>Sa:th_a*flood*pmax(0,1-Nas/K_a); Qa>Ia:(same); Pa,Qa>:mu_ea; Sa>Ea:lam_a; Ea>Ia:eps_a; Sa,Ea,Ia>:mu_a; >Pc:b_c*Ncs[abs]; Pc>Sc:th_c*flood*pmax(0,1-Ncs/K_c); Pc>:mu_ec; Sc>Ec:lam_c; Ec>Ic:eps_c; Sc,Ec,Ic>:mu_c; >Sl:mu_l*Nl[abs]; Sl>El:lam_l; El>Il:eps_l; Il>Rl:gam_l; Il>:d_l; livestock>:mu_l
plague_keeling_gilligan: >Srat:r_R*Srat[abs]; >Srat:r_R*(1-p_inh)*Rrat[abs]; >Rrat:r_R*p_inh*Rrat[abs]; Srat,Rrat>:d_R+r_R*Nrat/K_R; Srat>Irat:lam_rat; Irat>Rrat:m_R*g_R; Irat>:m_R*(1-g_R)+d_R; >Nf:r_F*Nf[abs]; Nf>:r_F*Nf/K_F; >Ff:(d_R+m_R*(1-g_R))*Irat*Nf[abs]; Ff>:d_F; >Shum:mu_hum*Nhum[abs]; Shum>Ihum:lam_hum; Ihum>Rhum:m_H*g_H; Ihum>:m_H*(1-g_H); humans>:mu_hum
trypanosomiasis_rogers: hosts SEIRS with births = deaths; Ih>:d_h; >Sf:mu_f*Nf[abs]; Sf>Ef1:lam_f; Ef1>Ef2:2/eip; Ef2>If:2/eip; flies>:mu_f
bluetongue_gubbins2008: Scow>Icow:lam_cow; Icow>Rcow:r_cow; Ssh>Ish:lam_sh; Ish>Rsh:r_sh; Ish>:d_sh; >Sv:mu_v*Nv[abs]; Sv>Ev1:lam_v; Ev1>Ev2>Ev3>Iv:3*nu_eip each; midges>:mu_v
```

Of the reconstructed models, Stauch et al. found that asymptomatic infections maintain visceral leishmaniasis (VL) transmission and that vector control beats case treatment ([DOAJ](https://doaj.org/article/090474c305d14de69ead51f19b49887d)). The 2014 follow-up puts the threshold for elimination at a 67% cut in sandfly density ([DOAJ](https://doaj.org/article/c75fd498252a425ab32ca195c75984f8)). The Gaff 2007 rift valley fever (RVF) values and the Keeling & Gilligan plague values are open access but were not transcribed this session. The EJDE fetch failed and the plague values are recalled. Both should be pulled before release. The Lyme, VL, Chagas, trypanosomiasis and RVF entries are composites, and their numbers are demonstration values, not calibrations.

## 5. Sexually transmitted and blood-borne infections: 12 specifications, with Whittles 2017 and Martin 2013 as anchors

Two models have parameter tables read from the paper.

- **Whittles, White & Didelot 2017 (gonorrhoea AMR).** Posterior means are θ = 5.2/yr, ψ = 0.63, σ = 77/yr, ν = 1.8/yr, a resistant-strain clearance multiplier of 1.8, μ_t = 136/yr, ρ = 53/yr and cefixime failure of 0.83, in a closed population of 0.6 million MSM in England ([PMC5663337](https://pmc.ncbi.nlm.nih.gov/articles/PMC5663337)).
- **Martin et al. 2013 (HCV in PWID).** Full Edinburgh, Melbourne and Vancouver settings are available. Halving prevalence in 15 years takes 15, 40 or 76 treatments per 1000 PWID per year respectively ([PMC3933734](https://pmc.ncbi.nlm.nih.gov/articles/PMC3933734)).

For both, the ODEs sit in supplements that were not readable. The flows are reconstructed from the text, and which Whittles compartments are infectious is an assumption. The syphilis stage durations come from Garnett 1997 via Feldman & Mishra 2019 ([PMC6806446](https://pmc.ncbi.nlm.nih.gov/articles/PMC6806446)): incubation 25 days, infectious period 154 days and protective immunity 5 years. The chlamydia asymptomatic duration of 433 days (95% CI 420–447) is from Althaus et al. 2010 ([BORIS](https://boris.unibe.ch/1316)).

| ID | Disease / area | Compartments | Force of infection | Typical parameters (units) | Primary source | Prov. |
|---|---|---|---|---|---|---|
| hiv_art_tasp_granich | HIV, 4 untreated + 4 ART stages | S,I1–I4,A1–A4 | lam_eff·(ΣI + eps_art·ΣA)/N, with lam_eff = lambda0·exp(−alpha_h·P^n_h) | lambda0 0.6/yr; alpha_h 4; rho_p 0.4/yr; rho_a 0.1/yr; eps_art 0.01 (R); tau 1/yr from t_start 10; dropout 0.015/yr; mu 0.02/yr | Granich et al. 2009 Lancet ([LSHTM](https://researchonline.lshtm.ac.uk/id/eprint/312973); [PDF](https://mysite.science.uottawa.ca/rsmith43/mat4996/granichgilks.pdf)) | I/R |
| hiv_core_noncore_mixing | HIV, core/non-core with epsilon mixing | S_h,I_h,Y_h,A_h,S_l,I_l,Y_l,A_l | cp_x·beta_p·(eps_mix·own-group prevalence + (1−eps_mix)·prop_inf) | cp_h 10, cp_l 1 partners/yr; beta_p 0.05; eps_mix 0.5; nu_i 0.1/yr; delta_y 1/yr; tau 0.2/yr; eps_art 0.04 (HPTN 052, R); f_h 0.05 | Hethcote & Yorke 1984; Garnett & Anderson 1993 (R) | I |
| hiv_prep_msm | HIV among MSM, PrEP + treatment as prevention | S,Sp,Ia,Ic,Dx,Av | cp·beta_c·(ra·Ia+Ic+Dx+eps_art·Av)/N | cp 5/yr; beta_c 0.1; ra 10; sig_a 4/yr; eff_prep 0.86 (PROUD, R); kappa_p 0.2/yr; omega_p 0.5/yr; tau_d 0.5/yr; tau_a 2/yr | Composite; McCormack et al. 2016 Lancet (R) | I |
| gonorrhoea_core_group_hethcote_yorke | Gonorrhoea, two sexes × core/non-core SIS | Sm_c,Im_c,Sm_n,Im_n,Sf_c,If_c,Sf_n,If_n | epsilon mixing across sexes | cp_c 20, cp_n 1 /yr; b_mf 0.5; b_fm 0.25; dur_m 0.05 yr; dur_f 0.25 yr; eps_mix 0.6 | Hethcote & Yorke 1984 (R) | I |
| gonorrhoea_amr_whittles | Gonorrhoea, cefixime resistance, MSM England | S,Us,Es,As,Tsc,Tso,Ur,Er,Ar,Trc,Tro | theta·Inf_strain/N | theta 5.2/yr; psi 0.63; sig_inc 77/yr; nu 1.8/yr; alpha_r 1.8; mu_t 136/yr; rho 53/yr; phi_f 0.83; pi_cef 0.3 (I) | Whittles et al. 2017 PLoS Med ([PMC5663337](https://pmc.ncbi.nlm.nih.gov/articles/PMC5663337)) | **V** |
| syphilis_staged_garnett | Syphilis, staged with waning immunity | S,E,P1,P2,Le,Ll,R | cp·beta_s·(P1+P2)/N | cp 15/yr (V); beta_s 0.3 (I; per-act 0.09–0.64 V); sig_inc 365/25/yr (V); d_inf 154/365 yr (V); tau 0.4/yr (range 0.1–0.8 V); w_imm 1/5/yr (V); mu 5831/366279/yr (V); rel 0.25/yr (I) | Garnett et al. 1997 STD; Feldman & Mishra 2019 ([PMC6806446](https://pmc.ncbi.nlm.nih.gov/articles/PMC6806446)) | partial V |
| chlamydia_seirs_screening_althaus | Chlamydia, two sexes, with screening | Sf,Ef,Iaf,Isf,Rf,Sm,Em,Iam,Ism,Rm | cp·beta·(infected of other sex)/N_other | d_af 433 d (V); cp 1/365/d; beta 0.3; sig_inc 1/14/d; p_sym_f 0.1; p_sym_m 0.2; d_s 30 d; w_imm 1/365/d; scr_f 0.2/365/d | Althaus et al. 2010 Epidemics ([BORIS](https://boris.unibe.ch/1316)) | partial V |
| hpv_vaccination_simplified | HPV 16/18, two sexes, vaccination | Sf,If,Rf,Vf,Sm,Im,Rm,Vm | cp·beta·I_other/N_other | cp 2/yr; beta_mf 0.8; beta_fm 0.7; clr_f 1/yr; clr_m 1.5/yr; p_imm_f 0.5; w_nat 0.1/yr; vc_f 0.7; ve 0.9; mu 1/30/yr | Elbasha, Dasbach & Insinga 2007 EID ([CDC](https://www.cdc.gov/eid/article/13/1/pdfs/06-0438-combined.pdf)) | I |
| hsv2_latent_reactivation | HSV-2 with latency and recurrence | S,Ip,Lt,Rc | cp·(beta_p·Ip+beta_r·Rc+beta_l·Lt)/N | cp 2/yr; beta_p 0.5; beta_r 0.3; beta_l 0.05; d_p 3/52 yr; react 4/yr; d_r 5/365 yr | Blower, Porco & Darby 1998 Nat Med (R) | I |
| hbv_carrier_vaccination | Hepatitis B with carriers, perinatal transmission, vaccination | S,Ac,Cr,R,V | beta_h·(Ac+eps_c·Cr)/N | beta_h 0.95/yr; eps_c 0.16; gam_a 4/yr; q_c 0.0886; gam_c 0.025/yr; pv 0.11; pc0 0.885; mu 0.00693/yr; w_v 0.1/yr (all R, Zou et al. 2010) | Medley et al. 2001 Nat Med; Zou, Zhang & Ruan 2010 JTB; Edmunds 1996 ([PMC2271713](https://pmc.ncbi.nlm.nih.gov/articles/PMC2271713/)) | R |
| hcv_pwid_treatment_martin | HCV in PWID, DAA treatment as prevention | X0,A0,Ch0,Tx0,Fail0,X1,A1,Ch1,Tx1,Fail1 | foi0 = pi_inf·Infectious/Ntot; foi1 = Pi_h·foi0 | Pi_h 3.6; kappa 12/14/yr; delta_c 0.26; psi_a 2/yr; mu1 1/11/yr; mu2 0.01/yr; phi_tx 8 per 1000/yr; omega 52/12/yr; alpha_svr 0.9 (all V); pi_inf 0.2/yr (calibrate) | Martin et al. 2013 Hepatology ([PMC3933734](https://pmc.ncbi.nlm.nih.gov/articles/PMC3933734)) | **V** |
| htlv1_sexual_vertical | HTLV-1, sexual + mother-to-child | Sf,Hf,Sm,Hm,Dz | cp·b·H_other/N_other | b_mf 0.06; b_fm 0.01; p_mtct 0.2; dis 0.001/yr; mu 1/50/yr; mu_d 0.5/yr | Composite; Gessain & Cassar 2012 (R) | I |

```text
hiv_art_tasp_granich:  >S:mu*N+rho_p*I4+rho_a*A4[abs]; S>I1:foi; I1>I2>I3>I4:rho_p; I4>:rho_p; Ik>Ak:tau_t; A1>A2>A3>A4:rho_a; A4>:rho_a; Ak>Ik:dropout; all>:mu
hiv_core_noncore_mixing: >S_h:f_h*(mu*N+delta_y*(Y_h+Y_l))[abs]; >S_l:(1-f_h)*(...)[abs]; S_g>I_g:foi_g; I_g>Y_g:nu_i; I_g>A_g:tau; Y_g>:delta_y; all>:mu
hiv_prep_msm:          >S:mu*N+nu_d*(Ic+Dx)[abs]; S>Sp:kappa_p; Sp>S:omega_p; S>Ia:foi; Sp>Ia:(1-eff_prep)*foi; Ia>Ic:sig_a; Ic>Dx:tau_d; Dx>Av:tau_a; Ic,Dx>:nu_d; all>:mu
gonorrhoea_core_group_hy: S*_g>I*_g:foi_*_g; Im_g>Sm_g:1/dur_m; If_g>Sf_g:1/dur_f
gonorrhoea_amr_whittles: S>Us:foi_s; S>Ur:foi_r; Ux>Ex:psi*sig_inc; Ux>Ax:(1-psi)*sig_inc; As>S:nu; Ar>S:alpha_r*nu; Ex>Txc:pi_cef*mu_t; Ex>Txo:(1-pi_cef)*mu_t; Tsc,Tso,Tro>S:rho; Trc>S:(1-phi_f)*rho; Trc>Ar:phi_f*rho
syphilis_staged_garnett: >S:mu*N[abs]; S>E:foi; E>P1:sig_inc; P1>P2:1/(f_prim*d_inf); P2>Le:1/((1-f_prim)*d_inf); Le>P2:rel; Le>Ll:el; P1,P2,Le>R:tau; Ll>R:tau_l; R>S:w_imm; all>:mu
chlamydia_althaus:     Sx>Ex:foi_x; Ex>Isx:p_sym_x*sig_inc; Ex>Iax:(1-p_sym_x)*sig_inc; Iax>Rx:1/d_ax; Iax>Sx:scr_x; Isx>Sx:1/d_s; Rx>Sx:w_imm
hpv_vaccination:       >Sx:mu*Nx*(1-vc_x)[abs]; >Vx:mu*Nx*vc_x[abs]; Sx>Ix:foi_x; Vx>Ix:(1-ve)*foi_x; Ix>Rx:p_imm_x*clr_x; Ix>Sx:(1-p_imm_x)*clr_x; Rx>Sx:w_nat; Vx>Sx:w_vac; all>:mu
hsv2_latent_reactivation: >S:mu*N[abs]; S>Ip:foi; Ip>Lt:1/d_p; Lt>Rc:react*(1-tx); Rc>Lt:1/d_r; all>:mu
hbv_carrier_vaccination: >S:births_ok*(1-vax)+mu_c*Cr[abs]; >V:births_ok*vax[abs]; >Cr:mu*pv*Cr*pc0[abs]; >Ac:mu*pv*Cr*(1-pc0)[abs]; S>Ac:foi; Ac>Cr:q_c*gam_a; Ac>R:(1-q_c)*gam_a; Cr>R:gam_c; V>S:w_v; S,Ac,R,V>:mu; Cr>:mu+mu_c
hcv_pwid_martin:       >X0:(1-f_high)*exit_r*Ntot[abs]; >X1:f_high*exit_r*Ntot[abs]; Xk>Ak:fook; Ak>Xk:delta_c*psi_a; Ak>Chk:(1-delta_c)*psi_a; Chk>Txk:tx_pc; Txk>Xk:alpha_svr*omega; Txk>Failk:(1-alpha_svr)*omega; (X1,A1,Ch1,Tx1,Fail1)>(...0):kappa; all>:exit_r
   (tx_pc = pmin(phi_tx/1000*Ntot/max(Chtot,1e-6), 52))
htlv1_sexual_vertical: >Sf,Sm:0.5*births*(1-pinf_birth)[abs]; >Hf,Hm:0.5*births*pinf_birth[abs]; Sf>Hf:foi_f; Sm>Hm:foi_m; Hf,Hm>Dz:dis; Dz>:mu_d+mu; Sf,Hf,Sm,Hm>:mu
```

Format limits show up most clearly in this group. Without age, the HBV model cannot reproduce Medley's catastrophic bistability, in which starting carrier prevalence is the most influential predictor of intervention outcomes ([record](https://acnpsearch.tweb-dev.unibo.it/singlejournalindex/82970255)). The HPV model drops cohort vaccination, CIN stages and cancer stages. Whittles's original is stochastic and Turner 2006's chlamydia model is individual-based, so the first is approximated by its mean-field analogue and the second was skipped. Fitted transmission coefficients (`lambda0` for Granich, `pi_inf` for Martin) must be calibrated by the user. The Granich exp(−αPⁿ) saturation form is itself recalled and needs checking against the PDF.

## 6. Water-borne, enteric, environmental and neglected tropical diseases: 13 specifications, with Codeço 2001 checked by back-calculation

The Codeço cholera entry is the best-verified model in this group. Table 2 for "Community 2" was read: n = 1e-4/d, a = 1/d, K = 10⁶ cells/ml, r = 0.2/d, nb − mb = −0.33/d, e = 10 cells/ml/d/person. The **B/(K+B) dose–response is confirmed by calculation**: S_crit = K·r·(mb−nb)/(e·a) gives the paper's published 6600 and 13,200 ([Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC29087/fullTextXML)). Two other papers supplied partial values. **Pitzer et al. 2014** give typhoid R0 = 2.8 with 90% from short-cycle transmission, carrier relative infectiousness 0.01, 10% of infections clinical, a reporting fraction of 0.0052 and age-specific carriage of 0.3–10.1% ([PLoS NTD](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0002642)). **Simmons et al. 2013** give norovirus disease immunity of 4.1–8.7 years and R0 = 1.64–1.88 for most model variants ([CDC EID](https://wwwnc.cdc.gov/eid/article/19/8/13-0472_article)). The Simmons paper is internally inconsistent: its abstract gives the Model C immunity CI as 6.8–11.3 years, while the Results give 4.0–11.3.

| ID | Disease / area | Compartments | Force of infection | Typical parameters (units) | Primary source | Prov. |
|---|---|---|---|---|---|---|
| cholera_codeco2001 | Cholera with aquatic reservoir | S,I,R; B (not in N) | a·B/(K+B) | mu_h 1e-4/d; a 1/d (Community 1: 0.5); K 1e6 cells/ml; r 0.2/d; nb_mb −0.33/d; e 10 | Codeço 2001 BMC Infect Dis ([Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC29087/fullTextXML)) | **V** |
| cholera_andrews_basu2011 | Haiti cholera, OCV + water intervention | S,I,A,R,V; W | beta_w·water_mult·W/(kap_w+W) + beta_h·(I+rel_a·A)/N | beta_w 1/d; beta_h 0.1/d; kap_w 1e6; p_sym 0.2; gam_c 0.2/d; del_w 0.033/d; red_w −log(0.99)/7/d (V scenario) | Andrews & Basu 2011 Lancet ([PMC3172163](https://pmc.ncbi.nlm.nih.gov/articles/PMC3172163)) | I/R |
| typhoid_pitzer2014 | Typhoid, carriers, short + long cycle | S1,I1,Car,R,S2,I2; W | beta_s·prev_w + beta_l·(1+amp_w·cos(2πt/52))·W | R0_ty 2.8; frac_short 0.9; r_car 0.01; p_clin 0.1; rep_f 0.0052 (V); theta_c 0.03 (I average of 0.3–10.1% V range); delta_i 0.25/wk; omega_r 0.0096/wk; xi_w 0.333/wk (R) | Pitzer et al. 2014 PLoS NTD ([link](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0002642)) | partial V |
| polio_wpv_opv | WPV, Sabin OPV, cVDPV | S,IW,IV,IVD,R | bw·IW/N; Sabin and cVDPV terms | gam_p 1/28/d; R0_w 6; rel_v 0.4; rel_vd 0.9; p_rev 0.01; p_opv 0.8; take_v 0.5; par_ratio 1/200 | Fine & Carneiro 1999 AJE; Pons-Salort/Grassly (R) | I |
| rotavirus_pitzer2009 | Rotavirus, graded reinfections | M,S0,I1,R1,S1,I2,R2,S2,I3,R3,S3 | b0·(1+amp_s·cos(2π(t−phi)/52))·(I1+rho_inf·(I2+I3))/N | omega_m 1/13/wk (V via secondary); sus_1 0.62; sus_2 0.37 (V, Velázquez via [PMC9746532](https://pmc.ncbi.nlm.nih.gov/articles/PMC9746532)); gam_1 1/wk; gam_2 2/wk; R0 15; amp 0.15 | Pitzer et al. 2009 Science ([Princeton](https://collaborate.princeton.edu/en/publications/demographic-variability-vaccination-and-the-spatiotemporal-dynami/)) | partial V |
| norovirus_simmons2013 | Norovirus SEIAR with asymptomatic reinfection | S,E,I,A,R | b_nv·(I+eps_ea·(E+A))/N | inc_d 1 d; sym_d 2 d; shed_d 10 d; eps_ea 0.05 (Model B); theta_d 5.1 yr; R0 1.7 (all V) | Simmons et al. 2013 EID ([CDC](https://wwwnc.cdc.gov/eid/article/19/8/13-0472_article)) | **V** (Model B) |
| hepatitisA_seir_vacc | Hepatitis A with vaccination | S,E,IA,IS,R,V | (R0_ha/dur_h)·(IS+rel_ia·IA)/N | inc_h 28 d; dur_h 21 d; p_sym_h 0.3; R0 2.5; eff_hv 0.95 | Van Effelterre et al. 2006 Clin Infect Dis (R) | I |
| schisto_snail_mwb | Schistosomiasis, mean worm burden + SEI snails | Wm,Sn,En,In (all not in N; H constant) | snails: beta_sn·0.5·Wm·H·phi_m; worms: lam_h·In [abs] | r_s 0.1/d; K_s 1e4; mu_n 1/60/d; sig_n 1/40/d; mu_i 1/10/d; mu_wm 1/(3.3·365)/d; k_agg 0.2; praziquantel kill 0.86 (R) | Macdonald 1965; Anderson & May 1991; Feng, Li & Milner 2002 ([Purdue](https://www.math.purdue.edu/~fengz/documents/publications/mathbio02.pdf)) | I |
| sth_mean_worm | Ascaris / Trichuris / hookworm | Mw, Lr (not in N) | beta_l·Lr [abs] | R0 3; mu_w 1/yr; mu_l 6/yr; k_agg 0.5; z_dd 0.93; fdd = (1+Mw(1−z)/k)^−(k+1) | Anderson & May 1991; Truscott et al. 2014 (R) | I |
| trachoma_sis_mda | Ocular chlamydia, mass azithromycin | S,I | (R0_tr/dur_tr)·I/N | dur 26 wk; R0 2; cover 0.8; eff 0.95; int 52 or 26 wk | Lietman et al. 1999 Nat Med ([abstract](https://www.curehunter.com/public/pubmed10229236.do)) | I |
| crypto_oocyst | Cryptosporidium (Giardia variant) | S,I,R; Oc | exp_c·(1−exp(−r_dose·vol_d·Oc·(1−rem_t))) | r_dose 0.004/oocyst (R, Haas 1996; Giardia 0.0199); gam_c 1/14/d; mu_o 0.05/d; sh_vol 2e-4 | Composite; Haas et al. 1996 (R) | I |
| lepto_rodent_env | Leptospirosis, rodent reservoir | Sh,Ih,Rh; Sr,Ir,Bl | humans: be_h·Bl/(K_b+Bl); rodents: bd_r·Ir/Nr + be_r·envB | b_r 1/180/d; K_r 1000; del_b 0.1/d; be_h 1e-4/d; gam_h 1/14/d | Holt, Davis & Leirs 2006 Acta Trop (R) | I |
| scabies_mda | Scabies, two-dose ivermectin | S,I,Ie | (R0_sc/dur_sc)·I/N | R0 1.5; dur 365 d; hatch 1/7/d; cover 0.85; eff 0.95; doses on days 100 and 110 | Lydeamore et al. 2019 Math Biosci (R) | I |

```text
cholera_codeco2001:    >S:mu_h*N[abs]; S,I,R>:mu_h; S>I:a*B/(K+B); I>R:r; >B:e*I[abs]; B>:-nb_mb
cholera_andrews_basu:  S>I:p_sym*lam; S>A:(1-p_sym)*lam; I>R:gam_c; I>:cfr_rate; A>R:gam_c; R>S:omega_r; S>V:nu_v*eff_v; V>S:omega_v; >W:xi_w*(I+rel_a*A)[abs]; W>:del_w
typhoid_pitzer2014:    >S1:mu_h*N[abs]; S1>I1:lam; I1>Car:theta_c*delta_i; I1>R:(1-theta_c)*delta_i; R>S2:omega_r; S2>I2:lam; I2>R:delta_i; humans>:mu_h; >W:xi_w*prev_w[abs]; W>:xi_w
polio_wpv_opv:         >S:mu_h*N*(1-p_opv*take_v)[abs]; >IV:mu_h*N*p_opv*take_v[abs]; S>IW:bw*IW/N; S>IV:bv*IV/N*(1-p_rev); S>IVD:bv*IV/N*p_rev+bvd*IVD/N; IW,IV,IVD>R:gam_p; R>S:omega_p; all>:mu_h
rotavirus_pitzer2009:  >M:mu_h*N[abs]; M>S0:omega_m; S0>I1:lam; I1>R1:gam_1; R1>S1:omega_r; S1>I2:sus_1*lam; I2>R2:gam_2; R2>S2:omega_r; S2>I3:sus_2*lam; I3>R3:gam_2; R3>S3:omega_r; S3>I3:sus_2*lam; all>:mu_h
norovirus_simmons2013: >S:mu_h*N[abs]; S>E:lam; E>I:1/inc_d; I>A:1/sym_d; A>R:1/shed_d; R>A:lam; R>S:1/theta_d; all>:mu_h
hepatitisA_seir_vacc:  >S:mu_h*N*(1-p_vacc*eff_hv)[abs]; >V:mu_h*N*p_vacc*eff_hv[abs]; S>E:lam; E>IA:(1-p_sym_h)/inc_h; E>IS:p_sym_h/inc_h; IA,IS>R:1/dur_h; all>:mu_h
schisto_snail_mwb:     >Sn:r_s*(1-Nsn/K_s)*(Sn+En)[abs]; Sn>En:beta_sn*miracidia; En>In:sig_n; Sn,En>:mu_n; In>:mu_i; >Wm:lam_h*In[abs]; Wm>:mu_wm+mu_h+tau_t
sth_mean_worm:         >Mw:beta_l*Lr[abs]; Mw>:mu_w+mu_h+tau_m; >Lr:lam_e*Mw*fdd[abs]; Lr>:mu_l
trachoma_sis_mda:      S>I:b_tr*I/N; I>S:1/dur_tr+mda_rate
crypto_oocyst:         >S:mu_h*N[abs]; S>I:lam; I>R:gam_c; R>S:omega_c; S,I,R>:mu_h; >Oc:sh_vol*I[abs]; Oc>:mu_o
lepto_rodent_env:      >Sr:max(0,b_r*Nr*(1-Nr/K_r))[abs]; Sr>Ir:bd_r*Ir/Nr+be_r*envB; Sr,Ir>:mu_r; >Bl:sh_r*Ir[abs]; Bl>:del_b; >Sh:mu_h*N[abs]; Sh>Ih:be_h*envB; Ih>Rh:gam_h; Ih>:cfr_h; Rh>Sh:omega_h; humans>:mu_h
scabies_mda:           S>I:b_sc*I/N; I>S:1/dur_sc; I>Ie:mda_k; Ie>I:hatch; Ie>S:mda_k
```

Three exclusions are deliberate. Lymphatic filariasis and onchocerciasis models (EPIFIL, LYMFASIM, ONCHOSIM, EPIONCHO) are vector-based, age-structured or individual-based. No standard published ODE model of Guinea worm could be verified. The worm-burden models (schistosomiasis, STH) track a per-host mean, not a population, so their compartments carry no N. The generator needs to accept a model whose N is an explicit parameter (`H`).

## 7. Zoonoses, wildlife disease and within-host dynamics: 17 specifications, where values are thinnest

This group spans multi-host spillover, environmental reservoirs fed by dead hosts, and the target-cell-limited family. It is also the group with the least verified numbers.

- **Zinsstag et al. 2009 (dog rabies).** Dog-to-dog β_d ≈ 0.0807 and dog-to-human β_dh ≈ 0.0002 were verified, both in km²·dog⁻¹·week⁻¹, with Re ≈ 1.01 at the start of observation. A single campaign reaching ≥70% coverage was reported to interrupt transmission to humans for at least 6 years ([PMC2728111](https://pmc.ncbi.nlm.nih.gov/articles/PMC2728111)). The area unit must be kept, or the verified rates lose their meaning.
- **Mitchell et al. 2008 (chytrid).** The exp(−φZ) density limit on saprobic zoospore growth was verified, along with the finding that longer zoospore persistence raises the risk of host extinction ([PMC2593721](https://pmc.ncbi.nlm.nih.gov/articles/PMC2593721)).
- **Within-host models.** Perelson 1996, Neumann 1998 and Baccam 2006 are the entries where a wrong number is most likely. Their values are widely quoted, and Baccam's have been rescaled inconsistently. A secondary table reports β ≈ 1.13e-6, δ = 0.1 h⁻¹, p ≈ 9.9e4 and c = 0.125 h⁻¹. It explains that the original δ (0.167 h⁻¹) implied an implausible 6-hour infected-cell lifetime ([Smith et al. parameter table](https://www.math.utah.edu/~adler/myreprints/smith2009atp.pdf)).

Two citation corrections came up. **Iwami, Takeuchi & Liu 2007 is in *Mathematical Biosciences***, not J Theor Biol ([review](https://www.math.miami.edu/~ruan/MyPapers/LiuRuanZhang-TheoryBiosci2015.pdf)). **Chowell et al. 2014 on MERS is in *Epidemics* 9:40–51** and treats the camel reservoir as exogenous introductions ([PMC4161258](https://pmc.ncbi.nlm.nih.gov/articles/PMC4161258)).

| ID | Disease / area | Compartments | Force of infection | Typical parameters (units) | Primary source | Prov. |
|---|---|---|---|---|---|---|
| ebola_chowell_2004 | Ebola, Congo 1995 / Uganda 2000 | S,E,I,R,D*,CC* | beta_t·I/N, with beta_t = beta0, then beta0·e^(−q(t−tau)) after tau | beta0 0.33/d; q 0.98/d; tau 130 d; k 1/5.3/d; gamma 1/5.61/d; cfr 0.67; R0≈1.83 | Chowell et al. 2004 JTB (doi 10.1016/j.jtbi.2004.03.006) | R |
| ebola_rivers_2014_intervention | Ebola, West Africa 2014 + vaccination | S,E,I,H,Fn,R,V | (beta_i·I+beta_h·H+beta_f·Fn)/N | beta_i 0.16; beta_h 0.062; beta_f 0.489 (/d); alpha 1/11; theta 0.248; gamma_h 1/4.12; gamma_i 1/10; gamma_d 1/8; gamma_dh 1/6; gamma_f 1/2.01; cfr 0.7 | Rivers et al. 2014 PLOS Currents (R); overlaps the existing Legrand model | R |
| nipah_bat_pig_human | Nipah, bat → pig → human | Sb,Ib,Rb,Sp,Ep,Ip,Rp (not in N); Sh,Eh,Ih,Hh,Rh,Dh* | pig: beta_pp·Ip/Np + beta_bp·Ib; human: beta_ph·Ip + (beta_hh·Ih+beta_hosp·Hh)/N | beta_pp 1.2/d; beta_hosp 0.15/d; sigma_h 1/9/d; cfr_h 0.7 (R); others I | Pulliam et al. 2012 J R Soc Interface; Nikolay et al. 2019 NEJM | I |
| mers_camel_human_hospital | MERS, camel spillover + hospital | Sc,Ic,Rc (not in N); S,E,Ii,Hh,R,D*,CZ*,CS* | zoonotic beta_ch·Ic + (beta_com·Ii+beta_hosp·Hh)/N | beta_com 0.04; beta_hosp 0.55; sigma 1/5.2; eta 1/5; gamma_i 1/10; gamma_h 1/13 (/d); cfr 0.4 (R) | Chowell et al. 2014 Epidemics ([PMC4161258](https://pmc.ncbi.nlm.nih.gov/articles/PMC4161258)) | R/I |
| h5n1_poultry_human_iwami | H5N1 with a mutant human-transmissible strain | Sa,Ia (not in N); S,Ii,Im,R,D* | birds beta_a·Ia; spillover beta_ah·Ia; mutant beta_m·Im (mass action) | ca 1000/d; mu_a 1/500/d; beta_a 2e-5; delta_a 0.4/d; beta_m 2e-7; nu 1e-4/d (all I) | Iwami, Takeuchi & Liu 2007 Math Biosci | I |
| mpox_human_animal_reservoir | Mpox, rodent reservoir + smallpox-vaccine immunity | Sr,Er,Ir (not in N); S,V,E,Ii,R,D* | beta_rh·Ir + beta_hh·Ii/N (V: ×(1−eps_v)) | beta_hh 0.035/d; sigma_h 1/8.5/d; gamma_h 1/21/d; cfr 0.035; eps_v 0.85 (R) | Bhunu & Mushayabasa 2011 IJAM ([PDF](https://www.iaeng.org/IJAM/issues_v41/issue_2/IJAM_41_2_09.pdf)) | I |
| smallpox_emergency_vaccination_response | Smallpox, traced vs mass vaccination | S,E,Ii,Q,V,R,D*,CV* | beta·Ii/N; vaccination capped by pmin(want, cap) | beta 0.105/d; sigma 1/14/d (R); theta 1/3/d; cfr 0.3 (R); eps_pre 0.95; eps_post 0.5; cap 2e5 doses/d | Kaplan, Craft & Wein 2002 PNAS (discrete-time original; ODE simplification) | I |
| brucellosis_sheep_human_environment | Brucellosis, sheep + environment → humans | Ss,Es,Is,Vs,W (not in N); Sh,Eh,Iah,Ich,CHh* | sheep beta_ss·Is + beta_sw·W; humans beta_hs·Is + beta_hw·W | sigma_h 4/yr; mu_h 0.006/yr (V, related literature); others I | Hou et al. 2013 Math Biosci ([record](https://read.qxmd.com/read/23313258/modeling-the-transmission-dynamics-of-sheep-brucellosis-in-inner-mongolia-autonomous-region-china)); Li et al. 2014 MBE | I |
| rabies_dog_human_vaccination | Dog rabies, N'Djaména | Sd,Ed,Id,Vd,Dd* (not in N); Sh,Eh,Ih,Dh* | beta_d·Id/area; beta_dh·Id/area | beta_d 0.0807, beta_dh 0.0002 km²/(dog·wk) (V); area 100 km²; rd 0.49; mu_d 1/0.6/wk; sigma_d 1/3.1/wk | Zinsstag et al. 2009 PNAS ([PMC2728111](https://pmc.ncbi.nlm.nih.gov/articles/PMC2728111)) | partial V |
| btb_cattle_badger | Bovine TB, cattle ↔ badger, test-and-slaughter | Sc,Ec,Ic,CRc*; Sb,Eb,Ib (not in N) | cattle beta_cc·Ic/N + beta_bc·(1+perturb·(cull_b>0))·Ib | beta_cc 0.4/yr; beta_bb 0.6/yr; test 1/yr; mu_c 0.25/yr; mu_b 0.3/yr (all I) | Cox et al. 2005 PNAS ([PMC1308899](https://pmc.ncbi.nlm.nih.gov/articles/PMC1308899)) | I |
| asf_wild_boar_carcass | African swine fever, wild boar | S,E,Ii,R; Cr (not in N) | beta_dir·Ii/N + beta_car·Cr | cfr 0.95 (R); decay_c 1/90/d; sigma 1/4/d; gamma 1/5/d; Kb 2000 | O'Neill et al. 2020 Sci Rep; Pepin et al. 2020 | I |
| cwd_deer_environmental | Chronic wasting disease, prions | S,Ii; Ev (not in N) | beta_dir·Ii + beta_env·Ev | alpha 0.55/yr (R); decay_e 0.05/yr; Kd 1500; bdeer 0.35/yr | Miller, Hobbs & Tavener 2006; Almberg et al. 2011 ([PMC3094393](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3094393/)) | I |
| chytrid_amphibian_zoospore | Bd chytrid in Bufo bufo | S,Ii; Zs (not in N) | beta_hh·Ii + beta_z·Zs | sapro = rz·Zs·exp(−phi·Zs) (form V); rz 0.1/d; clear_z 0.1/d; prod 100/d | Mitchell et al. 2008 Proc R Soc B ([PMC2593721](https://pmc.ncbi.nlm.nih.gov/articles/PMC2593721)) | form V, values I |
| hiv_target_cell_therapy | HIV-1 within host, RT inhibitor/PI | Tc,Ii,V | beta·(1−er)·V | delta 0.5/d; cV 3.07/d (R); lam 10; dT 0.01; beta 2.4e-5; p 1000 | Perelson et al. 1996 Science | R |
| influenza_within_host_teiv | Influenza A, eclipse phase | Tc,I1,I2,V | beta·V | Tc0 4e8; k 4/d; delta 2.4/d; cV 3/d; beta 2.7e-5; p 1.2e-2 (all R, conflicting) | Baccam et al. 2006 J Virol ([PMC1563736](https://pmc.ncbi.nlm.nih.gov/articles/PMC1563736)) | R (low) |
| hcv_neumann_interferon | HCV within host, interferon-α | Tc,Ii,V | beta·(1−ei)·V | cV 6.2/d; delta 0.14/d (range 0.01–0.5); eta 0.95 (R); tau_lag 0.3 d | Neumann et al. 1998 Science | R |
| hiv_ctl_nowak_bangham | Virus–CTL dynamics | Tc,Ii,V,Zc | beta·V; CTL killing pz·Zc | lam 10; dT 0.1; beta 0.002; aI 0.5; pz 1; kv 50; uV 5; cz 0.1; bz 0.05 (theory paper, I) | Nowak & Bangham 1996 Science | I |

```text
ebola_chowell_2004:    S>E:lambda; E>I:k; >CC:k*E[abs]  (see bug note); I>D:gamma*cfr; I>R:gamma*(1-cfr)
ebola_rivers_2014:     S>E:lambda; S>V:nu*eff; E>I:alpha; I>H:gamma_h*theta; I>Fn:gamma_d*cfr*(1-theta); I>R:gamma_i*(1-cfr)*(1-theta); H>Fn:gamma_dh*cfr; H>R:gamma_dh*(1-cfr); Fn>R:gamma_f
nipah_bat_pig_human:   >Sb:mu_b*Nb[abs]; Sb>Ib:lambda_b; Ib>Rb:gamma_b; Rb>Sb:wane_b; bats>:mu_b; >Sp:mu_p*Np[abs]; Sp>Ep:lambda_p; Ep>Ip:sigma_p; Ip>Rp:gamma_p; pigs>:mu_p; Sh>Eh:lambda_h; Eh>Ih:sigma_h; Ih>Hh:eta_h; Ih>Rh:gamma_h; Hh>Dh:gamma_h*cfr_h; Hh>Rh:gamma_h*(1-cfr_h)
mers_camel_human:      >Sc:mu_c*Nc[abs]; Sc>Ic:lambda_c; Ic>Rc:gamma_c; Rc>Sc:wane_c; camels>:mu_c; S>E:lambda_z+lambda_s; >CZ:lambda_z*S[abs]; >CS:lambda_s*S[abs]; E>Ii:sigma; Ii>Hh:eta; Ii>R:gamma_i; Hh>D:gamma_h*cfr; Hh>R:gamma_h*(1-cfr)
h5n1_iwami:            >Sa:ca[abs]; Sa>Ia:lambda_a; Sa>:mu_a+cull; Ia>:mu_a+delta_a+cull; >S:ch[abs]; S>Ii:lambda_ah; S>Im:lambda_m; Ii>Im:nu; Ii>:mu_h+delta_i; Im>R:gamma_m; Im>D:delta_m; Im,S,R>:mu_h
mpox_reservoir:        >Sr:br*Nr[abs]; Sr>Er:lambda_r; Er>Ir:sigma_r; Ir>:gamma_r; rodents>:mu_r; >S:bh[abs]; S>E:lambda_h; S>V:nu; V>E:lambda_h*(1-eps_v); V>S:wane_v; E>Ii:sigma_h; Ii>R:gamma_h*(1-cfr); Ii>D:gamma_h*cfr; S,V,E,Ii,R>:mu_h
smallpox_response:     S>E:lambda; S>V:rate_v*eps_pre; E>V:rate_v*eps_post; E>Ii:sigma; Ii>Q:theta; Q>D:gamma*cfr; Q>R:gamma*(1-cfr); >CV:pmin(want,cap)[abs]
brucellosis:           >Ss:As[abs]; Ss>Es:lambda_s; Ss>Vs:nu_s; Vs>Ss:wane_s; Es>Is:sigma_s; sheep>:mu_s; >W:shed*Is[abs]; W>:decay+disinf; Sh>Eh:lambda_h; Eh>Iah:sigma_h; Iah>Ich:chron; Iah>Sh:gamma_a; Ich>Sh:gamma_c; >CHh:sigma_h*Eh[abs]; humans>:mu_h
rabies_dog_human:      >Sd:bd*Nd[abs]; Sd>Ed:lambda_d; Sd>Vd:nu_d; Vd>Sd:wane_d; Ed>Id:sigma_d*rd; Ed>Sd:sigma_d*(1-rd); Id>Dd:mu_d; Sd,Ed,Vd>:md+cull_d; Sh>Eh:lambda_h; Eh>Ih:sigma_h*rh*(1-pep); Eh>Sh:sigma_h*(1-rh*(1-pep)); Ih>Dh:mu_h
btb_cattle_badger:     >Sc:mu_c*N+test*Ic[abs]; Sc>Ec:lambda_c; Ec>Ic:sigma_c; Ic>CRc:test; cattle>:mu_c; >Sb:mu_b*Nb+alpha_b*Ib+cull_b*Nb[abs]; Sb>Eb:lambda_b; Eb>Ib:sigma_b; Sb,Eb>:mu_b+cull_b; Ib>:mu_b+alpha_b+cull_b
asf_wild_boar:         >S:pmax(dens,0)*N[abs]; S>E:lambda; E>Ii:sigma; Ii>Cr:gamma*cfr; Ii>R:gamma*(1-cfr); S,E,Ii,R>:mub+hunt; Cr>:decay_c+remove_c
cwd_deer:              >S:births[abs]; S>Ii:lambda; S>:mu+harvest; Ii>:mu+alpha+harvest; >Ev:shed*Ii+carc*alpha*Ii[abs]; Ev>:decay_e
chytrid:               >S:bz*N*pmax(1-N/Kz,0)[abs]; S>Ii:lambda; Ii>S:rec; S>:mu_z; Ii>:mu_z+alpha_z; >Zs:prod*Ii+sapro[abs]; Zs>:clear_z
hiv_target_cell:       >Tc:lam[abs]; Tc>:dT; Tc>Ii:beta*(1-er)*V; Ii>:delta; >V:p*(1-ep)*Ii[abs]; V>:cV
influenza_teiv:        Tc>I1:beta*V; I1>I2:k; I2>:delta; >V:p*I2[abs]; V>:cV
hcv_neumann:           >Tc:sh[abs]; Tc>:dh; Tc>Ii:beta*(1-ei)*V; Ii>:delta; >V:p*(1-en)*Ii[abs]; V>:cV
hiv_ctl_nowak_bangham: >Tc:lam[abs]; Tc>:dT; Tc>Ii:beta*V; Ii>:aI+pz*Zc; >V:kv*Ii[abs]; V>:uV; >Zc:cz*Ii*Zc[abs]; Zc>:bz
```

The environmental-reservoir pattern is the reusable lesson here. ASF, CWD, chytrid, brucellosis, cholera, leptospirosis and Cryptosporidium all route shedding or dead hosts into a long-lived pathogen pool outside N, and the ratio of reservoir decay to infectious period controls persistence. Q fever and anthrax (Saad-Roy, van den Driessche & Yakubu 2017) are the cheapest additions because they share this structure. They were left out for budget reasons, not because they fail the format. The Keeling et al. 2001 FMD model was excluded on principle: it is a stochastic, spatial, farm-level simulation, and a small ODE would not be that model.

## 8. Specification bugs, duplicates and source conflicts to resolve before code generation

A cross-check of the flow listings found **three accumulator flows written as transfers**. These would drain a real compartment twice. In `ebola_chowell_2004` the note writes `E > CC : k` alongside `E > I : k`. In `mers_camel_human_hospital` it writes `S > CZ : lambda_z` and `S > CS : lambda_s` alongside `S > E : lambda_z + lambda_s`. In `brucellosis_sheep_human_environment` it writes `Eh > CHh : sigma_h` alongside `Eh > Iah : sigma_h`. Each must become an absolute inflow to the counter (` > CC : k*E`, ` > CZ : lambda_z*S`, ` > CS : lambda_s*S`, ` > CHh : sigma_h*Eh`). The listings above already show the corrected form. In the same brucellosis entry, the derived `R0s` expression is garbled (it contains a `shed * 0 + 1e-12` placeholder) and should be re-derived symbolically. In `covid_vax_waning`, `vh_red` is declared but unused. The `sidr_diagnosis` R0 works out to about 1.49 with the stated values, not the 1.53 in the note.

There are **three overlapping pairs**:

- `two_strain_cross_immunity` and `flu_two_strain_cross` are the same structure, differing only by waning. Ship one model with an `omega` switch.
- `hiv_test_and_treat` and `hiv_art_tasp_granich` are a 1-stage and a 4-stage Granich variant.
- `sis_core_group_gonorrhea` and `gonorrhoea_core_group_hethcote_yorke` are a single-sex and a two-sex Hethcote–Yorke variant.

`ebola_rivers_2014_intervention` reuses the already-implemented Legrand 2007 skeleton. It should ship only as "Legrand structure, 2014 parameters, plus a vaccination lever", or be dropped.

Five source conflicts affect validation targets:

- SIDARTHE's arXiv v1 values differ from the *Nature Medicine* version.
- Chowell's EID 2004 R0 values (0.58–1.17) are probably reproduction numbers under control.
- Baccam 2006 parameters appear in incompatible per-hour and per-day rescalings.
- Simmons 2013 reports two different Model C confidence intervals.
- The Arino 2008 Table 1 assignments of θ and σ rest on row order.

Several checkable R0 values come out right.

- Bowman's mosquito–bird threshold is about 2.9 with Π_M = 20,000 and μ_M = 0.05.
- Yakob & Clements gives about 4.1.
- Chowell Ebola gives β0/γ = 0.33 × 5.61 ≈ 1.85, against the reported 1.83.
- The saturated-treatment model gives about 1.49.
- Codeço's S_crit reproduces exactly.

These make good automated regression tests for the generated R code.

## Conclusion

The catalogue changes what the hard part of the library is. Expressing models is not the bottleneck: every family tried, from saturating biting and superinfection to dose–response, worm burdens and CTL bilinearity, fits a flows-plus-derived-variables grammar. The few requirements are absolute flows, accumulators kept out of N, `ifelse`-based pulses and Erlang chains. The binding constraint is **parameter provenance**. Roughly one model in seven has a source-read parameter set, and the most-cited within-host numbers are the least reliable. The library therefore needs a first-class metadata field per parameter (`verified` / `recalled` / `illustrative`, with source page or table). That field matters as much as the model code, and the verified subset should be the default "reference" models, with published R0 values as unit tests.

The cheapest next gains are open-access tables that were identified but not transcribed. These are Gaff 2007 (EJDE), Stauch 2011 (PMC3226461), Keeling & Gilligan 2000 (PMC1690796), Gumel 2004 Table 1, the full texts of Perelson 1996, Neumann 1998 and Baccam 2006, and the supplements of Whittles 2017, Martin 2013, Pitzer 2009/2014 and Simmons 2013. Pulling these would roughly double the verified share without adding a single new structure.
