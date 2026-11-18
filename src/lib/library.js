// Library of standard compartmental models.
//
// Each entry is written in a compact notation and expanded into the editor's
// model format by def(). Flows are written "FROM > TO : rate | options":
//   "S > I : lam"          per-capita flow  (lam * S)
//   " > S : mu * N"        inflow from outside (absolute)
//   "I >   : mu"           outflow to outside (per capita)
//   options: out=<name> (record flow as output), curve=<px>, abs (absolute rate),
//            at=<dx>,<dy> (position of the outside end relative to the compartment)
// Compartments: [name, description, init, col, row, {inN, color}]
// Parameters:   [name, value, description, symbol?]
// Variables:    [name, expr, description, output?, symbol?]
//
// Parameter values are typical textbook values chosen for illustration; check
// them against your own setting before using any model for inference.

const COLORS = {
  S: '#FFFFFF', M: '#F3C6D3', E: '#E8C468', L: '#E8C468', P: '#F0A35E', I: '#B71C1C', A: '#E57373',
  C: '#8E0E0E', R: '#4E8A3E', W: '#A5D6A7', V: '#2457C5', Q: '#00838F', H: '#6A2C9C', F: '#5D4037', D: '#111111',
};

function colorFor(name) {
  return COLORS[name] || COLORS[name[0]] || '#DDDDDD';
}

function def(s) {
  const comps = s.comps.map((c) => (Array.isArray(c) ? c : [c.name, c.description, c.init, c.col, c.row, { inN: c.inN, color: c.color }]));
  const compartments = comps.map(([name, description, init, col, row, opts = {}], i) => ({
    id: 'c' + (i + 1),
    name,
    description,
    init: String(init),
    inN: opts.inN !== false,
    color: opts.color || colorFor(name),
    x: 120 + col * 200,
    y: 140 + row * 160,
  }));
  const byName = Object.fromEntries(compartments.map((c) => [c.name, c]));
  const nIn = {};
  const nOut = {};
  const flows = s.flows.map((spec, i) => {
    if (typeof spec === 'object') {
      spec = `${spec.from || ''} > ${spec.to || ''} : ${spec.rate.replace(/\|/g, '')} |${spec.absolute ? ' abs' : ''}${spec.curve ? ' curve=' + spec.curve : ''}`;
    }
    const [main0, optStr = ''] = spec.split('|');
    const main = main0.replace(//g, '|');
    const m = /^\s*([A-Za-z0-9_.]*)\s*>\s*([A-Za-z0-9_.]*)\s*:\s*(.+?)\s*$/.exec(main);
    if (!m) throw new Error(`${s.id}: bad flow "${spec}"`);
    const [, a, b, rate] = m;
    if ((a && !byName[a]) || (b && !byName[b])) throw new Error(`${s.id}: unknown compartment in "${spec}"`);
    const opts = {};
    optStr.trim().split(/\s+/).filter(Boolean).forEach((o) => { const [k, v = true] = o.split('='); opts[k] = v; });
    const f = {
      id: 'f' + (i + 1),
      from: a ? byName[a].id : null,
      to: b ? byName[b].id : null,
      kind: !a || opts.abs ? 'absolute' : 'per_capita',
      rate,
      output: opts.out || '',
      curve: Number(opts.curve || 0),
    };
    if (!a || !b) {
      const node = byName[a || b];
      let dx;
      let dy;
      if (opts.at) [dx, dy] = opts.at.split(',').map(Number);
      else if (!a) { const k = (nIn[b] = (nIn[b] || 0) + 1) - 1; dx = k * 50; dy = -130; }
      else { const k = (nOut[a] = (nOut[a] || 0) + 1) - 1; dx = k * 50; dy = 130; }
      f.ox = node.x + dx;
      f.oy = node.y + dy;
    }
    return f;
  });
  return {
    version: 1,
    id: s.id,
    category: s.category,
    name: s.name,
    description: s.description,
    reference: s.reference || '',
    referenceUrl: s.url || '',
    keywords: s.keywords || [],
    time: { start: String(s.time[0]), stop: String(s.time[1]), step: String(s.time[2]), unit: s.time[3] },
    population: { name: s.pop || 'N' },
    ...(s.hmax ? { solver: { hmax: String(s.hmax) } } : {}),
    compartments,
    parameters: (s.params || []).map((p) => (Array.isArray(p) ? p : [p.name, p.value, p.description, p.symbol]))
      .map(([name, value, description, symbol = '']) => ({ name, symbol: symbol || '', value: String(value), description })),
    variables: (s.vars || []).map((v) => (Array.isArray(v) ? v : [v.name, v.expr, v.description, v.output, v.symbol]))
      .map(([name, expr, description, output = false, symbol = '']) => ({ name, symbol: symbol || '', expr, description, output: !!output })),
    flows,
  };
}

import RESEARCH from './library-research.js';
import { validate } from './generator.js';

const KR = 'Keeling & Rohani (2008), Modeling Infectious Diseases in Humans and Animals';
const VW = 'Vynnycky & White (2010), An Introduction to Infectious Disease Modelling';
const AM = 'Anderson & May (1991), Infectious Diseases of Humans';

export const LIBRARY = [
  // -------------------------------------------------------------------------
  // Basic
  // -------------------------------------------------------------------------
  def({
    id: 'si', category: 'Basic', name: 'SI',
    description: 'Infection with no recovery (e.g. untreated chronic infection) in a closed population.',
    reference: KR + ', ch. 2',
    time: [0, 100, 1, 'days'],
    comps: [['S', 'Susceptible', 999, 0, 0], ['I', 'Infected', 1, 1, 0]],
    params: [['beta', 0.2, 'transmission rate per day']],
    vars: [['lam', 'beta * I / N', 'force of infection']],
    flows: ['S > I : lam | out=incidence'],
  }),
  def({
    id: 'sis', category: 'Basic', name: 'SIS',
    description: 'Infection without lasting immunity (e.g. many bacterial STIs): recovered people become susceptible again.',
    reference: KR + ', §2.4',
    time: [0, 200, 1, 'days'],
    comps: [['S', 'Susceptible', 999, 0, 0], ['I', 'Infectious', 1, 1, 0]],
    params: [['beta', 0.3, 'transmission rate per day'], ['gamma', '1 / 10', 'recovery rate = 1/(infectious period)']],
    vars: [['lam', 'beta * I / N', 'force of infection']],
    flows: ['S > I : lam | curve=40', 'I > S : gamma | curve=40'],
  }),
  def({
    id: 'sir', category: 'Basic', name: 'SIR',
    description: 'Basic SIR epidemic in a closed population.',
    reference: 'Kermack & McKendrick (1927); ' + KR + ', §2.1',
    time: [0, 160, 1, 'days'],
    comps: [['S', 'Susceptible', 990, 0, 0], ['I', 'Infectious', 10, 1, 0], ['R', 'Recovered', 0, 2, 0]],
    params: [['beta', 0.3, 'transmission rate per day'], ['gamma', '1 / 10', 'recovery rate = 1/(infectious period)']],
    vars: [['lam', 'beta * I / N', 'force of infection']],
    flows: ['S > I : lam | out=incidence', 'I > R : gamma'],
  }),
  def({
    id: 'sir_vital', category: 'Basic', name: 'SIR with births and deaths',
    description: 'Endemic SIR with demography; damped oscillations to an endemic equilibrium (measles-like parameters).',
    reference: KR + ', §2.1.2 (Program 2.2)',
    time: [0, 60, 0.02, 'years'],
    comps: [['S', 'Susceptible', 10000, 0, 0], ['I', 'Infectious', 10, 1, 0], ['R', 'Recovered', 89990, 2, 0]],
    params: [
      ['beta', 520, 'transmission rate per year'],
      ['gamma', '365 / 7', 'recovery rate per year (7-day infectious period)'],
      ['mu', '1 / 70', 'birth = death rate per year (70-year life expectancy)'],
    ],
    vars: [['lam', 'beta * I / N', 'force of infection']],
    flows: [' > S : mu * N', 'S > I : lam | out=incidence', 'I > R : gamma', 'S > : mu', 'I > : mu', 'R > : mu'],
  }),
  def({
    id: 'sirs', category: 'Basic', name: 'SIRS',
    description: 'SIR with waning immunity: recovered individuals return to the susceptible class.',
    reference: KR + ', §2.4',
    time: [0, 1000, 1, 'days'],
    comps: [['S', 'Susceptible', 990, 0, 0], ['I', 'Infectious', 10, 1, 0], ['R', 'Recovered', 0, 2, 0]],
    params: [
      ['beta', 0.3, 'transmission rate per day'],
      ['gamma', '1 / 10', 'recovery rate'],
      ['omega', '1 / 180', 'rate of loss of immunity = 1/(duration of immunity)'],
    ],
    vars: [['lam', 'beta * I / N', 'force of infection']],
    flows: ['S > I : lam', 'I > R : gamma', 'R > S : omega | curve=90'],
  }),
  def({
    id: 'seir', category: 'Basic', name: 'SEIR',
    description: 'Epidemic with a latent (infected, not yet infectious) period, closed population.',
    reference: KR + ', §2.5',
    time: [0, 200, 1, 'days'],
    comps: [['S', 'Susceptible', 9999, 0, 0], ['E', 'Exposed (latent)', 1, 1, 0], ['I', 'Infectious', 0, 2, 0], ['R', 'Recovered', 0, 3, 0]],
    params: [
      ['beta', 0.5, 'transmission rate per day'],
      ['sigma', '1 / 5', 'rate of becoming infectious = 1/(latent period)'],
      ['gamma', '1 / 7', 'recovery rate = 1/(infectious period)'],
    ],
    vars: [['lam', 'beta * I / N', 'force of infection']],
    flows: ['S > E : lam', 'E > I : sigma | out=onsets', 'I > R : gamma'],
  }),
  def({
    id: 'seir_vital', category: 'Basic', name: 'SEIR with births and deaths',
    description: 'Endemic SEIR parameterised by R0 (measles-like): transmission rate derived from R0.',
    reference: KR + ', §2.5 (Program 2.6)',
    time: [0, 50, 0.01, 'years'],
    comps: [['S', 'Susceptible', 6000, 0, 0], ['E', 'Exposed', 10, 1, 0], ['I', 'Infectious', 10, 2, 0], ['R', 'Recovered', 93980, 3, 0]],
    params: [
      ['R0', 13, 'basic reproduction number'],
      ['sigma', '365 / 8', 'rate of becoming infectious per year (8-day latent period)'],
      ['gamma', '365 / 5', 'recovery rate per year (5-day infectious period)'],
      ['mu', '1 / 70', 'birth = death rate per year'],
    ],
    vars: [
      ['beta', 'R0 * (gamma + mu) * (sigma + mu) / sigma', 'transmission rate implied by R0'],
      ['lam', 'beta * I / N', 'force of infection'],
    ],
    flows: [' > S : mu * N', 'S > E : lam', 'E > I : sigma | out=incidence', 'I > R : gamma', 'S > : mu', 'E > : mu', 'I > : mu', 'R > : mu'],
  }),
  def({
    id: 'seirs', category: 'Basic', name: 'SEIRS',
    description: 'SEIR with waning immunity, giving recurrent epidemics.',
    reference: KR + ', §2.5',
    time: [0, 1500, 1, 'days'],
    comps: [['S', 'Susceptible', 9999, 0, 0], ['E', 'Exposed', 1, 1, 0], ['I', 'Infectious', 0, 2, 0], ['R', 'Recovered', 0, 3, 0]],
    params: [
      ['beta', 0.5, 'transmission rate per day'],
      ['sigma', '1 / 5', '1/(latent period)'],
      ['gamma', '1 / 7', '1/(infectious period)'],
      ['omega', '1 / 365', '1/(duration of immunity)'],
    ],
    vars: [['lam', 'beta * I / N', 'force of infection']],
    flows: ['S > E : lam', 'E > I : sigma', 'I > R : gamma', 'R > S : omega | curve=110'],
  }),
  def({
    id: 'sird', category: 'Basic', name: 'SIRD',
    description: 'SIR with disease-induced deaths. D is excluded from the living population N.',
    reference: KR + ', §2.2',
    time: [0, 200, 1, 'days'],
    comps: [['S', 'Susceptible', 9990, 0, 0], ['I', 'Infectious', 10, 1, 0], ['R', 'Recovered', 0, 2, 0], ['D', 'Dead', 0, 2, 1, { inN: false }]],
    params: [
      ['beta', 0.3, 'transmission rate per day'],
      ['gamma', '1 / 10', 'rate of leaving the infectious state'],
      ['cfr', 0.02, 'case fatality ratio'],
    ],
    vars: [['lam', 'beta * I / N', 'force of infection']],
    flows: ['S > I : lam', 'I > R : gamma * (1 - cfr)', 'I > D : gamma * cfr | out=deaths'],
  }),
  def({
    id: 'seird', category: 'Basic', name: 'SEIRD',
    description: 'SEIR with disease-induced deaths.',
    reference: KR + ', §2.5',
    time: [0, 250, 1, 'days'],
    comps: [['S', 'Susceptible', 99990, 0, 0], ['E', 'Exposed', 10, 1, 0], ['I', 'Infectious', 0, 2, 0], ['R', 'Recovered', 0, 3, 0], ['D', 'Dead', 0, 3, 1, { inN: false }]],
    params: [
      ['beta', 0.4, 'transmission rate per day'],
      ['sigma', '1 / 5', '1/(latent period)'],
      ['gamma', '1 / 8', 'rate of leaving the infectious state'],
      ['cfr', 0.01, 'case fatality ratio'],
    ],
    vars: [['lam', 'beta * I / N', 'force of infection']],
    flows: ['S > E : lam', 'E > I : sigma | out=incidence', 'I > R : gamma * (1 - cfr)', 'I > D : gamma * cfr | out=deaths'],
  }),

  // -------------------------------------------------------------------------
  // Interventions
  // -------------------------------------------------------------------------
  def({
    id: 'sir_newborn_vacc', category: 'Interventions', name: 'SIR with newborn vaccination',
    description: 'A proportion p of newborns is vaccinated (enters V). Elimination requires p > 1 - 1/R0.',
    reference: KR + ', §8.1; ' + AM + ', ch. 5',
    time: [0, 80, 0.02, 'years'],
    comps: [['S', 'Susceptible', 10000, 0, 0], ['I', 'Infectious', 10, 1, 0], ['R', 'Recovered', 89990, 2, 0], ['V', 'Vaccinated', 0, 0, 1]],
    params: [
      ['R0', 15, 'basic reproduction number'],
      ['gamma', '365 / 7', 'recovery rate per year'],
      ['mu', '1 / 70', 'birth = death rate per year'],
      ['p', 0.8, 'proportion of newborns vaccinated'],
    ],
    vars: [
      ['beta', 'R0 * (gamma + mu)', 'transmission rate implied by R0'],
      ['lam', 'beta * I / N', 'force of infection'],
    ],
    flows: [' > S : (1 - p) * mu * N', ' > V : p * mu * N | at=-150,0', 'S > I : lam | out=incidence', 'I > R : gamma', 'S > : mu | at=-60,110', 'I > : mu', 'R > : mu', 'V > : mu | at=0,130'],
  }),
  def({
    id: 'sirv_campaign', category: 'Interventions', name: 'SIR with a vaccination campaign',
    description: 'Mass vaccination of susceptibles at rate v starting on day t_vac (assumes a perfect vaccine).',
    reference: KR + ', §8.1',
    time: [0, 200, 1, 'days'],
    comps: [['S', 'Susceptible', 9990, 0, 0], ['I', 'Infectious', 10, 1, 0], ['R', 'Recovered', 0, 2, 0], ['V', 'Vaccinated', 0, 0, 1]],
    params: [
      ['beta', 0.3, 'transmission rate per day'],
      ['gamma', '1 / 10', 'recovery rate'],
      ['v', 0.02, 'vaccination rate per day once the campaign starts'],
      ['t_vac', 30, 'day the campaign starts'],
    ],
    vars: [
      ['lam', 'beta * I / N', 'force of infection'],
      ['vacc', 'v * (t >= t_vac)', 'vaccination rate at time t'],
    ],
    flows: ['S > I : lam | out=incidence', 'I > R : gamma', 'S > V : vacc'],
  }),
  def({
    id: 'leaky_vaccine', category: 'Interventions', name: 'SVIR with a leaky vaccine',
    description: 'Vaccinated individuals are infected at a reduced rate (1 - eff) * lam.',
    reference: VW + ', ch. 7',
    time: [0, 300, 1, 'days'],
    comps: [['S', 'Susceptible', 3999, 0, 0], ['V', 'Vaccinated', 6000, 0, 1], ['I', 'Infectious', 1, 1, 0], ['R', 'Recovered', 0, 2, 0]],
    params: [
      ['beta', 0.4, 'transmission rate per day'],
      ['gamma', '1 / 7', 'recovery rate'],
      ['eff', 0.6, 'vaccine efficacy against infection'],
    ],
    vars: [['lam', 'beta * I / N', 'force of infection']],
    flows: ['S > I : lam', 'V > I : (1 - eff) * lam | out=breakthrough', 'I > R : gamma'],
  }),
  def({
    id: 'siqr', category: 'Interventions', name: 'SIQR (isolation)',
    description: 'Infectious individuals are detected and isolated (Q) at rate q and no longer transmit.',
    reference: 'Hethcote, Zhien & Shengbing (2002), Math. Biosci. 180:141-160',
    time: [0, 200, 1, 'days'],
    comps: [['S', 'Susceptible', 9990, 0, 0], ['I', 'Infectious', 10, 1, 0], ['Q', 'Isolated', 0, 1, 1], ['R', 'Recovered', 0, 2, 0]],
    params: [
      ['beta', 0.4, 'transmission rate per day'],
      ['gamma', '1 / 7', 'recovery rate'],
      ['q', 0.1, 'rate of detection and isolation'],
    ],
    vars: [['lam', 'beta * I / N', 'force of infection']],
    flows: ['S > I : lam', 'I > R : gamma', 'I > Q : q', 'Q > R : gamma'],
  }),
  def({
    id: 'seir_lockdown', category: 'Interventions', name: 'SEIR with a temporary lockdown',
    description: 'Contacts fall by a proportion eff between days t_on and t_off.',
    reference: 'Flaxman et al. (2020), Nature 584:257-261 (concept)',
    time: [0, 365, 1, 'days'],
    comps: [['S', 'Susceptible', 999990, 0, 0], ['E', 'Exposed', 10, 1, 0], ['I', 'Infectious', 0, 2, 0], ['R', 'Recovered', 0, 3, 0]],
    params: [
      ['R0', 2.5, 'basic reproduction number'],
      ['sigma', '1 / 5', '1/(latent period)'],
      ['gamma', '1 / 5', '1/(infectious period)'],
      ['eff', 0.6, 'proportional reduction in contacts during lockdown'],
      ['t_on', 50, 'day lockdown starts'],
      ['t_off', 120, 'day lockdown ends'],
    ],
    vars: [
      ['contact', '1 - eff * (t >= t_on & t < t_off)', 'relative contact rate'],
      ['lam', 'contact * R0 * gamma * I / N', 'force of infection'],
      ['Rt', 'contact * R0 * S / N', 'effective reproduction number', true],
    ],
    flows: ['S > E : lam', 'E > I : sigma | out=incidence', 'I > R : gamma'],
  }),

  // -------------------------------------------------------------------------
  // Natural history
  // -------------------------------------------------------------------------
  def({
    id: 'msir', category: 'Natural history', name: 'MSIR (maternal immunity)',
    description: 'Newborns are protected by maternal antibodies (M) before becoming susceptible.',
    reference: 'Hethcote (2000), SIAM Review 42:599-653',
    time: [0, 60, 0.02, 'years'],
    comps: [['M', 'Maternally immune', 500, 0, 0], ['S', 'Susceptible', 9500, 1, 0], ['I', 'Infectious', 10, 2, 0], ['R', 'Recovered', 89990, 3, 0]],
    params: [
      ['beta', 700, 'transmission rate per year'],
      ['delta', 2, 'rate of losing maternal immunity per year (6 months)'],
      ['gamma', '365 / 7', 'recovery rate per year'],
      ['mu', '1 / 70', 'birth = death rate per year'],
    ],
    vars: [['lam', 'beta * I / N', 'force of infection']],
    flows: [' > M : mu * N', 'M > S : delta', 'S > I : lam | out=incidence', 'I > R : gamma', 'M > : mu', 'S > : mu', 'I > : mu', 'R > : mu'],
  }),
  def({
    id: 'seiar', category: 'Natural history', name: 'SEIAR (asymptomatic infection)',
    description: 'A proportion of infections is asymptomatic (A) and less infectious.',
    reference: 'Arino et al. (2008), J. R. Soc. Interface 5(Suppl 2)',
    time: [0, 200, 1, 'days'],
    comps: [['S', 'Susceptible', 99990, 0, 0.5], ['E', 'Exposed', 10, 1, 0.5], ['I', 'Symptomatic', 0, 2, 0], ['A', 'Asymptomatic', 0, 2, 1], ['R', 'Recovered', 0, 3, 0.5]],
    params: [
      ['beta', 0.5, 'transmission rate per day'],
      ['rA', 0.5, 'relative infectiousness of asymptomatics', 'r_A'],
      ['sigma', '1 / 4', '1/(latent period)'],
      ['ps', 0.6, 'proportion of infections that are symptomatic', 'p_s'],
      ['gamma', '1 / 6', 'recovery rate'],
    ],
    vars: [['lam', 'beta * (I + rA * A) / N', 'force of infection']],
    flows: ['S > E : lam', 'E > I : sigma * ps | out=cases', 'E > A : sigma * (1 - ps)', 'I > R : gamma', 'A > R : gamma'],
  }),
  def({
    id: 'covid_sepiar', category: 'Natural history', name: 'COVID-19 style SEPIR',
    description: 'Latent (E), pre-symptomatic infectious (P), then symptomatic (Is) or asymptomatic (Ia) infection.',
    reference: 'Davies et al. (2020), Nature Medicine 26:1205-1211 (structure)',
    time: [0, 300, 0.5, 'days'],
    comps: [
      ['S', 'Susceptible', 999990, 0, 0.5], ['E', 'Latent', 10, 1, 0.5], ['P', 'Pre-symptomatic', 0, 2, 0.5],
      ['Is', 'Symptomatic', 0, 3, 0, { color: '#B71C1C' }], ['Ia', 'Asymptomatic', 0, 3, 1, { color: '#E57373' }], ['R', 'Recovered', 0, 4, 0.5],
    ],
    params: [
      ['R0', 2.5, 'basic reproduction number'],
      ['sigma', '1 / 3', '1/(latent period)'],
      ['delta', '1 / 2', '1/(pre-symptomatic period)'],
      ['gamma', '1 / 5', '1/(infectious period after the pre-symptomatic stage)'],
      ['ps', 0.6, 'proportion symptomatic', 'p_s'],
      ['rP', 1, 'relative infectiousness, pre-symptomatic', 'r_P'],
      ['rA', 0.5, 'relative infectiousness, asymptomatic', 'r_A'],
    ],
    vars: [
      ['beta', 'R0 / (rP / delta + (ps + (1 - ps) * rA) / gamma)', 'transmission rate implied by R0'],
      ['lam', 'beta * (rP * P + Is + rA * Ia) / N', 'force of infection'],
    ],
    flows: ['S > E : lam', 'E > P : sigma', 'P > Is : delta * ps | out=onsets', 'P > Ia : delta * (1 - ps)', 'Is > R : gamma', 'Ia > R : gamma'],
  }),
  def({
    id: 'seihrd', category: 'Natural history', name: 'SEIHRD (hospitalisation)',
    description: 'A proportion of cases is hospitalised; deaths occur in hospital.',
    reference: VW + ', ch. 4 (structure)',
    time: [0, 300, 1, 'days'],
    comps: [
      ['S', 'Susceptible', 99990, 0, 0], ['E', 'Exposed', 10, 1, 0], ['I', 'Infectious', 0, 2, 0], ['R', 'Recovered', 0, 3, 0],
      ['H', 'Hospitalised', 0, 2, 1], ['D', 'Dead', 0, 3, 1, { inN: false }],
    ],
    params: [
      ['beta', 0.45, 'transmission rate per day'],
      ['sigma', '1 / 5', '1/(latent period)'],
      ['gamma', '1 / 7', 'rate of leaving the infectious state'],
      ['h', 0.05, 'proportion hospitalised'],
      ['tauH', '1 / 10', 'rate of leaving hospital'],
      ['m', 0.2, 'in-hospital fatality'],
    ],
    vars: [['lam', 'beta * I / N', 'force of infection (hospitalised patients do not transmit)']],
    flows: ['S > E : lam', 'E > I : sigma', 'I > R : gamma * (1 - h)', 'I > H : gamma * h | out=admissions', 'H > R : tauH * (1 - m)', 'H > D : tauH * m | out=deaths'],
  }),
  def({
    id: 'seachr', category: 'Natural history', name: 'SEACHR with re-exposure',
    description: 'Asymptomatic (A) and clinical (C) infection, hospitalisation (H), death (D),\nwaning immunity and re-exposure of recovereds (ER). Parameter values are illustrative;\nthe force of infection is an assumption (the diagram does not define it).',
    reference: 'Custom diagram',
    time: [0, 365, 1, 'days'],
    comps: [
      ['S', 'Susceptible', '100000 - E', 0, 1], ['E', 'Exposed', 10, 0.95, 1], ['A', 'Asymptomatic', 0, 2, 0.4, { color: '#C62828' }],
      ['C', 'Clinical', 0, 2, 1.6, { color: '#8E0E0E' }], ['R', 'Recovered', 0, 3.2, 1], ['ER', 'Exposed recovered', 0, 3.2, -0.2, { color: '#D9A441' }],
      ['H', 'Hospitalised', 0, 2, 2.85], ['D', 'Dead', 0, 3.2, 2.85, { inN: false }],
    ],
    params: [
      ['beta', 0.35, 'transmission rate per day'],
      ['rA', 0.5, 'relative infectiousness of asymptomatics', 'r_A'],
      ['nu', '1 / 9', '1/(latent period)'],
      ['gamma', 0.7, 'proportion of infections that are clinical'],
      ['tauA', '1 / 7', 'recovery rate of asymptomatics'],
      ['tauC', '1 / 5', 'rate of leaving the clinical stage'],
      ['theta', 0.3, 'proportion of clinical cases hospitalised'],
      ['tauH', '1 / 10', 'rate of leaving hospital'],
      ['omega', 0.4, 'case fatality among hospitalised'],
      ['sigma', 0.2, 'relative susceptibility of recovereds'],
      ['alpha', '1 / 365', 'rate of waning immunity'],
    ],
    vars: [['lam', 'beta * (rA * A + C) / N', 'force of infection (assumed)']],
    flows: [
      'S > E : lam', 'E > A : nu * (1 - gamma)', 'E > C : nu * gamma', 'A > R : tauA', 'C > R : tauC * (1 - theta)',
      'C > H : tauC * theta | out=admissions', 'H > R : tauH * (1 - omega)', 'H > D : tauH * omega | out=deaths', 'R > ER : sigma * lam',
      'ER > A : nu * (1 - gamma / 2)', 'ER > C : nu * gamma / 2', 'R > S : alpha | curve=-190',
    ],
  }),
  def({
    id: 'seeiir', category: 'Natural history', name: 'SEEIIR (gamma-distributed periods)',
    description: 'Two latent and two infectious stages give gamma-distributed (Erlang k = 2) latent and infectious\nperiods with the same means as the SEIR model ("linear chain trick").',
    reference: KR + ', §3.3',
    time: [0, 200, 1, 'days'],
    comps: [
      ['S', 'Susceptible', 9999, 0, 0], ['E1', 'Latent, stage 1', 1, 1, 0], ['E2', 'Latent, stage 2', 0, 2, 0],
      ['I1', 'Infectious, stage 1', 0, 3, 0], ['I2', 'Infectious, stage 2', 0, 4, 0], ['R', 'Recovered', 0, 5, 0],
    ],
    params: [
      ['beta', 0.5, 'transmission rate per day'],
      ['sigma', '1 / 5', '1/(mean latent period)'],
      ['gamma', '1 / 7', '1/(mean infectious period)'],
    ],
    vars: [
      ['lam', 'beta * (I1 + I2) / N', 'force of infection'],
      ['prevalence', 'I1 + I2', 'number infectious', true],
    ],
    flows: ['S > E1 : lam', 'E1 > E2 : 2 * sigma', 'E2 > I1 : 2 * sigma | out=onsets', 'I1 > I2 : 2 * gamma', 'I2 > R : 2 * gamma'],
  }),
  def({
    id: 'sirws', category: 'Natural history', name: 'SIRWS (waning and boosting)',
    description: 'Immunity wanes through a W stage; re-exposure while in W boosts immunity back to R (pertussis-like).',
    reference: 'Lavine, King & Bjornstad (2011), PNAS 108:7259-7264',
    time: [0, 100, 0.02, 'years'],
    comps: [['S', 'Susceptible', 5000, 0, 0], ['I', 'Infectious', 10, 1, 0], ['R', 'Immune', 70000, 2, 0], ['W', 'Waning immunity', 24990, 1, 1]],
    params: [
      ['beta', 260, 'transmission rate per year'],
      ['gamma', '365 / 21', 'recovery rate per year (21-day infectious period)'],
      ['kappa', '1 / 10', 'rate of waning per year (two stages, so mean 2 / (2 kappa) = 10 years)'],
      ['boost', 20, 'boosting factor relative to the force of infection'],
      ['mu', '1 / 70', 'birth = death rate per year'],
    ],
    vars: [['lam', 'beta * I / N', 'force of infection']],
    flows: [
      ' > S : mu * N', 'S > I : lam | out=incidence', 'I > R : gamma', 'R > W : 2 * kappa | curve=30', 'W > S : 2 * kappa',
      'W > R : boost * lam | curve=30', 'S > : mu | at=-60,110', 'I > : mu | at=0,-130', 'R > : mu', 'W > : mu',
    ],
  }),
  def({
    id: 'carrier', category: 'Natural history', name: 'SICR (chronic carriers)',
    description: 'A proportion q of acute infections become long-term carriers (hepatitis B / typhoid-like).',
    reference: AM + ', ch. 6; ' + KR + ', §3.2',
    time: [0, 150, 0.1, 'years'],
    comps: [['S', 'Susceptible', 99000, 0, 0.5], ['I', 'Acute infection', 10, 1, 0.5], ['C', 'Chronic carrier', 990, 2, 1], ['R', 'Recovered', 0, 2, 0]],
    params: [
      ['beta', 4, 'transmission rate per year'],
      ['eps', 0.1, 'relative infectiousness of carriers', 'ε'],
      ['gamma', 6, 'rate of leaving acute infection per year (2 months)'],
      ['q', 0.1, 'proportion becoming carriers'],
      ['gammaC', 0.025, 'rate of carriers clearing infection per year', 'γ_C'],
      ['mu', '1 / 70', 'birth = death rate per year'],
    ],
    vars: [['lam', 'beta * (I + eps * C) / N', 'force of infection']],
    flows: [' > S : mu * N', 'S > I : lam | out=incidence', 'I > R : gamma * (1 - q)', 'I > C : gamma * q', 'C > R : gammaC', 'S > : mu', 'I > : mu', 'C > : mu', 'R > : mu | at=0,-130'],
  }),
  def({
    id: 'tb', category: 'Natural history', name: 'Tuberculosis (fast and slow progression)',
    description: 'New infections progress rapidly to disease with probability p, otherwise enter long-term latency (L)\nwith slow reactivation; treated/self-cured cases can relapse.',
    reference: VW + ', ch. 3; Blower et al. (1995), Nature Medicine 1:815-821',
    time: [0, 200, 0.1, 'years'],
    comps: [['S', 'Susceptible', 99990, 0, 0.5], ['L', 'Latent infection', 0, 1, 1], ['I', 'Active TB', 10, 2, 0.5], ['R', 'Recovered', 0, 3, 0.5]],
    params: [
      ['beta', 10, 'effective contacts per case per year'],
      ['p', 0.1, 'proportion progressing rapidly to disease'],
      ['v', 0.001, 'reactivation rate per year'],
      ['r', 0.5, 'rate of cure / self-cure per year'],
      ['relapse', 0.002, 'relapse rate per year'],
      ['mu', '1 / 70', 'background death rate per year'],
      ['muTB', 0.2, 'TB mortality rate per year', 'μ_TB'],
    ],
    vars: [['lam', 'beta * I / N', 'force of infection']],
    flows: [
      ' > S : mu * N + muTB * I', 'S > I : p * lam | out=fast_cases', 'S > L : (1 - p) * lam', 'L > I : v', 'I > R : r', 'R > I : relapse | curve=40',
      'S > : mu', 'L > : mu', 'I > : mu | at=-30,-130', 'I > : muTB | at=30,-130', 'R > : mu',
    ],
  }),
  def({
    id: 'hiv_stages', category: 'Natural history', name: 'HIV with disease stages',
    description: 'Acute (I1), chronic (I2) and AIDS (A) stages with stage-specific infectiousness; AIDS deaths leave the population.',
    reference: AM + ', ch. 11',
    time: [0, 50, 0.05, 'years'],
    comps: [['S', 'Susceptible', 99900, 0, 0], ['I1', 'Acute HIV', 100, 1, 0], ['I2', 'Chronic HIV', 0, 2, 0], ['A', 'AIDS', 0, 3, 0], ['D', 'AIDS deaths', 0, 3, 1, { inN: false }]],
    params: [
      ['cp', 2, 'partner change rate per year', 'c'],
      ['b1', 0.1, 'transmission probability per partnership, acute', 'β_1'],
      ['b2', 0.02, 'transmission probability per partnership, chronic', 'β_2'],
      ['b3', 0.05, 'transmission probability per partnership, AIDS', 'β_3'],
      ['k1', 4, 'rate of leaving acute stage per year (3 months)', 'κ_1'],
      ['k2', '1 / 9', 'rate of progression to AIDS per year', 'κ_2'],
      ['alpha', 0.5, 'AIDS death rate per year'],
      ['mu', '1 / 35', 'rate of entering/leaving the sexually active population'],
    ],
    vars: [['lam', 'cp * (b1 * I1 + b2 * I2 + b3 * A) / N', 'force of infection']],
    flows: [' > S : mu * N', 'S > I1 : lam | out=incidence', 'I1 > I2 : k1', 'I2 > A : k2', 'A > D : alpha | out=aids_deaths', 'S > : mu', 'I1 > : mu', 'I2 > : mu', 'A > : mu | at=0,-130'],
  }),
  def({
    id: 'ebola', category: 'Natural history', name: 'Ebola (community, hospital, funeral)',
    description: 'Transmission in the community (I), in hospital (H) and at funerals (F); structure after Legrand et al.\nParameter values are illustrative.',
    reference: 'Legrand et al. (2007), Epidemiology and Infection 135:610-621',
    time: [0, 300, 0.5, 'days'],
    comps: [
      ['S', 'Susceptible', 99999, 0, 0.5], ['E', 'Exposed', 1, 1, 0.5], ['I', 'Infectious (community)', 0, 2, 0.5],
      ['H', 'Hospitalised', 0, 3, 1.3], ['Fu', 'Dead, not yet buried (funeral)', 0, 4, 0.5, { color: '#5D4037' }], ['R', 'Removed', 0, 3, -0.3],
    ],
    params: [
      ['betaI', 0.3, 'transmission rate in the community', 'β_I'],
      ['betaH', 0.05, 'transmission rate in hospital', 'β_H'],
      ['betaF', 0.4, 'transmission rate at funerals', 'β_F'],
      ['alpha', '1 / 7', '1/(incubation period)'],
      ['theta', 0.5, 'proportion of cases hospitalised'],
      ['delta', 0.7, 'case fatality ratio'],
      ['gammaH', '1 / 5', '1/(onset to hospitalisation)', 'γ_h'],
      ['gammaD', '1 / 9.6', '1/(onset to death)', 'γ_d'],
      ['gammaI', '1 / 10', '1/(onset to recovery)', 'γ_i'],
      ['gammaDH', '1 / 4.6', '1/(hospitalisation to death)', 'γ_dh'],
      ['gammaIH', '1 / 5', '1/(hospitalisation to recovery)', 'γ_ih'],
      ['gammaF', '1 / 2', '1/(death to burial)', 'γ_f'],
    ],
    vars: [['lam', '(betaI * I + betaH * H + betaF * Fu) / N', 'force of infection']],
    flows: [
      'S > E : lam', 'E > I : alpha | out=onsets', 'I > H : gammaH * theta', 'I > Fu : gammaD * (1 - theta) * delta | curve=-40',
      'I > R : gammaI * (1 - theta) * (1 - delta)', 'H > Fu : gammaDH * delta', 'H > R : gammaIH * (1 - delta)', 'Fu > : gammaF | out=burials',
    ],
  }),

  // -------------------------------------------------------------------------
  // Transmission routes & structure
  // -------------------------------------------------------------------------
  def({
    id: 'cholera_siwr', category: 'Transmission routes', name: 'Cholera SIWR (water reservoir)',
    description: 'Direct transmission plus indirect transmission through a contaminated water compartment W\n(W is a normalised pathogen concentration, so it is excluded from N).',
    reference: 'Tien & Earn (2010), Bull. Math. Biol. 72:1506-1533',
    time: [0, 200, 0.5, 'days'],
    comps: [['S', 'Susceptible', 9990, 0, 0], ['I', 'Infectious', 10, 1, 0], ['R', 'Recovered', 0, 2, 0], ['W', 'Water (pathogen)', 0, 1, 1, { inN: false, color: '#4FC3F7' }]],
    params: [
      ['betaI', 0.1, 'direct transmission rate per day', 'β_I'],
      ['betaW', 0.4, 'water-borne transmission rate per day', 'β_W'],
      ['gamma', '1 / 5', 'recovery rate'],
      ['xi', '1 / 30', 'shedding / decay rate of pathogen in water', 'ξ'],
    ],
    vars: [['lam', 'betaI * I / N + betaW * W', 'force of infection']],
    flows: ['S > I : lam | out=incidence', 'I > R : gamma', ' > W : xi * I / N | at=-160,0', 'W > : xi'],
  }),
  def({
    id: 'ross_macdonald', category: 'Transmission routes', name: 'Ross–Macdonald (malaria)',
    description: 'SIS humans and SI mosquitoes. Two separate populations, so totals are defined as variables.',
    reference: 'Smith et al. (2012), PLoS Pathogens 8:e1002588; ' + AM + ', ch. 14',
    time: [0, 730, 1, 'days'],
    comps: [
      ['Sh', 'Susceptible humans', 9990, 0, 0, { inN: false }], ['Ih', 'Infected humans', 10, 1, 0, { inN: false }],
      ['Sv', 'Susceptible mosquitoes', 50000, 0, 1.4, { inN: false }], ['Iv', 'Infectious mosquitoes', 0, 1, 1.4, { inN: false }],
    ],
    params: [
      ['a', 0.3, 'human bites per mosquito per day'],
      ['b', 0.5, 'probability an infectious bite infects a human'],
      ['cv', 0.5, 'probability a mosquito is infected biting an infected human', 'c'],
      ['r', '1 / 100', 'human recovery rate per day'],
      ['g', '1 / 10', 'mosquito death rate per day'],
    ],
    vars: [
      ['Nh', 'Sh + Ih', 'human population'],
      ['Nv', 'Sv + Iv', 'mosquito population'],
      ['lam_h', 'a * b * Iv / Nh', 'force of infection on humans', false, 'λ_h'],
      ['lam_v', 'a * cv * Ih / Nh', 'force of infection on mosquitoes', false, 'λ_v'],
    ],
    flows: ['Sh > Ih : lam_h | curve=30 out=human_incidence', 'Ih > Sh : r | curve=30', ' > Sv : g * Nv | at=-150,0', 'Sv > Iv : lam_v', 'Sv > : g', 'Iv > : g'],
  }),
  def({
    id: 'dengue_vector', category: 'Transmission routes', name: 'Dengue host SEIR – vector SEI',
    description: 'Human SEIR coupled to mosquito SEI with an extrinsic incubation period.',
    reference: 'Manore et al. (2014), J. Theor. Biol. 356:174-191 (structure)',
    time: [0, 365, 0.5, 'days'],
    comps: [
      ['Sh', 'Susceptible humans', 99990, 0, 0, { inN: false }], ['Eh', 'Exposed humans', 0, 1, 0, { inN: false }],
      ['Ih', 'Infectious humans', 10, 2, 0, { inN: false }], ['Rh', 'Recovered humans', 0, 3, 0, { inN: false }],
      ['Sv', 'Susceptible mosquitoes', 200000, 0, 1.5, { inN: false }], ['Ev', 'Exposed mosquitoes', 0, 1, 1.5, { inN: false }],
      ['Iv', 'Infectious mosquitoes', 0, 2, 1.5, { inN: false }],
    ],
    params: [
      ['a', 0.5, 'bites per mosquito per day'],
      ['bh', 0.4, 'mosquito-to-human transmission probability', 'b_h'],
      ['bv', 0.4, 'human-to-mosquito transmission probability', 'b_v'],
      ['sigmaH', '1 / 5.5', '1/(intrinsic incubation period)', 'σ_h'],
      ['gammaH', '1 / 5', 'human recovery rate', 'γ_h'],
      ['sigmaV', '1 / 10', '1/(extrinsic incubation period)', 'σ_v'],
      ['muV', '1 / 14', 'mosquito death rate', 'μ_v'],
    ],
    vars: [
      ['Nh', 'Sh + Eh + Ih + Rh', 'human population'],
      ['Nv', 'Sv + Ev + Iv', 'mosquito population'],
      ['lam_h', 'a * bh * Iv / Nh', 'force of infection on humans', false, 'λ_h'],
      ['lam_v', 'a * bv * Ih / Nh', 'force of infection on mosquitoes', false, 'λ_v'],
    ],
    flows: [
      'Sh > Eh : lam_h', 'Eh > Ih : sigmaH | out=human_cases', 'Ih > Rh : gammaH', ' > Sv : muV * Nv | at=-150,0',
      'Sv > Ev : lam_v', 'Ev > Iv : sigmaV', 'Sv > : muV', 'Ev > : muV', 'Iv > : muV',
    ],
  }),
  def({
    id: 'rabies_density', category: 'Transmission routes', name: 'Rabies SEI (density-dependent)',
    description: 'Wildlife rabies with density-dependent transmission (lam = beta * I, not divided by N)\nand logistic population regulation through density-dependent mortality.',
    reference: 'Anderson et al. (1981), Nature 289:765-771',
    time: [0, 30, 0.01, 'years'],
    comps: [['S', 'Susceptible foxes', 900, 0, 0], ['E', 'Incubating', 0, 1, 0], ['I', 'Rabid', 5, 2, 0]],
    params: [
      ['b', 1, 'per-capita birth rate per year'],
      ['mu', 0.5, 'baseline death rate per year'],
      ['K', 1000, 'carrying capacity'],
      ['beta', 0.08, 'density-dependent transmission coefficient per year'],
      ['sigma', 13, '1/(incubation period) per year (~4 weeks)'],
      ['alpha', 73, 'rabies death rate per year (~5 days)'],
    ],
    vars: [
      ['d', 'mu + (b - mu) * N / K', 'density-dependent death rate'],
      ['lam', 'beta * I', 'force of infection (density dependent)'],
    ],
    flows: [' > S : b * N', 'S > E : lam', 'E > I : sigma | out=cases', 'S > : d', 'E > : d', 'I > : d | at=-30,130', 'I > : alpha | at=40,130 out=rabies_deaths'],
  }),
  def({
    id: 'sir_two_age', category: 'Transmission routes', name: 'SIR with two age groups',
    description: 'Children and adults with a who-acquires-infection-from-whom (WAIFW) matrix and ageing.',
    reference: AM + ', ch. 9; ' + KR + ', §3.1',
    time: [0, 50, 0.02, 'years'],
    comps: [
      ['Sc', 'Susceptible children', 2000, 0, 0], ['Ic', 'Infectious children', 10, 1, 0], ['Rc', 'Recovered children', 17990, 2, 0],
      ['Sa', 'Susceptible adults', 1000, 0, 1.2], ['Ia', 'Infectious adults', 0, 1, 1.2], ['Ra', 'Recovered adults', 79000, 2, 1.2],
    ],
    params: [
      ['b_cc', 1000, 'transmission child -> child per year', 'β_cc'],
      ['b_ca', 300, 'transmission adult -> child per year', 'β_ca'],
      ['b_ac', 300, 'transmission child -> adult per year', 'β_ac'],
      ['b_aa', 400, 'transmission adult -> adult per year', 'β_aa'],
      ['gamma', '365 / 7', 'recovery rate per year'],
      ['ageing', '1 / 15', 'rate of moving from child to adult per year'],
      ['mu', '1 / 70', 'birth rate per year (deaths in adults keep N constant)'],
    ],
    vars: [
      ['lam_c', '(b_cc * Ic + b_ca * Ia) / N', 'force of infection on children', false, 'λ_c'],
      ['lam_a', '(b_ac * Ic + b_aa * Ia) / N', 'force of infection on adults', false, 'λ_a'],
      ['mu_a', 'mu * N / (Sa + Ia + Ra)', 'adult death rate balancing births'],
    ],
    flows: [
      ' > Sc : mu * N', 'Sc > Ic : lam_c', 'Ic > Rc : gamma', 'Sa > Ia : lam_a', 'Ia > Ra : gamma',
      'Sc > Sa : ageing', 'Ic > Ia : ageing', 'Rc > Ra : ageing', 'Sa > : mu_a', 'Ia > : mu_a', 'Ra > : mu_a',
    ],
  }),
  def({
    id: 'seir_seasonal', category: 'Transmission routes', name: 'SEIR with seasonal forcing',
    description: 'Sinusoidal seasonality in transmission produces multi-annual epidemic cycles (measles-like).',
    reference: KR + ', §5.2',
    time: [0, 100, 0.01, 'years'],
    comps: [['S', 'Susceptible', 6000, 0, 0], ['E', 'Exposed', 10, 1, 0], ['I', 'Infectious', 10, 2, 0], ['R', 'Recovered', 93980, 3, 0]],
    params: [
      ['beta0', 1250, 'mean transmission rate per year', 'β_0'],
      ['beta1', 0.25, 'amplitude of seasonality', 'β_1'],
      ['sigma', '365 / 8', '1/(latent period) per year'],
      ['gamma', '365 / 5', '1/(infectious period) per year'],
      ['mu', '1 / 50', 'birth = death rate per year'],
    ],
    vars: [
      ['beta', 'beta0 * (1 + beta1 * cos(2 * pi * t))', 'seasonal transmission rate'],
      ['lam', 'beta * I / N', 'force of infection'],
    ],
    flows: [' > S : mu * N', 'S > E : lam', 'E > I : sigma | out=incidence', 'I > R : gamma', 'S > : mu', 'E > : mu', 'I > : mu', 'R > : mu'],
  }),
];

// Models imported from the research notes (see cli/import-research.mjs). A model is only listed
// if it builds and passes validation; QUARANTINE holds ids that validate but failed to solve in R.
// REJECTED records every model left out and why (printed by `npm run examples`).
export const QUARANTINE = new Set([]);
export const REJECTED = [];
const known = new Set(LIBRARY.map((m) => m.id));
for (const spec of RESEARCH) {
  if (known.has(spec.id)) continue;
  if (QUARANTINE.has(spec.id)) { REJECTED.push({ id: spec.id, reasons: ['failed to solve in R (see QUARANTINE)'] }); continue; }
  let m;
  try { m = def(spec); } catch (e) { REJECTED.push({ id: spec.id, reasons: [e.message] }); continue; }
  const errors = validate(m).errors;
  if (errors.length) { REJECTED.push({ id: spec.id, reasons: [...new Set(errors.map((x) => x.message))] }); continue; }
  LIBRARY.push(m);
  known.add(m.id);
}

export const BY_ID = Object.fromEntries(LIBRARY.map((m) => [m.id, m]));

export const CATEGORIES = [...new Set(LIBRARY.map((m) => m.category))];
