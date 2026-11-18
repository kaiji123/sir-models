#!/usr/bin/env node
// Converts the structured research notes (research_notes/<topic>/*.md, one "### id — name"
// block per model) into src/lib/library-research.js, which library.js appends to the library.
//
// Usage: node cli/import-research.mjs [--skip id1,id2] [notes-dir]
// Prints a report of every model: imported, or rejected with the validator's reasons.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { analyze } from '../src/lib/generator.js';

const root = path.join(path.dirname(fileURLToPath(import.meta.url)), '..');
const args = process.argv.slice(2);
const skipIdx = args.indexOf('--skip');
// Duplicates found in the literature review: keep the more detailed version of each pair.
const DUPLICATES = [
  'hiv_test_and_treat',      // simplified Granich 2009; hiv_art_tasp_granich has the 4-stage structure
  'sis_core_group_gonorrhea', // same as gonorrhoea_core_group_hethcote_yorke (two sexes)
  'flu_two_strain_cross',    // same structure as two_strain_cross_immunity
];
const skip = new Set([...DUPLICATES, ...(skipIdx >= 0 ? args.splice(skipIdx, 2)[1].split(',') : [])]);
const notesDir = args[0] || path.join(root, 'research_notes', 'Compartmental model library catalogue');
const outFile = path.join(root, 'src', 'lib', 'library-research.js');

// Category per notes file, used when a model's own category line is missing or unknown.
const FILE_CATEGORY = {
  'classics_interventions.md': 'Classic extensions',
  'vector_borne.md': 'Vector-borne',
  'respiratory.md': 'Respiratory',
  'sti_bloodborne.md': 'Sexually transmitted & blood-borne',
  'enteric_ntd.md': 'Water, enteric & NTDs',
  'zoonoses_animal_withinhost.md': 'Zoonoses & emerging',
};
const KNOWN_CATEGORIES = ['Vector-borne', 'Respiratory', 'Sexually transmitted & blood-borne', 'Water, enteric & NTDs',
  'Zoonoses & emerging', 'Animal & wildlife', 'Within-host'];

const stripMd = (s) => String(s || '').replace(/\[([^\]]*)\]\([^)]*\)/g, '$1').replace(/[*`]/g, '').replace(/\s+/g, ' ').trim();
const dashSplit = (s) => s.split(/\s+[—–]\s+/);

function firstNumber(s) {
  s = String(s).replace(/(\d),(?=\d{3}b)/g, '$1');
  const m = /-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?(?:\s*\/\s*\d+(?:\.\d+)?)?/.exec(String(s));
  return m ? m[0].replace(/\s+/g, ' ') : null;
}

function exprOrNumber(s) {
  const t = String(s).trim().replace(/`/g, '').replace(/(\d),(?=\d{3}b)/g, '$1');
  if (analyze(t).ok) return t;
  const head = t.split(/\s+\(|\s+[—–]|;|,\s/)[0].trim();
  if (head && analyze(head).ok) return head;
  return firstNumber(t);
}

/** Split the notes file into "### id — name" blocks. */
function blocks(md) {
  const out = [];
  const parts = md.split(/^#{3,4} /m).slice(1);
  for (const part of parts) {
    const nl = part.indexOf('\n');
    const header = part.slice(0, nl).trim();
    const m = /^`?([a-z][a-z0-9_]*)`?\s*[—–-]+\s*(.+)$/.exec(header);
    if (!m) continue;
    out.push({ id: m[1], name: stripMd(m[2]), body: part.slice(nl + 1) });
  }
  return out;
}

/** Collect "- Key: ..." sections with their nested "  - item" lines. */
function sections(body) {
  const secs = {};
  let cur = null;
  for (const line of body.split(/\r?\n/)) {
    const top = /^- ([A-Za-z][A-Za-z /&-]*?)(\s*\([^)]*\))?:\s*(.*)$/.exec(line);
    if (top) { cur = top[1].trim().toLowerCase(); secs[cur] = { head: top[3], items: [] }; continue; }
    const item = /^\s{2,}[-*] (.*)$/.exec(line);
    if (item && cur) { secs[cur].items.push(item[1]); continue; }
    if (/^#{1,4} /.test(line)) cur = null;
  }
  // Inline lists ("- Flows: `A > B : x`; `B > C : y`") become items too.
  for (const s of Object.values(secs)) {
    if (!s.items.length && s.head.includes('`')) {
      const ticks = s.head.match(/`[^`]+`[^`;]*/g);
      if (ticks && ticks.length > 1) s.items = ticks;
    }
  }
  return secs;
}

function parseModel(b, file) {
  const s = sections(b.body);
  const problems = [];
  const get = (k) => Object.keys(s).find((x) => x.startsWith(k));

  // --- category, disease
  const dis = s[get('disease')]?.head || '';
  const catMatch = /Category:\s*([^|]+)/.exec(dis);
  let category = catMatch ? stripMd(catMatch[1]) : '';
  if (!KNOWN_CATEGORIES.includes(category)) category = FILE_CATEGORY[file] || 'Other';
  const disease = stripMd(dis.split('|')[0]);

  // --- source
  // The citation is on the "- Source:" line itself, or (when that line is empty)
  // in the first nested item. For a nested item keep just the citation: drop a
  // leading label ("Model:") or explanation ("... See Smith DL ..."), and
  // everything from its first link onwards (the link itself goes into url).
  const srcSec = s[get('source')];
  let srcRaw = (srcSec?.head || '').trim();
  let linkSrc = srcRaw;
  if (!srcRaw && srcSec?.items.length) {
    linkSrc = srcSec.items[0].replace(/^(?:model|structural basis|structure motivated by)\s*:?\s*/i, '');
    srcRaw = linkSrc.replace(/^.*?\bSee (?=[A-Z])/, '').split(/\s*\(?(?:https?:\/\/|PubMed:)/)[0];
  }
  // DOIs may contain balanced brackets and ";" inside them, e.g.
  // 10.1016/S0022-5193(03)00228-5 or 10.1890/1051-0761(2006)016[2208:DOPDTI]2.0.CO;2.
  const doi = /\b(10\.\d{4,9}\/(?:[^\s;,()[\]]|\([^\s()]*\)|\[[^\s[\]]*\]|;(?=\w))+)/.exec(linkSrc);
  const url = doi ? 'https://doi.org/' + doi[1].replace(/[.]$/, '') : (/(https?:\/\/[^\s)\]]+)/.exec(linkSrc) || [])[1] || '';
  let reference = stripMd(srcRaw.replace(/\(\[[^\]]*\]\([^)]*\)[^)]*\)/g, '').split(/\s*DOI[:\s]/i)[0]).replace(/[.;,\s]+$/, '');
  if (reference.length > 260) reference = reference.slice(0, 257).replace(/\s+\S*$/, '') + '…';

  // --- time
  const timeHead = (s[get('time')]?.head || 'days').toLowerCase();
  const unit = (/(day|week|month|year|hour)s?/.exec(timeHead) || ['days', 'day'])[1] + 's';
  const lenMatch = /(\d+(?:\.\d+)?)\s*(?:day|week|month|year|hour)/.exec(timeHead.split(';').slice(1).join(';')) || /(\d+(?:\.\d+)?)/.exec(timeHead.split(';').slice(1).join(';'));
  const stop = lenMatch ? +lenMatch[1] : { days: 365, weeks: 104, years: 50, months: 60, hours: 240 }[unit] || 365;
  const raw = stop / 1000;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 5, 10].map((k) => k * mag).find((v) => v >= raw * 0.999);
  const stepStr = step >= 1 ? String(step) : String(+step.toPrecision(3));

  // --- compartments ("`A` — desc — init" or grouped "`A`, `B` — desc — initA, initB")
  const comps = [];
  for (const it of s[get('compartment')]?.items || []) {
    const m = /^((?:`[A-Za-z][A-Za-z0-9_.]*`\s*(?:,|and|\/)?\s*)+)[—–-]+\s*(.*)$/.exec(it.trim());
    if (!m) continue;
    const names = [...m[1].matchAll(/`([A-Za-z][A-Za-z0-9_.]*)`/g)].map((x) => x[1]);
    const segs = dashSplit(m[2]);
    const notInN = /not in N/i.test(it);
    let initSeg = null;
    for (let i = segs.length - 1; i >= 1 && initSeg == null; i--) if (!/not in N/i.test(segs[i]) && firstNumber(segs[i]) != null) initSeg = segs[i];
    const inits = names.length > 1 && initSeg ? initSeg.replace(/(\d),(?=\d{3}b)/g, '$1').split(/\s*,\s*/) : [initSeg];
    names.forEach((name, k) => {
      const raw = inits[k] ?? inits[inits.length - 1];
      comps.push({ name, description: stripMd(segs[0]).slice(0, 120), init: (raw != null && exprOrNumber(raw)) || '0', inN: !notInN });
    });
  }

  // --- parameters ("`a` = v (unit) — meaning — provenance", also "`a` = `b` = v" and "`a` = v1, `b` = v2")
  const params = [];
  let verified = 0;
  const addParam = (name, valueText, segs) => {
    if (params.some((p) => p.name === name)) return;
    const value = exprOrNumber(valueText.replace(/\s*\([^)]*\)\s*$/, '')) ?? exprOrNumber(valueText);
    if (value == null) { problems.push(`parameter ${name}: cannot read value "${valueText}"`); return; }
    const unitM = /\(([^)]*)\)\s*$/.exec(valueText);
    const prov = segs.slice(2).join(' — ');
    const isVerified = /\[(?:verified|v)\b/i.test(prov) || /^\s*\[?verified/i.test(prov);
    if (isVerified) verified++;
    const tag = isVerified ? 'verified' : /illustrative|\[i\]/i.test(prov) ? 'illustrative' : /recalled|memory|\[m\]|unverified/i.test(prov) ? 'recalled' : '';
    const desc = stripMd(segs[1] || '') + (unitM && unitM[1] && !/^[–-]$/.test(unitM[1]) ? ` (${stripMd(unitM[1])})` : '') + (tag ? ` [${tag}]` : '');
    params.push({ name, value, description: desc.slice(0, 160) });
  };
  for (const it of s[get('parameter')]?.items || []) {
    if (!/^`[A-Za-z]/.test(it.trim())) continue;
    const segs = dashSplit(it.trim());
    // list form: "`a`, `b`, `c` = 1, 2, 3 (unit)"
    const list = /^((?:`[A-Za-z][A-Za-z0-9_.]*`\s*,\s*)+`[A-Za-z][A-Za-z0-9_.]*`)\s*=\s*(.+)$/.exec(segs[0]);
    if (list) {
      const names = [...list[1].matchAll(/`([A-Za-z][A-Za-z0-9_.]*)`/g)].map((x) => x[1]);
      const values = list[2].replace(/\s*\([^)]*\)\s*$/, '').split(/\s*,\s*/);
      names.forEach((n, k) => addParam(n, values[k] ?? values[values.length - 1], segs));
      continue;
    }
    for (const piece of segs[0].split(/,\s*(?=`)/)) {
      const names = [...piece.matchAll(/`([A-Za-z][A-Za-z0-9_.]*)`\s*=/g)].map((x) => x[1]);
      const valueText = piece.slice(piece.lastIndexOf('=') + 1).trim();
      names.forEach((n) => addParam(n, valueText, segs));
    }
  }

  // --- variables
  const vars = [];
  for (const it of s[get('derived')]?.items || []) {
    const m = /^`([A-Za-z][A-Za-z0-9_.]*)`\s*=\s*(.+)$/.exec(it.trim());
    if (!m) continue;
    const segs = dashSplit(m[2]);
    let expr = segs[0].replace(/`/g, '').trim();
    // drop a trailing prose comment such as "(used for both beta and delta)"
    while (!analyze(expr).ok && /\s*\([^()]*[A-Za-z]{3,}\s+[A-Za-z][^()]*\)\s*$/.test(expr)) expr = expr.replace(/\s*\([^()]*\)\s*$/, '');
    vars.push({ name: m[1], expr, description: stripMd(segs.slice(1).join(' — ')).slice(0, 140) });
  }

  // --- flows ("`A > B : rate`", optionally "(absolute)"; "each of `A`, `B` `> : rate`" expands)
  const flows = [];
  for (const it of s[get('flow')]?.items || []) {
    const m = /`\s*([A-Za-z0-9_.]*)\s*>\s*([A-Za-z0-9_.]*)\s*:\s*([^`]+)`(.*)$/.exec(it);
    if (!m) continue;
    const absolute = (!m[1] && !!m[2]) || /\babsolute\b/i.test(m[3] + m[4]);
    const rate = m[3].replace(/\s*[([][^()[\]]*\babsolute\b[^()[\]]*[)\]]/gi, '').replace(/\s*\babsolute\b.*$/i, '').trim();
    if (!m[1] && !m[2]) {
      const before = it.slice(0, m.index);
      const srcs = [...before.matchAll(/`([A-Za-z][A-Za-z0-9_.]*)`/g)].map((x) => x[1]);
      srcs.forEach((src) => flows.push({ from: src, to: null, rate, absolute: false }));
      continue;
    }
    flows.push({ from: m[1] || null, to: m[2] || null, rate, absolute });
  }

  // --- names defined only in the prose ("x = 5" in notes) become parameters
  const declared = new Set([...comps.map((c) => c.name), ...params.map((p) => p.name), ...vars.map((v) => v.name), 't', 'pi', 'N']);
  const used = new Set();
  [...comps.map((c) => c.init), ...vars.map((v) => v.expr), ...flows.map((f) => f.rate)].forEach((e) => (analyze(e).idents || []).forEach((x) => used.add(x)));
  for (const name of used) {
    if (declared.has(name)) continue;
    const esc = name.replace(/[.]/g, '\.');
    const m = new RegExp(String.raw`\b${esc}\b` + String.raw`` + '`?' + String.raw`\s*=\s*` + '`?' + String.raw`(-?\d[\d,]*(?:\.\d+)?(?:[eE][+-]?\d+)?)`).exec(b.body);
    if (m) params.push({ name, value: m[1].replace(/,/g, ''), description: 'value taken from the notes text [recalled]' });
    else if (/0$/.test(name) && comps.some((c) => c.name === name.slice(0, -1) && !analyze(c.init).idents.length)) {
      const c = comps.find((x) => x.name === name.slice(0, -1));
      params.push({ name, value: c.init, description: `initial value of ${c.name}` });
    } else if (name === 'N0') params.push({ name, value: '100000', description: 'total initial population [illustrative default]' });
    else problems.push(`undefined name ${name}`);
  }

  const notes = stripMd(s[get('notes')]?.head || '');
  const firstSentence = (notes.match(/^.{20,240}?[.!?](\s|$)/) || [notes.slice(0, 240)])[0].trim();
  const provenance = params.length
    ? `Parameter values: ${verified} of ${params.length} checked against the source; the rest are recalled or illustrative.`
    : '';
  const description = [disease ? disease.replace(/\.$/, '') + '.' : '', firstSentence, provenance].filter(Boolean).join('\n');
  const keywords = disease.split(/[,;()/]| and /).map((w) => w.trim().toLowerCase()).filter((w) => w && w.length < 40);

  const ownN = [...comps, ...params, ...vars].some((x) => x.name === 'N');
  return {
    spec: { id: b.id, pop: ownN ? 'N_total' : undefined, category, name: b.name, description, reference, url, keywords, time: [0, stop, stepStr, unit], comps, params, vars, flows },
    problems,
  };
}

// Corrections applied after parsing, so they survive re-imports. Keep each one justified.
const OVERRIDES = {
  // Illustrative values gave one huge epidemic and a ~1e-22 trough for decades, where solver noise
  // turns I negative. Gentler illustrative values for a chronic, vertically transmitted infection.
  sir_vertical_transmission: { params: { beta: '1.5', gamma: '0.5', b_birth: '1 / 20', mu: '1 / 20' } },
};

const COUNTER_RE = /cumulat|counter|bookkeep|accumulat|tally/i;

/**
 * A "counter" compartment (cumulative cases, doses...) fed from a living compartment that already
 * has its real outflow would drain that compartment twice. Replace such flows by a recorded output
 * (the rate), and drop the counter compartment if nothing else uses it. Returns notes for the report.
 */
function convertCounters(spec) {
  const notes = [];
  const byName = Object.fromEntries(spec.comps.map((c) => [c.name, c]));
  const keep = [];
  for (const f of spec.flows) {
    const tgt = f.to && byName[f.to];
    const src = f.from && byName[f.from];
    // Only a counter if its rate duplicates (all of, or one summed term of) another flow from the
    // same source: e.g. E > CC : k next to E > I : k, or S > CZ : lam_z next to S > E : lam_z + lam_s.
    // A real transfer such as I > D : gamma * cfr (next to I > R : gamma * (1 - cfr)) is kept.
    const canon = (e) => analyze(e).code || e;
    const terms = (e) => { const c = canon(e).replace(/^\((.*)\)$/, '$1'); return c.split(/ \+ /).map(canon); };
    const isCounter = tgt && src && !tgt.inN && src.inN && COUNTER_RE.test(tgt.description)
      && spec.flows.some((g) => g !== f && g.from === f.from && g.to && byName[g.to] && byName[g.to].inN
        && (canon(g.rate) === canon(f.rate) || terms(g.rate).includes(canon(f.rate))));
    if (!isCounter) { keep.push(f); continue; }
    const usesSrc = (analyze(f.rate).idents || []).includes(src.name);
    const term = usesSrc ? f.rate : `(${f.rate}) * ${src.name}`;
    const name = `${tgt.name}_rate`;
    spec.vars.push({ name, expr: term, description: `rate of new ${tgt.description} (integrate with cumsum for the cumulative count)`, output: true });
    notes.push(`${f.from} > ${f.to} counter replaced by output ${name}`);
  }
  spec.flows = keep;
  spec.comps = spec.comps.filter((c) => {
    if (!COUNTER_RE.test(c.description) || c.inN) return true;
    const used = spec.flows.some((f) => f.from === c.name || f.to === c.name)
      || [...spec.vars.map((v) => v.expr), ...spec.flows.map((f) => f.rate), ...spec.comps.map((x) => x.init)].some((e) => (analyze(e).idents || []).includes(c.name));
    if (!used) notes.push(`removed unused counter ${c.name}`);
    return used;
  });
  return notes;
}

function applyOverrides(spec) {
  const o = OVERRIDES[spec.id];
  if (!o) return [];
  for (const [name, value] of Object.entries(o.params || {})) {
    const p = spec.params.find((x) => x.name === name);
    if (p) { p.value = value; p.description = p.description.replace(/\s*\[[a-z]+\]$/, '') + ' [illustrative, adjusted for numerical stability]'; }
  }
  return [`overrides applied: ${Object.keys(o.params || {}).join(', ')}`];
}

/** Time-dependent schedules (t comparisons, floor) get hmax = output step so no window is skipped. */
function solverOptions(spec) {
  const exprs = [...spec.vars.map((v) => v.expr), ...spec.flows.map((f) => f.rate)];
  const timed = exprs.some((e) => /\bfloor\(/.test(e) || /\bt\s*[<>]=?|[<>]=?\s*t\b/.test(e) || /\bifelse\s*\(\s*\(?\s*t\b/.test(e));
  if (timed) spec.hmax = spec.time[2];
}

/** Layered layout: compartments in columns by longest path from sources, components stacked. */
function layout(spec) {
  const names = spec.comps.map((c) => c.name);
  const idx = Object.fromEntries(names.map((n, i) => [n, i]));
  const edges = spec.flows.filter((f) => f.from && f.to && idx[f.from] != null && idx[f.to] != null).map((f) => [idx[f.from], idx[f.to]]);
  const n = names.length;
  // connected components (undirected)
  const comp = new Array(n).fill(-1);
  let nc = 0;
  for (let i = 0; i < n; i++) {
    if (comp[i] >= 0) continue;
    const stack = [i];
    comp[i] = nc;
    while (stack.length) {
      const u = stack.pop();
      for (const [a, b] of edges) {
        for (const [x, y] of [[a, b], [b, a]]) if (x === u && comp[y] < 0) { comp[y] = nc; stack.push(y); }
      }
    }
    nc++;
  }
  // DFS from nodes in declaration order to find back edges, then longest-path layering on the DAG
  const state = new Array(n).fill(0);
  const back = new Set();
  const visit = (u) => {
    state[u] = 1;
    for (const [a, b] of edges) {
      if (a !== u) continue;
      if (state[b] === 1) back.add(a + '>' + b);
      else if (state[b] === 0) visit(b);
    }
    state[u] = 2;
  };
  for (let i = 0; i < n; i++) if (!state[i]) visit(i);
  const dag = edges.filter(([a, b]) => !back.has(a + '>' + b));
  const layer = new Array(n).fill(0);
  for (let pass = 0; pass < n; pass++) for (const [a, b] of dag) if (layer[b] < layer[a] + 1) layer[b] = layer[a] + 1;
  // rows: within each component, stack nodes sharing a layer
  let rowBase = 0;
  const pos = new Array(n);
  for (let k = 0; k < nc; k++) {
    const members = names.map((_, i) => i).filter((i) => comp[i] === k);
    const perLayer = {};
    members.forEach((i) => { (perLayer[layer[i]] ||= []).push(i); });
    const height = Math.max(...Object.values(perLayer).map((l) => l.length));
    Object.values(perLayer).forEach((list) => list.forEach((i, r) => { pos[i] = [layer[i], rowBase + r + (height - list.length) / 2]; }));
    rowBase += height + 0.35;
  }
  spec.comps.forEach((c, i) => { c.col = pos[i][0]; c.row = +pos[i][1].toFixed(2); });
  // curve back edges and edges between nodes that are not adjacent columns on the same row
  spec.flows.forEach((f) => {
    if (!f.from || !f.to || idx[f.from] == null || idx[f.to] == null) return;
    const [a, b] = [idx[f.from], idx[f.to]];
    if (back.has(a + '>' + b)) f.curve = 60;
    else if (pos[a][1] === pos[b][1] && Math.abs(pos[a][0] - pos[b][0]) > 1) f.curve = 50;
  });
}

// ---------------------------------------------------------------------------
const files = fs.readdirSync(notesDir).filter((f) => f.endsWith('.md')).sort();
const specs = [];
const report = [];
const seen = new Set();
for (const file of files) {
  for (const b of blocks(fs.readFileSync(path.join(notesDir, file), 'utf8'))) {
    if (skip.has(b.id) || seen.has(b.id)) continue;
    const { spec, problems } = parseModel(b, file);
    if (!spec.comps.length || !spec.flows.length) continue; // not a model block
    seen.add(b.id);
    problems.push(...convertCounters(spec), ...applyOverrides(spec));
    solverOptions(spec);
    layout(spec);
    specs.push(spec);
    report.push({ file, id: b.id, problems });
  }
}

const header = `// GENERATED by cli/import-research.mjs from research_notes/ — do not edit by hand.
// Each entry is expanded by def() in library.js. Parameter descriptions carry the provenance
// tag recorded by the researcher: [verified], [recalled] or [illustrative].
`;
fs.writeFileSync(outFile, header + 'export default ' + JSON.stringify(specs, null, 1) + ';\n');
console.log(`${specs.length} model specs written to ${path.relative(root, outFile)}`);
report.filter((r) => r.problems.length).forEach((r) => console.log(`  ${r.id}: ${r.problems.join('; ')}`));
