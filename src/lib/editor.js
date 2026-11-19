// Interactive diagram editor. All R code comes from generator.generate(model);
// this file only edits the model and renders it.
import * as G from './generator.js';
import { LIBRARY, BY_ID, CATEGORIES } from './library.js';
import { indexModel, search } from './search.js';

const NW = 96;
const NH = 72;
const STORAGE_KEY = 'sir-diagram-to-r:model';
const PALETTE = ['#FFFFFF', '#E8C468', '#B71C1C', '#4E8A3E', '#6A2C9C', '#2457C5', '#D9A441', '#111111', '#00838F', '#8E0E0E'];
const SVGNS = 'http://www.w3.org/2000/svg';

const $ = (s) => document.querySelector(s);
const $$ = (s) => Array.from(document.querySelectorAll(s));
const clone = (o) => JSON.parse(JSON.stringify(o));
const esc = (s) => String(s ?? '').replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------
let model = null;
let selection = null;          // { kind: 'node' | 'flow', id }
let tool = 'select';
let tab = 'inspector';
let view = { x: 0, y: 0, k: 1 };
let pending = null;            // flow being drawn: { from: id|null, x, y }
let pointer = { x: 0, y: 0 };  // last pointer position (world coords)
let drag = null;
let undoStack = [];
let redoStack = [];
let result = { ok: false, code: null, errors: [], warnings: [] };

const svg = $('#canvas');
const world = $('#world');

// ---------------------------------------------------------------------------
// Model helpers
// ---------------------------------------------------------------------------
function normalize(m) {
  m = clone(m);
  m.version = G.FORMAT_VERSION;
  m.name ??= 'Model';
  m.description ??= '';
  m.reference ??= '';
  m.referenceUrl ??= '';
  m.time = Object.assign({ start: '0', stop: '100', step: '1', unit: 'days' }, m.time || {});
  m.population = Object.assign({ name: 'N' }, m.population || {});
  m.compartments = (m.compartments || []).map((c) => Object.assign({ description: '', init: '0', inN: true, color: '#FFFFFF', x: 0, y: 0 }, c));
  m.parameters = (m.parameters || []).map((p) => Object.assign({ symbol: '', value: '0', description: '' }, p));
  m.variables = (m.variables || []).map((v) => Object.assign({ symbol: '', expr: '0', description: '', output: false }, v));
  m.flows = (m.flows || []).map((f) => Object.assign({ kind: 'per_capita', rate: '0', output: '', curve: 0, from: null, to: null }, f));
  return m;
}

const comp = (id) => model.compartments.find((c) => c.id === id);
const flow = (id) => model.flows.find((f) => f.id === id);

function nextId(prefix, list) {
  let n = 0;
  list.forEach((x) => { const m = new RegExp('^' + prefix + '(\\d+)$').exec(x.id); if (m) n = Math.max(n, +m[1]); });
  return prefix + (n + 1);
}

function allNames() {
  return new Set([
    model.population.name,
    ...model.compartments.map((c) => c.name),
    ...model.parameters.map((p) => p.name),
    ...model.variables.map((v) => v.name),
    ...model.flows.map((f) => f.output).filter(Boolean),
  ]);
}

function uniqueName(base) {
  const taken = allNames();
  if (!taken.has(base)) return base;
  for (let i = 2; ; i++) if (!taken.has(base + i)) return base + i;
}

function symbolMap() {
  const map = {};
  model.parameters.forEach((p) => { if (p.symbol) map[p.name] = p.symbol; });
  model.variables.forEach((v) => { if (v.symbol) map[v.name] = v.symbol; });
  return map;
}

// ---------------------------------------------------------------------------
// History & persistence
// ---------------------------------------------------------------------------
function snapshot() {
  undoStack.push(JSON.stringify(model));
  if (undoStack.length > 200) undoStack.shift();
  redoStack = [];
}

function mutate(fn, { panels = true } = {}) {
  snapshot();
  fn();
  refresh({ panels });
}

function undo() {
  if (!undoStack.length) return;
  redoStack.push(JSON.stringify(model));
  model = JSON.parse(undoStack.pop());
  fixSelection();
  refresh();
}

function redo() {
  if (!redoStack.length) return;
  undoStack.push(JSON.stringify(model));
  model = JSON.parse(redoStack.pop());
  fixSelection();
  refresh();
}

function fixSelection() {
  if (selection && !(selection.kind === 'node' ? comp(selection.id) : flow(selection.id))) selection = null;
}

function persist() {
  try { localStorage.setItem(STORAGE_KEY, G.serialize(model)); } catch (e) { /* storage unavailable */ }
}

function encodeModel(m) {
  const bytes = new TextEncoder().encode(JSON.stringify(JSON.parse(G.serialize(m))));
  let bin = '';
  bytes.forEach((b) => { bin += String.fromCharCode(b); });
  return btoa(bin).replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, '');
}

function decodeModel(s) {
  const bin = atob(s.replace(/-/g, '+').replace(/_/g, '/'));
  const bytes = Uint8Array.from(bin, (ch) => ch.charCodeAt(0));
  return JSON.parse(new TextDecoder().decode(bytes));
}

function loadInitial() {
  const hash = location.hash.match(/^#m=(.+)$/);
  if (hash) {
    try { return normalize(decodeModel(hash[1])); } catch (e) { console.warn('Bad model link', e); }
  }
  const ex = location.hash.match(/^#example=([A-Za-z0-9_]+)$/);
  if (ex && BY_ID[ex[1]]) {
    history.replaceState(null, '', location.pathname);
    // keep the previously saved model reachable with Undo
    try { const saved = localStorage.getItem(STORAGE_KEY); if (saved) undoStack.push(JSON.stringify(normalize(JSON.parse(saved)))); } catch (e) { /* ignore */ }
    return normalize(BY_ID[ex[1]]);
  }
  try {
    const saved = localStorage.getItem(STORAGE_KEY);
    if (saved) return normalize(JSON.parse(saved));
  } catch (e) { /* ignore */ }
  return normalize(BY_ID.sir);
}

function replaceModel(m) {
  snapshot();
  model = normalize(m);
  selection = null;
  pending = null;
  refresh();
  fit();
}

// ---------------------------------------------------------------------------
// Refresh pipeline
// ---------------------------------------------------------------------------
let panelTimer = null;

function refresh({ panels = true } = {}) {
  result = G.generate(model);
  persist();
  renderCanvas();
  renderCode();
  renderStatus();
  renderToolbar();
  if (panels) schedulePanels(); else decorateErrors();
}

// Panels are re-rendered after the current event so focus can move first (Tab),
// then focus is restored to the same logical field.
function schedulePanels() {
  clearTimeout(panelTimer);
  panelTimer = setTimeout(() => {
    const active = document.activeElement;
    const key = active && active.dataset ? active.dataset.key : null;
    const selStart = active && 'selectionStart' in active ? active.selectionStart : null;
    renderPanels();
    if (key) {
      const el = document.querySelector(`[data-key="${CSS.escape(key)}"]`);
      if (el) {
        el.focus();
        if (selStart != null && 'setSelectionRange' in el) { try { el.setSelectionRange(selStart, selStart); } catch (e) { /* not a text input */ } }
      }
    }
  }, 0);
}

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------
function toWorld(evt) {
  const r = svg.getBoundingClientRect();
  return { x: (evt.clientX - r.left - view.x) / view.k, y: (evt.clientY - r.top - view.y) / view.k };
}

function rectExit(c, toward, gap = 6) {
  const dx = toward.x - c.x;
  const dy = toward.y - c.y;
  if (!dx && !dy) return { x: c.x, y: c.y };
  const s = Math.min((NW / 2 + gap) / Math.abs(dx || 1e-9), (NH / 2 + gap) / Math.abs(dy || 1e-9));
  return { x: c.x + dx * s, y: c.y + dy * s };
}

function flowGeom(f) {
  const A = f.from ? comp(f.from) : null;
  const B = f.to ? comp(f.to) : null;
  const a = A ? { x: A.x, y: A.y } : { x: f.ox ?? 0, y: f.oy ?? 0 };
  const b = B ? { x: B.x, y: B.y } : { x: f.ox ?? 0, y: f.oy ?? 0 };
  const dx = b.x - a.x;
  const dy = b.y - a.y;
  const len = Math.hypot(dx, dy) || 1;
  const n = { x: dy / len, y: -dx / len }; // left-hand normal (screen coords)
  const curve = +f.curve || 0;
  const ctrl = { x: (a.x + b.x) / 2 + n.x * curve * 2, y: (a.y + b.y) / 2 + n.y * curve * 2 };
  const s = A ? rectExit(a, ctrl) : a;
  const e = B ? rectExit(b, ctrl) : b;
  const mid = { x: 0.25 * s.x + 0.5 * ctrl.x + 0.25 * e.x, y: 0.25 * s.y + 0.5 * ctrl.y + 0.25 * e.y };
  return { s, e, ctrl, mid, n, a, b, len };
}

// ---------------------------------------------------------------------------
// Canvas rendering
// ---------------------------------------------------------------------------
function textSegments(segs, cls = 'sub') {
  // dy-based subscripts (works in every browser, unlike baseline-shift)
  let html = '';
  let shifted = false;
  segs.forEach((sg) => {
    if (sg.sub && !shifted) { html += `<tspan class="${cls}" dy="5">${esc(sg.text)}</tspan>`; shifted = true; }
    else if (!sg.sub && shifted) { html += `<tspan dy="-5">${esc(sg.text)}</tspan>`; shifted = false; }
    else if (sg.sub) html += `<tspan class="${cls}">${esc(sg.text)}</tspan>`;
    else html += `<tspan>${esc(sg.text)}</tspan>`;
  });
  return html;
}

function contrastText(hex) {
  const m = /^#?([0-9a-f]{6})$/i.exec(hex || '');
  if (!m) return '#111';
  const n = parseInt(m[1], 16);
  const lum = (0.2126 * ((n >> 16) & 255) + 0.7152 * ((n >> 8) & 255) + 0.0722 * (n & 255)) / 255;
  return lum > 0.6 ? '#111' : '#fff';
}

function errorIds() {
  const nodes = new Set();
  const flows = new Set();
  result.errors.forEach((e) => {
    if (e.kind === 'compartment') nodes.add(e.id);
    if (e.kind === 'flow') flows.add(e.id);
  });
  return { nodes, flows };
}

function renderCanvas() {
  world.setAttribute('transform', `translate(${view.x} ${view.y}) scale(${view.k})`);
  const bad = errorIds();
  const syms = symbolMap();
  let h = '';

  model.flows.forEach((f) => {
    if ((f.from && !comp(f.from)) || (f.to && !comp(f.to))) return;
    const g = flowGeom(f);
    const sel = selection?.kind === 'flow' && selection.id === f.id;
    const err = bad.flows.has(f.id);
    const marker = err ? 'arrow-err' : sel ? 'arrow-sel' : 'arrow';
    const d = `M${g.s.x},${g.s.y} Q${g.ctrl.x},${g.ctrl.y} ${g.e.x},${g.e.y}`;
    // label placement
    const t = { x: g.e.x - g.s.x, y: g.e.y - g.s.y };
    const tl = Math.hypot(t.x, t.y) || 1;
    let nx = t.y / tl;
    let ny = -t.x / tl;
    if (+f.curve) { const sgn = Math.sign(+f.curve); nx *= sgn; ny *= sgn; }
    else if (Math.abs(t.x) >= Math.abs(t.y)) { if (ny > 0) { nx = -nx; ny = -ny; } }
    else if (nx < 0) { nx = -nx; ny = -ny; }
    const lx = g.mid.x + nx * 16;
    const ly = g.mid.y + ny * 16 + (ny > 0.5 ? 12 : 0);
    const anchor = Math.abs(nx) > 0.6 ? (nx > 0 ? 'start' : 'end') : 'middle';
    const segs = G.prettySegments(f.rate, syms);
    const tracked = f.output ? `<tspan class="tracked" dx="4">[${esc(f.output)}]</tspan>` : '';
    h += `<g class="flow${sel ? ' selected' : ''}${err ? ' error' : ''}" data-flow="${esc(f.id)}">`
      + `<path class="hit" d="${d}"/>`
      + `<path class="line" d="${d}" marker-end="url(#${marker})"/>`
      + (!f.from ? `<circle class="origin" cx="${g.s.x}" cy="${g.s.y}" r="5"/>` : '')
      + `<text class="flow-label" x="${lx}" y="${ly}" text-anchor="${anchor}" dominant-baseline="middle">${textSegments(segs)}${tracked}</text>`
      + `</g>`;
  });

  model.compartments.forEach((c) => {
    const sel = selection?.kind === 'node' && selection.id === c.id;
    const err = bad.nodes.has(c.id);
    const src = pending && pending.from === c.id;
    const fg = contrastText(c.color);
    const segs = G.prettySegments(c.name, { [c.name]: c.name.replace(/^([A-Za-z]+)(\d+)$/, '$1_$2') });
    h += `<g class="node${sel ? ' selected' : ''}${err ? ' error' : ''}${src ? ' pending-src' : ''}" data-node="${esc(c.id)}" transform="translate(${c.x - NW / 2} ${c.y - NH / 2})">`
      + `<rect class="box" width="${NW}" height="${NH}" rx="3" fill="${esc(c.color)}"/>`
      + `<text x="${NW / 2}" y="${NH / 2}" fill="${fg}">${textSegments(segs)}</text>`
      + (c.inN === false ? `<text class="excluded" x="${NW / 2}" y="${NH + 12}">not in ${esc(model.population.name)}</text>` : '')
      + `</g>`;
  });

  // handles for the selected flow
  if (selection?.kind === 'flow' && tool === 'select') {
    const f = flow(selection.id);
    if (f) {
      const g = flowGeom(f);
      h += `<circle class="handle" data-handle="curve" cx="${g.mid.x}" cy="${g.mid.y}" r="6"><title>Drag to bend</title></circle>`;
      if (!f.from || !f.to) {
        const p = !f.from ? g.s : g.e;
        h += `<circle class="handle" data-handle="outside" cx="${p.x}" cy="${p.y}" r="6"><title>Drag to move</title></circle>`;
      }
    }
  }

  if (pending) {
    const from = pending.from ? comp(pending.from) : { x: pending.x, y: pending.y };
    const s = pending.from ? rectExit(from, pointer) : from;
    h += `<line class="rubber" x1="${s.x}" y1="${s.y}" x2="${pointer.x}" y2="${pointer.y}"/>`;
  }

  world.innerHTML = h;
  svg.classList.toggle('tool-compartment', tool === 'compartment');
  svg.classList.toggle('tool-flow', tool === 'flow');
  renderHint();
}

function renderHint() {
  let msg = '';
  if (tool === 'compartment') msg = 'Click to place a compartment · Esc to stop';
  else if (tool === 'flow' && !pending) msg = 'Click the source compartment (or empty space for an inflow from outside)';
  else if (tool === 'flow' && pending) msg = 'Click the target compartment (or empty space for an outflow) · Esc to cancel';
  $('#hint').textContent = msg;
}

function fit() {
  const pts = [];
  model.compartments.forEach((c) => { pts.push([c.x - NW / 2, c.y - NH / 2], [c.x + NW / 2, c.y + NH / 2 + 14]); });
  model.flows.forEach((f) => { if (!f.from || !f.to) pts.push([f.ox ?? 0, f.oy ?? 0]); const g = flowGeom(f); pts.push([g.mid.x, g.mid.y]); });
  const r = svg.getBoundingClientRect();
  if (!pts.length || !r.width) { view = { x: 0, y: 0, k: 1 }; renderCanvas(); return; }
  const xs = pts.map((p) => p[0]);
  const ys = pts.map((p) => p[1]);
  const minX = Math.min(...xs) - 50, maxX = Math.max(...xs) + 80;
  const minY = Math.min(...ys) - 50, maxY = Math.max(...ys) + 50;
  const statusH = $('#status').offsetHeight + 20;
  const k = Math.min(1.4, r.width / (maxX - minX), (r.height - statusH) / (maxY - minY));
  view = { k, x: (r.width - (maxX - minX) * k) / 2 - minX * k, y: (r.height - statusH - (maxY - minY) * k) / 2 - minY * k };
  renderCanvas();
}

// ---------------------------------------------------------------------------
// Canvas interaction
// ---------------------------------------------------------------------------
function addCompartment(x, y) {
  const preferred = ['S', 'E', 'I', 'R', 'V', 'A', 'C', 'H', 'D', 'Q'];
  const taken = allNames();
  let name = preferred.find((n) => !taken.has(n));
  if (!name) name = uniqueName('X1');
  const id = nextId('c', model.compartments);
  mutate(() => {
    model.compartments.push({
      id, name, description: '', init: '0', inN: true,
      color: PALETTE[model.compartments.length % PALETTE.length],
      x: Math.round(x / 10) * 10, y: Math.round(y / 10) * 10,
    });
    selection = { kind: 'node', id };
  });
  setTab('inspector');
}

function addFlow(from, to, pt) {
  const A = from ? comp(from) : null;
  const B = to ? comp(to) : null;
  const kind = from ? 'per_capita' : 'absolute';
  const pname = uniqueName(`r_${A ? A.name : 'in'}_${B ? B.name : 'out'}`);
  const id = nextId('f', model.flows);
  const twin = model.flows.find((f) => from && to && ((f.from === from && f.to === to) || (f.from === to && f.to === from)));
  mutate(() => {
    model.parameters.push({ name: pname, symbol: '', value: kind === 'absolute' ? '0' : '0.1', description: `rate ${A ? A.name : 'outside'} -> ${B ? B.name : 'outside'}` });
    const f = { id, from, to, kind, rate: pname, output: '', curve: twin ? 30 : 0 };
    if (!from || !to) { f.ox = Math.round(pt.x); f.oy = Math.round(pt.y); }
    model.flows.push(f);
    selection = { kind: 'flow', id };
  });
  setTab('inspector');
}

function deleteSelection() {
  if (!selection) return;
  if (selection.kind === 'node') {
    const id = selection.id;
    mutate(() => {
      model.compartments = model.compartments.filter((c) => c.id !== id);
      model.flows = model.flows.filter((f) => f.from !== id && f.to !== id);
      selection = null;
    });
  } else {
    const id = selection.id;
    mutate(() => { model.flows = model.flows.filter((f) => f.id !== id); selection = null; });
  }
}

function setTool(t) {
  tool = t;
  pending = null;
  renderToolbar();
  renderCanvas();
}

svg.addEventListener('pointerdown', (evt) => {
  if (evt.button !== 0) return;
  const p = toWorld(evt);
  pointer = p;
  const nodeEl = evt.target.closest('[data-node]');
  const flowEl = evt.target.closest('[data-flow]');
  const handleEl = evt.target.closest('[data-handle]');

  if (tool === 'compartment') {
    if (!nodeEl) addCompartment(p.x, p.y);
    setTool('select');
    return;
  }

  if (tool === 'flow') {
    const id = nodeEl ? nodeEl.dataset.node : null;
    if (!pending) {
      pending = { from: id, x: p.x, y: p.y };
      renderCanvas();
    } else {
      const from = pending.from;
      const start = { x: pending.x, y: pending.y };
      pending = null;
      if (from === id && from) { renderCanvas(); return; }      // clicked source again
      if (!from && !id) { renderCanvas(); return; }             // outside -> outside
      setTool('select');
      addFlow(from, id, from ? p : start);
    }
    return;
  }

  // select tool
  if (handleEl && selection?.kind === 'flow') {
    drag = { type: handleEl.dataset.handle, id: selection.id, before: JSON.stringify(model), moved: false };
  } else if (nodeEl) {
    const c = comp(nodeEl.dataset.node);
    selection = { kind: 'node', id: c.id };
    drag = { type: 'node', id: c.id, dx: p.x - c.x, dy: p.y - c.y, before: JSON.stringify(model), moved: false };
    refresh();
  } else if (flowEl) {
    selection = { kind: 'flow', id: flowEl.dataset.flow };
    refresh();
  } else {
    if (selection) { selection = null; refresh(); }
    drag = { type: 'pan', sx: evt.clientX, sy: evt.clientY, vx: view.x, vy: view.y };
    svg.classList.add('panning');
  }
  if (drag) svg.setPointerCapture(evt.pointerId);
  if (selection) setTab('inspector');
});

svg.addEventListener('pointermove', (evt) => {
  const p = toWorld(evt);
  pointer = p;
  if (pending) { renderCanvas(); return; }
  if (!drag) return;
  if (drag.type === 'pan') {
    view.x = drag.vx + evt.clientX - drag.sx;
    view.y = drag.vy + evt.clientY - drag.sy;
    renderCanvas();
    return;
  }
  drag.moved = true;
  if (drag.type === 'node') {
    const c = comp(drag.id);
    c.x = Math.round((p.x - drag.dx) / 10) * 10;
    c.y = Math.round((p.y - drag.dy) / 10) * 10;
  } else if (drag.type === 'curve') {
    const f = flow(drag.id);
    const g = flowGeom(Object.assign({}, f, { curve: 0 }));
    // signed distance of pointer from the straight line, along the normal
    let d = (p.x - (g.a.x + g.b.x) / 2) * g.n.x + (p.y - (g.a.y + g.b.y) / 2) * g.n.y;
    if (Math.abs(d) < 8) d = 0;
    f.curve = Math.round(d);
  } else if (drag.type === 'outside') {
    const f = flow(drag.id);
    f.ox = Math.round(p.x / 10) * 10;
    f.oy = Math.round(p.y / 10) * 10;
  }
  renderCanvas();
});

function endDrag() {
  if (!drag) return;
  svg.classList.remove('panning');
  if (drag.moved && drag.before) {
    undoStack.push(drag.before);
    redoStack = [];
    refresh({ panels: drag.type !== 'node' });
  }
  drag = null;
}
svg.addEventListener('pointerup', endDrag);
svg.addEventListener('pointercancel', endDrag);

svg.addEventListener('wheel', (evt) => {
  evt.preventDefault();
  const r = svg.getBoundingClientRect();
  const mx = evt.clientX - r.left;
  const my = evt.clientY - r.top;
  const k = Math.min(3, Math.max(0.25, view.k * Math.exp(-evt.deltaY * 0.0015)));
  view.x = mx - ((mx - view.x) / view.k) * k;
  view.y = my - ((my - view.y) / view.k) * k;
  view.k = k;
  renderCanvas();
}, { passive: false });

svg.addEventListener('dblclick', (evt) => {
  if (evt.target.closest('[data-node]') || evt.target.closest('[data-flow]')) {
    setTab('inspector');
    setTimeout(() => {
      const el = document.querySelector('#panel-inspector [data-focus]');
      if (el) { el.focus(); el.select?.(); }
    }, 10);
  }
});

document.addEventListener('keydown', (evt) => {
  const typing = /^(INPUT|TEXTAREA|SELECT)$/.test(document.activeElement?.tagName || '');
  const mod = evt.ctrlKey || evt.metaKey;
  if (mod && evt.key.toLowerCase() === 'z' && !typing) { evt.preventDefault(); evt.shiftKey ? redo() : undo(); return; }
  if (mod && evt.key.toLowerCase() === 'y' && !typing) { evt.preventDefault(); redo(); return; }
  if (typing || mod) return;
  if (evt.key === 'Escape') { setTool('select'); return; }
  if (evt.key === '/') { evt.preventDefault(); openLibrary(); return; }
  if (evt.key === 'Delete' || evt.key === 'Backspace') { evt.preventDefault(); deleteSelection(); return; }
  const k = evt.key.toLowerCase();
  if (k === 'v') setTool('select');
  if (k === 'c') setTool('compartment');
  if (k === 'f') setTool('flow');
});

// ---------------------------------------------------------------------------
// Panels
// ---------------------------------------------------------------------------
function setTab(t) {
  tab = t;
  $$('.tab').forEach((b) => b.classList.toggle('active', b.dataset.tab === t));
  $$('.panel').forEach((p) => p.classList.toggle('active', p.dataset.panel === t));
}

function input(key, value, { cls = '', placeholder = '', focus = false } = {}) {
  return `<input type="text" data-key="${esc(key)}" class="${cls}" value="${esc(value)}" placeholder="${esc(placeholder)}" spellcheck="false" autocomplete="off"${focus ? ' data-focus' : ''}/>`;
}
function fieldErr(key) { return `<div class="field-error" data-err="${esc(key)}"></div>`; }
function field(label, key, value, opts = {}) {
  return `<label class="field"><span>${label}</span>${input(key, value, opts)}${fieldErr(key)}</label>`;
}

function derivativeLine(c) {
  if (!result.code) return '';
  const line = result.code.split('\n').find((l) => l.trim().startsWith(`d${c.name} <- `));
  return line ? line.trim() : '';
}

function renderInspector() {
  const el = $('#panel-inspector');
  if (!selection) {
    el.innerHTML = `<h3>Nothing selected</h3>
      <p class="help">Use <b>+ Compartment</b> (C) to add boxes and <b>+ Flow</b> (F) to connect them.
      Click empty canvas while drawing a flow to make an inflow/outflow from outside the system (births, deaths).</p>
      <p class="help">Each flow has a <b>rate</b>. A <i>per-capita</i> rate is multiplied by the source compartment
      (an arrow labelled λ from S gives λ·S), which is how diagrams are normally labelled.
      An <i>absolute</i> rate is used as-is (e.g. births μ·N).</p>
      <p class="help">Define constants under <b>Parameters</b> and time-varying quantities such as the force of infection under
      <b>Variables</b>. The R code updates as you edit.</p>`;
    return;
  }
  if (selection.kind === 'node') {
    const c = comp(selection.id);
    const k = (f) => `comp:${c.id}:${f}`;
    const ins = model.flows.filter((f) => f.to === c.id);
    const outs = model.flows.filter((f) => f.from === c.id);
    const flowItem = (f) => {
      const other = f.to === c.id ? (f.from ? comp(f.from)?.name : 'outside') : (f.to ? comp(f.to)?.name : 'outside');
      return `<li><a data-select-flow="${esc(f.id)}">${f.to === c.id ? '← from' : '→ to'} ${esc(other)}</a> <code>${esc(f.rate)}</code></li>`;
    };
    el.innerHTML = `<h3>Compartment</h3>
      <div class="row">${field('Name (R variable)', k('name'), c.name, { focus: true })}
        <label class="field" style="flex:0 0 auto"><span>Colour</span><input type="color" data-key="${k('color')}" value="${esc(c.color)}"/></label></div>
      ${field('Description', k('description'), c.description, { placeholder: 'e.g. Susceptible' })}
      ${field('Initial value', k('init'), c.init, { cls: 'expr', placeholder: 'number, or e.g. 1000 - I' })}
      <p class="help">May use parameters and other compartments' initial values.</p>
      <label class="check"><input type="checkbox" data-key="${k('inN')}" ${c.inN !== false ? 'checked' : ''}/> Counts towards total population <code>${esc(model.population.name)}</code></label>
      <h4>Flows</h4>
      ${ins.length || outs.length ? `<ul class="flowlist">${ins.map(flowItem).join('')}${outs.map(flowItem).join('')}</ul>` : '<p class="help">No flows yet.</p>'}
      <h4>Equation</h4>
      <div class="eq">${esc(derivativeLine(c) || '(fix errors to see the equation)')}</div>
      <h4></h4>
      <button class="danger" data-action="delete">Delete compartment</button>`;
    return;
  }
  const f = flow(selection.id);
  const k = (x) => `flow:${f.id}:${x}`;
  const opts = (cur) => `<option value="">(outside)</option>` + model.compartments.map((c) => `<option value="${esc(c.id)}"${cur === c.id ? ' selected' : ''}>${esc(c.name)}</option>`).join('');
  const src = f.from ? comp(f.from)?.name : null;
  const term = f.kind === 'absolute' ? f.rate : `(${f.rate}) × ${src ?? '?'}`;
  el.innerHTML = `<h3>Flow</h3>
    <div class="row">
      <label class="field"><span>From</span><select class="full" data-key="${k('from')}">${opts(f.from)}</select>${fieldErr(k('from'))}</label>
      <label class="field"><span>To</span><select class="full" data-key="${k('to')}">${opts(f.to)}</select>${fieldErr(k('to'))}</label>
    </div>
    <span class="small muted">Rate type</span>
    <div class="radio-row">
      <label><input type="radio" name="kind" data-key="${k('kind')}" value="per_capita" ${f.kind !== 'absolute' ? 'checked' : ''}/> per capita (× source)</label>
      <label><input type="radio" name="kind" data-key="${k('kind')}" value="absolute" ${f.kind === 'absolute' ? 'checked' : ''}/> absolute (total)</label>
    </div>
    ${fieldErr(k('kind'))}
    ${field('Rate', k('rate'), f.rate, { cls: 'expr', focus: true, placeholder: 'e.g. lam, nu * (1 - gamma)' })}
    <p class="help">Flow per unit time = <code>${esc(term)}</code></p>
    ${field('Record as output (optional)', k('output'), f.output, { placeholder: 'e.g. incidence' })}
    <p class="help">Named flows are returned by the model function and appear as columns in <code>out</code>.</p>
    <div class="row">
      <label class="field"><span>Curvature</span><input type="text" data-key="${k('curve')}" value="${esc(f.curve)}"/></label>
      <div></div>
    </div>
    <button data-action="reverse">Reverse direction</button>
    <button class="danger" data-action="delete">Delete flow</button>`;
}

function renderParameters() {
  const el = $('#panel-parameters');
  el.innerHTML = `<h3>Parameters</h3>
    <p class="help">Constants passed to <code>ode()</code>. Values are numbers or arithmetic (e.g. <code>1 / (50 * 52)</code>).
    The symbol is only used for diagram labels (<code>tau_A</code> shows as τ<sub>A</sub>).</p>
    ${model.parameters.length ? greekPad('param') : ''}
    ${model.parameters.map((p, i) => `<div class="card">
      <div class="row">
        ${input(`param:${i}:name`, p.name, { placeholder: 'name' })}
        ${input(`param:${i}:value`, p.value, { cls: 'expr', placeholder: 'value' })}
        <span class="order"><button data-move-param="${i}:-1" title="Move up">↑</button><button data-move-param="${i}:1" title="Move down">↓</button></span>
        <button class="del danger" data-del-param="${i}" title="Delete">✕</button>
      </div>
      <div class="row">
        ${input(`param:${i}:symbol`, p.symbol, { placeholder: 'symbol: ' + G.autoSymbol(p.name || '') })}
        ${input(`param:${i}:description`, p.description, { placeholder: 'description' })}
      </div>
      ${fieldErr(`param:${i}:name`)}${fieldErr(`param:${i}:value`)}
    </div>`).join('')}
    <button class="primary" data-action="add-param">+ Add parameter</button>`;
}

// Greek keypad for the Parameters and Variables panels. Letters go into the last-focused
// symbol/description field of that panel (diagram labels only: R names and expressions stay ASCII).
const GREEK_LOWER = 'αβγδεζηθικλμνξοπρστυφχψω';
const GREEK_UPPER = 'ΓΔΘΛΞΠΣΦΨΩ';
const GREEK_NAMES = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu',
  'nu', 'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', 'upsilon', 'phi', 'chi', 'psi', 'omega'];
const greekName = (ch) => {
  const i = GREEK_LOWER.indexOf(ch.toLowerCase());
  const n = GREEK_NAMES[i];
  return ch === GREEK_LOWER[i] ? n : n[0].toUpperCase() + n.slice(1);
};
const greekTarget = { param: null, var: null }; // data-key of the last focused symbol/description input per panel

function greekPad(scope) {
  const btn = (ch) => `<button type="button" data-greek="${ch}" data-name="${greekName(ch)}" aria-label="${greekName(ch)}">${ch}</button>`;
  return `<div class="greek-pad" data-greek-scope="${scope}" aria-label="Greek letters">
    <span class="small muted">Insert into the selected symbol or description:</span>
    <div class="greek-keys">${[...GREEK_LOWER].map(btn).join('')}</div>
    <div class="greek-keys">${[...GREEK_UPPER].map(btn).join('')}</div>
  </div>`;
}

function insertGreek(ch, scope) {
  const last = (scope === 'param' ? model.parameters : model.variables).length - 1;
  const saved = greekTarget[scope];
  const key = saved && document.querySelector(`[data-key="${CSS.escape(saved)}"]`) ? saved : `${scope}:${last}:symbol`;
  const el = document.querySelector(`[data-key="${CSS.escape(key)}"]`);
  if (!el) return;
  const start = el.selectionStart ?? el.value.length;
  const end = el.selectionEnd ?? start;
  el.value = el.value.slice(0, start) + ch + el.value.slice(end);
  el.focus();
  el.setSelectionRange(start + ch.length, start + ch.length);
  greekTarget[scope] = key;
  applyField(key, el);
}

function renderVariables() {
  const el = $('#panel-variables');
  el.innerHTML = `<h3>Variables (derived quantities)</h3>
    <p class="help">Recomputed at every time step, in dependency order, from compartments, parameters, other variables,
    <code>${esc(model.population.name)}</code> (total population) and <code>t</code> (time). Use them for forces of infection,
    seasonality and interventions, e.g. <code>(t &gt;= t_start) * coverage</code>.</p>
    ${model.variables.length ? greekPad('var') : ''}
    ${model.variables.map((v, i) => `<div class="card">
      <div class="row">
        ${input(`var:${i}:name`, v.name, { placeholder: 'name' })}
        ${input(`var:${i}:symbol`, v.symbol, { placeholder: 'symbol: ' + G.autoSymbol(v.name || '') })}
        <span class="order"><button data-move-var="${i}:-1" title="Move up">↑</button><button data-move-var="${i}:1" title="Move down">↓</button></span>
        <button class="del danger" data-del-var="${i}" title="Delete">✕</button>
      </div>
      <div class="row">${input(`var:${i}:expr`, v.expr, { cls: 'expr', placeholder: 'expression' })}</div>
      <div class="row">
        ${input(`var:${i}:description`, v.description, { placeholder: 'description' })}
        <label class="check" style="flex:0 0 auto;margin:0"><input type="checkbox" data-key="var:${i}:output" ${v.output ? 'checked' : ''}/> output</label>
      </div>
      ${fieldErr(`var:${i}:name`)}${fieldErr(`var:${i}:expr`)}
    </div>`).join('')}
    <button class="primary" data-action="add-var">+ Add variable</button>`;
}

function renderModelPanel() {
  const el = $('#panel-model');
  el.innerHTML = `<h3>Model</h3>
    ${field('Name', 'model:name', model.name)}
    <label class="field"><span>Description (becomes the header comment)</span>
      <textarea data-key="model:description" rows="3">${esc(model.description)}</textarea></label>
    ${field('Reference', 'model:reference', model.reference, { placeholder: 'e.g. Keeling & Rohani (2008), §2.1' })}
    <label class="field"><span>Reference link (DOI or URL)</span>${input('model:referenceUrl', model.referenceUrl, { placeholder: 'https://doi.org/…' })}
      ${/^https?:\/\//.test(model.referenceUrl) ? `<a class="small" href="${esc(model.referenceUrl)}" target="_blank" rel="noopener noreferrer">Open source ↗</a>` : ''}</label>
    <h4>Time</h4>
    <div class="row">
      ${field('Start', 'model:time.start', model.time.start, { cls: 'expr' })}
      ${field('Stop', 'model:time.stop', model.time.stop, { cls: 'expr' })}
      ${field('Step', 'model:time.step', model.time.step, { cls: 'expr' })}
      ${field('Unit', 'model:time.unit', model.time.unit)}
    </div>
    ${field('Maximum solver step (optional)', 'model:solver.hmax', model.solver?.hmax ?? '', { cls: 'expr', placeholder: 'e.g. 1, for pulsed treatment or campaigns' })}
    <p class="help">Stops the solver stepping over short events such as treatment rounds (deSolve <code>hmax</code>).</p>
    <h4>Population</h4>
    ${field('Name of the total population variable', 'model:population', model.population.name)}
    <p class="help">Sum of all compartments marked "counts towards total population". Untick it for compartments such as deaths.</p>
    ${fieldErr('model:compartments')}${fieldErr('model:init')}${fieldErr('model:variables')}`;
}

function renderPanels() {
  renderInspector();
  renderParameters();
  renderVariables();
  renderModelPanel();
  decorateErrors();
}

function errorKey(e) {
  if (e.kind === 'compartment') return `comp:${e.id}:${e.field}`;
  if (e.kind === 'flow') return `flow:${e.id}:${e.field}`;
  if (e.kind === 'parameter') return `param:${e.id}:${e.field}`;
  if (e.kind === 'variable') return `var:${e.id}:${e.field}`;
  return `model:${e.field === 'population' ? 'population' : e.field}`;
}

function decorateErrors() {
  const byKey = {};
  result.errors.forEach((e) => { const k = errorKey(e); (byKey[k] ||= []).push(e.message); });
  $$('[data-key]').forEach((el) => {
    const msgs = byKey[el.dataset.key];
    el.classList.toggle('invalid', !!msgs);
    el.title = msgs ? msgs.join('\n') : '';
  });
  $$('[data-err]').forEach((el) => { el.textContent = (byKey[el.dataset.err] || []).join(' · '); });
}

// Apply a field change identified by data-key.
function applyField(key, el) {
  const parts = key.split(':');
  const value = el.type === 'checkbox' ? el.checked : el.value;
  const text = typeof value === 'string' ? value.trim() : value;
  const [scope, id, prop] = parts;

  mutate(() => {
    if (scope === 'comp') {
      const c = comp(id);
      if (prop === 'name' && text !== c.name) { G.renameReferences(model, c.name, text); c.name = text; }
      else if (prop === 'inN') c.inN = value;
      else c[prop] = text;
    } else if (scope === 'flow') {
      const f = flow(id);
      if (prop === 'from' || prop === 'to') {
        const v = text || null;
        f[prop] = v;
        if (!f.from || !f.to) {
          const other = comp(f.from || f.to);
          if (other && (f.ox == null || f.oy == null)) { f.ox = other.x + (f.from ? 0 : -150); f.oy = other.y + (f.from ? 150 : 0); }
        }
      } else if (prop === 'kind') f.kind = value;
      else if (prop === 'curve') f.curve = Number(text) || 0;
      else f[prop] = text;
    } else if (scope === 'param') {
      const p = model.parameters[+id];
      if (prop === 'name' && text !== p.name) { G.renameReferences(model, p.name, text); p.name = text; }
      else p[prop] = text;
    } else if (scope === 'var') {
      const v = model.variables[+id];
      if (prop === 'name' && text !== v.name) { G.renameReferences(model, v.name, text); v.name = text; }
      else if (prop === 'output') v.output = value;
      else v[prop] = text;
    } else if (scope === 'model') {
      const k = parts[1];
      if (k === 'name') model.name = text;
      else if (k === 'description') model.description = String(value);
      else if (k === 'reference') model.reference = text;
      else if (k === 'referenceUrl') model.referenceUrl = text;
      else if (k === 'population') {
        if (text !== model.population.name) { G.renameReferences(model, model.population.name, text); model.population.name = text; }
      } else if (k.startsWith('time.')) model.time[k.slice(5)] = text;
      else if (k === 'solver.hmax') { if (text) model.solver = { hmax: text }; else delete model.solver; }
    }
  });
}

document.querySelector('.side').addEventListener('change', (evt) => {
  const el = evt.target;
  if (el.dataset && el.dataset.key) applyField(el.dataset.key, el);
});

document.querySelector('.side').addEventListener('focusin', (evt) => {
  const k = evt.target.dataset?.key;
  const m = k && /^(param|var):\d+:(symbol|description)$/.exec(k);
  if (m) greekTarget[m[1]] = k;
});
// Keep focus (and the caret) in the text field when a Greek key is pressed.
document.querySelector('.side').addEventListener('mousedown', (evt) => {
  if (evt.target.closest('[data-greek]')) evt.preventDefault();
});

document.querySelector('.side').addEventListener('click', (evt) => {
  const t = evt.target.closest('button, a');
  if (!t) return;
  if (t.dataset.greek) { insertGreek(t.dataset.greek, t.closest('[data-greek-scope]').dataset.greekScope); return; }
  if (t.dataset.tab) { setTab(t.dataset.tab); return; }
  if (t.dataset.selectFlow) { selection = { kind: 'flow', id: t.dataset.selectFlow }; refresh(); return; }
  if (t.dataset.delParam != null) { const i = +t.dataset.delParam; mutate(() => model.parameters.splice(i, 1)); return; }
  if (t.dataset.delVar != null) { const i = +t.dataset.delVar; mutate(() => model.variables.splice(i, 1)); return; }
  if (t.dataset.moveParam || t.dataset.moveVar) {
    const list = t.dataset.moveParam ? model.parameters : model.variables;
    const [i, d] = (t.dataset.moveParam || t.dataset.moveVar).split(':').map(Number);
    if (i + d < 0 || i + d >= list.length) return;
    mutate(() => { const [x] = list.splice(i, 1); list.splice(i + d, 0, x); });
    return;
  }
  const action = t.dataset.action;
  if (action === 'delete') deleteSelection();
  if (action === 'reverse') {
    const f = flow(selection.id);
    mutate(() => { [f.from, f.to] = [f.to, f.from]; if (!f.from && f.kind === 'per_capita') f.kind = 'absolute'; });
  }
  if (action === 'add-param') {
    mutate(() => model.parameters.push({ name: uniqueName('p1'), symbol: '', value: '0', description: '' }));
    setTimeout(() => document.querySelector(`[data-key="param:${model.parameters.length - 1}:name"]`)?.select(), 20);
  }
  if (action === 'add-var') {
    mutate(() => model.variables.push({ name: uniqueName('v1'), symbol: '', expr: '0', description: '', output: false }));
    setTimeout(() => document.querySelector(`[data-key="var:${model.variables.length - 1}:name"]`)?.select(), 20);
  }
});

// ---------------------------------------------------------------------------
// Status bar, code view, toolbar
// ---------------------------------------------------------------------------
function issueTarget(e) {
  if (e.kind === 'compartment') return () => { selection = { kind: 'node', id: e.id }; setTab('inspector'); refresh(); };
  if (e.kind === 'flow') return () => { selection = { kind: 'flow', id: e.id }; setTab('inspector'); refresh(); };
  if (e.kind === 'parameter') return () => setTab('parameters');
  if (e.kind === 'variable') return () => setTab('variables');
  return () => setTab('model');
}

function renderStatus() {
  const el = $('#status');
  const { errors, warnings } = result;
  if (!errors.length && !warnings.length) {
    el.innerHTML = `<span class="ok">✓ Model is valid.</span> <span class="muted">${model.compartments.length} compartments, ${model.flows.length} flows. R code is up to date.</span>`;
    return;
  }
  const items = errors.map((e) => ({ e, cls: 'error' })).concat(warnings.map((e) => ({ e, cls: 'warning' })));
  el.innerHTML = `<b>${errors.length ? `${errors.length} error${errors.length > 1 ? 's' : ''}` : '✓ No errors'}</b>${warnings.length ? `, ${warnings.length} warning${warnings.length > 1 ? 's' : ''}` : ''}
    <ul>${items.slice(0, 12).map(({ e, cls }, i) => `<li class="${cls}"><span class="jump" data-issue="${i}">${esc(e.message)}</span>${e.unknown ? `<button data-add-param="${esc(e.unknown)}">+ parameter ${esc(e.unknown)}</button><button data-add-var="${esc(e.unknown)}">+ variable</button>` : ''}</li>`).join('')}</ul>`;
  el.querySelectorAll('[data-issue]').forEach((s) => s.addEventListener('click', () => issueTarget(items[+s.dataset.issue].e)()));
  el.querySelectorAll('[data-add-param]').forEach((b) => b.addEventListener('click', () => {
    const name = b.dataset.addParam;
    mutate(() => model.parameters.push({ name, symbol: '', value: '0', description: '' }));
    setTab('parameters');
  }));
  el.querySelectorAll('[data-add-var]').forEach((b) => b.addEventListener('click', () => {
    const name = b.dataset.addVar;
    mutate(() => model.variables.push({ name, symbol: '', expr: '0', description: '', output: false }));
    setTab('variables');
  }));
}

let lastGoodCode = '';
function highlightR(code) {
  const re = /(#[^\n]*)|("(?:[^"\\]|\\.)*")|(<-)|\b(function|if|else|for|while|return|TRUE|FALSE|NULL)\b|\b([A-Za-z_.][A-Za-z0-9_.]*)(?=\()|\b(\d+\.?\d*(?:[eE][+-]?\d+)?)\b/g;
  let out = '';
  let i = 0;
  let m;
  while ((m = re.exec(code))) {
    out += esc(code.slice(i, m.index));
    const cls = m[1] ? 'com' : m[2] ? 'str' : m[3] ? 'assign' : m[4] ? 'kw' : m[5] ? 'fn' : 'num';
    out += `<span class="tok-${cls}">${esc(m[0])}</span>`;
    i = m.index + m[0].length;
  }
  return out + esc(code.slice(i));
}

function renderCode() {
  const pre = $('#code');
  const badge = $('#code-badge');
  if (result.ok) {
    lastGoodCode = result.code;
    pre.innerHTML = highlightR(result.code);
    pre.classList.remove('stale');
    $('#code-errors').innerHTML = '';
    badge.className = 'badge ok';
    badge.textContent = '✓';
  } else {
    pre.innerHTML = lastGoodCode ? highlightR(lastGoodCode) : '';
    pre.classList.add('stale');
    $('#code-errors').innerHTML = `<b>Fix ${result.errors.length} error${result.errors.length > 1 ? 's' : ''} to update the code</b> (last valid version shown):<br>`
      + result.errors.slice(0, 5).map((e) => esc(e.message)).join('<br>');
    badge.className = 'badge bad';
    badge.textContent = String(result.errors.length);
  }
  $('#btn-copy').disabled = !result.ok;
  $('#btn-download').disabled = !result.ok;
}

function renderToolbar() {
  $$('.tool').forEach((b) => b.classList.toggle('active', b.dataset.tool === tool));
  $('#btn-undo').disabled = !undoStack.length;
  $('#btn-redo').disabled = !redoStack.length;
}

function download(name, text, type) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const a = document.createElement('a');
  a.href = url;
  a.download = name;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function fileBase() { return (model.name || 'model').replace(/[^A-Za-z0-9_-]+/g, '_'); }

function flash(btn, text) {
  const old = btn.textContent;
  btn.textContent = text;
  setTimeout(() => { btn.textContent = old; }, 1300);
}

$$('.tool').forEach((b) => b.addEventListener('click', () => setTool(b.dataset.tool)));
$$('.tab').forEach((b) => b.addEventListener('click', () => setTab(b.dataset.tab)));
$('#btn-undo').addEventListener('click', undo);
$('#btn-redo').addEventListener('click', redo);
$('#btn-fit').addEventListener('click', fit);
// ---- Library picker ----------------------------------------------------------
const libDocs = LIBRARY.map(indexModel);
const libDialog = $('#lib-dialog');
let libCat = null;
let libHits = [];
let libActive = 0;

function renderLibCats() {
  const counts = {};
  LIBRARY.forEach((m) => { counts[m.category] = (counts[m.category] || 0) + 1; });
  $('#lib-cats').innerHTML = [`<button class="chip-btn${libCat ? '' : ' on'}" data-cat="">All ${LIBRARY.length}</button>`]
    .concat(CATEGORIES.map((c) => `<button class="chip-btn${libCat === c ? ' on' : ''}" data-cat="${esc(c)}">${esc(c)} ${counts[c]}</button>`)).join('');
}

function renderLibResults() {
  libHits = search(libDocs, $('#lib-q').value, { category: libCat });
  libActive = Math.min(libActive, Math.max(0, libHits.length - 1));
  $('#lib-count').textContent = `${libHits.length} model${libHits.length === 1 ? '' : 's'}`;
  $('#lib-results').innerHTML = libHits.length ? libHits.map(({ model: m }, i) => `
    <li class="lib-hit${i === libActive ? ' active' : ''}" role="option" data-i="${i}">
      <div class="lib-hit-top"><b>${esc(m.name)}</b><span class="muted small">${esc(m.category)}</span></div>
      <div class="lib-hit-comps">${m.compartments.map((c) => `<span class="chip" style="background:${esc(c.color)};color:${contrastText(c.color)}">${esc(c.name)}</span>`).join('')}</div>
      <div class="muted small lib-hit-desc">${esc(m.description.split(/\r?\n/)[0])}</div>
      ${m.reference ? `<div class="small lib-hit-ref">📖 ${esc(m.reference)}</div>` : ''}
    </li>`).join('') : '<li class="muted lib-empty">No models match. Try fewer words.</li>';
  $('#lib-results .lib-hit.active')?.scrollIntoView({ block: 'nearest' });
}

function openLibrary() {
  if (libDialog.open) return;
  libDialog.showModal();
  renderLibCats();
  renderLibResults();
  $('#lib-q').focus();
  $('#lib-q').select();
}

function pickLibrary(i) {
  const hit = libHits[i];
  if (!hit) return;
  libDialog.close();
  replaceModel(hit.model);
}

$('#btn-library').addEventListener('click', openLibrary);
$('#lib-close').addEventListener('click', () => libDialog.close());
libDialog.addEventListener('click', (e) => { if (e.target === libDialog) libDialog.close(); });
$('#lib-q').addEventListener('input', () => { libActive = 0; renderLibResults(); });
$('#lib-q').addEventListener('keydown', (e) => {
  if (e.key === 'ArrowDown') { e.preventDefault(); libActive = Math.min(libActive + 1, libHits.length - 1); renderLibResults(); }
  else if (e.key === 'ArrowUp') { e.preventDefault(); libActive = Math.max(libActive - 1, 0); renderLibResults(); }
  else if (e.key === 'Enter') { e.preventDefault(); pickLibrary(libActive); }
});
$('#lib-cats').addEventListener('click', (e) => {
  const b = e.target.closest('[data-cat]');
  if (!b) return;
  libCat = b.dataset.cat || null;
  libActive = 0;
  renderLibCats();
  renderLibResults();
  $('#lib-q').focus();
});
$('#lib-results').addEventListener('click', (e) => {
  const li = e.target.closest('[data-i]');
  if (li) pickLibrary(+li.dataset.i);
});
$('#btn-new').addEventListener('click', () => replaceModel({ name: 'Model', compartments: [], parameters: [], variables: [], flows: [] }));
$('#btn-open').addEventListener('click', () => $('#file-open').click());
$('#file-open').addEventListener('change', async (e) => {
  const file = e.target.files[0];
  e.target.value = '';
  if (!file) return;
  try { replaceModel(JSON.parse(await file.text())); } catch (err) { $('#status').textContent = 'Could not read that file: ' + err.message; }
});
$('#btn-save').addEventListener('click', () => download(fileBase() + '.json', G.serialize(model), 'application/json'));
$('#btn-download').addEventListener('click', () => { if (result.ok) download(fileBase() + '.R', result.code, 'text/plain'); });
$('#btn-copy').addEventListener('click', async (e) => {
  if (!result.ok) return;
  try { await navigator.clipboard.writeText(result.code); flash(e.target, 'Copied'); } catch (err) { flash(e.target, 'Copy failed'); }
});
$('#btn-link').addEventListener('click', async (e) => {
  const url = location.origin + location.pathname + '#m=' + encodeModel(model);
  history.replaceState(null, '', url);
  try { await navigator.clipboard.writeText(url); flash(e.target, 'Link copied'); } catch (err) { flash(e.target, 'Link in address bar'); }
});
window.addEventListener('resize', () => renderCanvas());

// ---------------------------------------------------------------------------
// Boot
// ---------------------------------------------------------------------------
model = loadInitial();
setTab('inspector');
refresh();
requestAnimationFrame(fit);
