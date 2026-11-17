/*
 * Compartmental-model diagram -> R (deSolve) code generator.
 *
 * Pure and deterministic: the same model JSON always yields byte-identical R
 * code (no timestamps, no randomness, ordering taken only from the model's
 * arrays). Plain ES module shared by the web editor, the CLI and the tests.
 *
 * Model format (version 1):
 * {
 *   version: 1,
 *   name: "SEIR",
 *   description: "", reference: "", referenceUrl: "",
 *   time: { start: "0", stop: "100", step: "1", unit: "days" },
 *   population: { name: "N" },
 *   solver: { hmax: "1" },            (optional: maximum solver step, for pulsed schedules)
 *   compartments: [{ id, name, description, init, inN, color, x, y }],
 *   parameters:   [{ name, symbol, value, description }],
 *   variables:    [{ name, symbol, expr, description, output }],
 *   flows:        [{ id, from, to, kind: "per_capita"|"absolute", rate,
 *                    output, curve, ox, oy }]
 * }
 * flows[].from / flows[].to are compartment ids, or null for "outside the
 * system" (births, immigration, deaths, emigration).
 *   per_capita flow  = rate * <source compartment>   (what diagrams usually label)
 *   absolute flow    = rate                          (e.g. births mu * N)
 */

var FORMAT_VERSION = 1;

var R_RESERVED = [
  'if', 'else', 'repeat', 'while', 'function', 'for', 'in', 'next', 'break',
  'TRUE', 'FALSE', 'NULL', 'Inf', 'NaN', 'NA', 'NA_integer_', 'NA_real_',
  'NA_character_', 'T', 'F',
];
// Names the generated script itself uses.
var INTERNAL = ['t', 'c', 'list', 'with', 'parameters', 'state', 'times', 'out', 'ode',
  'library', 'seq', 'as', 'matplot', 'legend'];
var R_FUNCS = ['exp', 'log', 'log10', 'log2', 'log1p', 'expm1', 'sqrt', 'abs', 'sin', 'cos', 'tan',
  'min', 'max', 'pmin', 'pmax', 'ifelse', 'floor', 'ceiling', 'round', 'sign'];
var R_CONSTS = ['pi'];
var BINARY_OPS = ['%/%', '%%', '<=', '>=', '==', '!=', '&&', '||', '+', '-', '*', '/', '^', '<', '>', '&', '|'];
var ALL_OPS = BINARY_OPS.concat(['!']);
// Binary operators binding looser than "*": an expression containing one of
// these at top level must be parenthesised before multiplying it.
var LOW_PRECEDENCE = ['+', '-', '<=', '>=', '==', '!=', '&&', '||', '<', '>', '&', '|'];
var IDENT_RE = /^[A-Za-z][A-Za-z0-9_.]*$/;

function has(arr, v) { return arr.indexOf(v) !== -1; }

// ---------------------------------------------------------------------------
// Expressions
// ---------------------------------------------------------------------------

function tokenize(src) {
  var toks = [];
  var i = 0;
  src = String(src == null ? '' : src);
  while (i < src.length) {
    var ch = src[i];
    if (/\s/.test(ch)) { i++; continue; }
    var rest = src.slice(i);
    var m = /^(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?/.exec(rest);
    if (m) { toks.push({ t: 'num', v: m[0] }); i += m[0].length; continue; }
    m = /^[A-Za-z][A-Za-z0-9_.]*/.exec(rest);
    if (m) { toks.push({ t: 'id', v: m[0] }); i += m[0].length; continue; }
    if (ch === '(') { toks.push({ t: 'lp', v: '(' }); i++; continue; }
    if (ch === ')') { toks.push({ t: 'rp', v: ')' }); i++; continue; }
    if (ch === ',') { toks.push({ t: 'comma', v: ',' }); i++; continue; }
    var op = null;
    for (var k = 0; k < ALL_OPS.length; k++) {
      if (rest.indexOf(ALL_OPS[k]) === 0) { op = ALL_OPS[k]; break; }
    }
    if (op) { toks.push({ t: 'op', v: op }); i += op.length; continue; }
    throw new Error('Unexpected character "' + ch + '"');
  }
  return toks;
}

/**
 * Parse-check an expression. Returns
 *   { ok, error, tokens, idents (non-call identifiers, in order, unique),
 *     funcs, code (canonical formatting), lowPrecedence }
 */
function analyze(src) {
  var res = { ok: false, error: null, tokens: [], idents: [], funcs: [], code: '', lowPrecedence: false };
  var toks;
  try { toks = tokenize(src); } catch (e) { res.error = e.message; return res; }
  res.tokens = toks;
  if (!toks.length) { res.error = 'Expression is empty'; return res; }

  var expectOperand = true;
  var stack = [];
  var out = '';
  for (var i = 0; i < toks.length; i++) {
    var tk = toks[i];
    if (expectOperand) {
      if (tk.t === 'num') {
        out += tk.v; expectOperand = false;
      } else if (tk.t === 'id') {
        var next = toks[i + 1];
        if (next && next.t === 'lp') {
          if (!has(res.funcs, tk.v)) res.funcs.push(tk.v);
          stack.push('call'); out += tk.v + '('; i++;
          if (toks[i + 1] && toks[i + 1].t === 'rp') { res.error = tk.v + '() needs an argument'; return res; }
        } else {
          if (!has(res.idents, tk.v)) res.idents.push(tk.v);
          out += tk.v; expectOperand = false;
        }
      } else if (tk.t === 'lp') {
        stack.push('paren'); out += '(';
      } else if (tk.t === 'op' && (tk.v === '-' || tk.v === '+' || tk.v === '!')) {
        if (stack.length === 0) res.lowPrecedence = true;
        out += tk.v;
      } else {
        res.error = 'Unexpected "' + tk.v + '"'; return res;
      }
    } else {
      if (tk.t === 'op' && tk.v !== '!') {
        if (stack.length === 0 && has(LOW_PRECEDENCE, tk.v)) res.lowPrecedence = true;
        out += tk.v === '^' ? '^' : ' ' + tk.v + ' ';
        expectOperand = true;
      } else if (tk.t === 'rp') {
        if (!stack.length) { res.error = 'Unmatched ")"'; return res; }
        stack.pop(); out += ')';
      } else if (tk.t === 'comma') {
        if (stack[stack.length - 1] !== 'call') { res.error = 'Unexpected ","'; return res; }
        out += ', '; expectOperand = true;
      } else {
        res.error = 'Missing operator before "' + tk.v + '"'; return res;
      }
    }
  }
  if (expectOperand) { res.error = 'Expression ends unexpectedly'; return res; }
  if (stack.length) { res.error = 'Unclosed "("'; return res; }
  res.ok = true;
  res.code = out;
  return res;
}

function paren(a) { return a.lowPrecedence ? '(' + a.code + ')' : a.code; }

/** Rename identifier `from` to `to` in an expression, preserving everything else. */
function renameInExpr(src, from, to) {
  var toks;
  try { toks = tokenize(src); } catch (e) { return src; }
  if (!toks.some(function (t) { return t.t === 'id' && t.v === from; })) return src;
  // Rebuild from source positions so user spacing is preserved.
  var s = String(src);
  var outStr = '';
  var i = 0;
  var re = /[A-Za-z][A-Za-z0-9_.]*|(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?/g;
  var m;
  while ((m = re.exec(s))) {
    outStr += s.slice(i, m.index);
    outStr += m[0] === from ? to : m[0];
    i = m.index + m[0].length;
  }
  return outStr + s.slice(i);
}

// ---------------------------------------------------------------------------
// Display helpers (used by the diagram, not by code generation)
// ---------------------------------------------------------------------------

var GREEK = {
  alpha: 'α', beta: 'β', gamma: 'γ', delta: 'δ', epsilon: 'ε', zeta: 'ζ', eta: 'η',
  theta: 'θ', iota: 'ι', kappa: 'κ', lambda: 'λ', mu: 'μ', nu: 'ν', xi: 'ξ', rho: 'ρ',
  sigma: 'σ', tau: 'τ', upsilon: 'υ', phi: 'φ', chi: 'χ', psi: 'ψ', omega: 'ω',
  Gamma: 'Γ', Delta: 'Δ', Theta: 'Θ', Lambda: 'Λ', Sigma: 'Σ', Phi: 'Φ', Psi: 'Ψ', Omega: 'Ω',
};
var GREEK_EXACT = { lam: 'λ', pi: 'π' };
var GREEK_KEYS = Object.keys(GREEK).sort(function (a, b) { return b.length - a.length || (a < b ? -1 : 1); });

/** Default display symbol for an identifier, e.g. "gamma" -> "γ", "tau_A" -> "τ_A", "R0" -> "R_0". */
function autoSymbol(name) {
  if (GREEK_EXACT[name]) return GREEK_EXACT[name];
  if (GREEK[name]) return GREEK[name];
  for (var i = 0; i < GREEK_KEYS.length; i++) {
    var k = GREEK_KEYS[i];
    if (name.indexOf(k) === 0) {
      var underscore = name[k.length] === '_';
      var rest = name.slice(k.length + (underscore ? 1 : 0));
      if (/^[A-Za-z0-9]+$/.test(rest) && (underscore || rest.length <= 2)) return GREEK[k] + '_' + rest;
    }
  }
  var m = /^([A-Za-z]+)(\d+)$/.exec(name);
  if (m && m[1].length <= 2) return m[1] + '_' + m[2];
  return name;
}

/**
 * Turn an expression into pretty display segments [{text, sub}] using the
 * given symbol map (identifier -> symbol, "x_y" renders y as a subscript).
 */
function prettySegments(src, symbols) {
  var toks;
  try { toks = tokenize(src); } catch (e) { return [{ text: String(src), sub: false }]; }
  var segs = [];
  function push(text, sub) {
    if (!text) return;
    var last = segs[segs.length - 1];
    if (last && last.sub === !!sub) last.text += text; else segs.push({ text: text, sub: !!sub });
  }
  function pushSymbol(sym) {
    var u = sym.indexOf('_');
    if (u > 0) { push(sym.slice(0, u), false); push(sym.slice(u + 1), true); } else push(sym, false);
  }
  var prevOperandEnd = false;
  for (var i = 0; i < toks.length; i++) {
    var tk = toks[i];
    var next = toks[i + 1];
    if (tk.t === 'id') {
      if (next && next.t === 'lp') { push(tk.v, false); prevOperandEnd = false; continue; }
      pushSymbol((symbols && symbols[tk.v]) || autoSymbol(tk.v));
      prevOperandEnd = true;
    } else if (tk.t === 'num') {
      push(tk.v, false); prevOperandEnd = true;
    } else if (tk.t === 'lp') {
      push('(', false); prevOperandEnd = false;
    } else if (tk.t === 'rp') {
      push(')', false); prevOperandEnd = true;
    } else if (tk.t === 'comma') {
      push(', ', false); prevOperandEnd = false;
    } else if (tk.t === 'op') {
      var binary = prevOperandEnd;
      if (tk.v === '*' && binary) push(next && next.t === 'num' ? '·' : '', false);
      else if (tk.v === '-') push('−', false);
      else if (tk.v === '>=') push(' ≥ ', false);
      else if (tk.v === '<=') push(' ≤ ', false);
      else if (tk.v === '!=') push(' ≠ ', false);
      else if (tk.v === '==') push(' = ', false);
      else if (binary && (tk.v === '<' || tk.v === '>' || tk.v === '&' || tk.v === '|' || tk.v === '&&' || tk.v === '||')) push(' ' + tk.v + ' ', false);
      else push(tk.v, false);
      prevOperandEnd = false;
    }
  }
  return segs;
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

function compartmentById(model) {
  var map = {};
  (model.compartments || []).forEach(function (c) { map[c.id] = c; });
  return map;
}

function flowTitle(f, byId) {
  var a = f.from && byId[f.from] ? byId[f.from].name : '∅';
  var b = f.to && byId[f.to] ? byId[f.to].name : '∅';
  return 'Flow ' + a + ' → ' + b;
}

/** Stable topological sort. items: [{name, deps:[names]}]; returns {order, cycle} */
function stableTopo(items) {
  var names = items.map(function (it) { return it.name; });
  var done = {};
  var order = [];
  var remaining = items.slice();
  while (remaining.length) {
    var idx = -1;
    for (var i = 0; i < remaining.length; i++) {
      var it = remaining[i];
      var ready = it.deps.every(function (d) { return d !== it.name && (done[d] || !has(names, d)); });
      if (ready) { idx = i; break; }
    }
    if (idx === -1) return { order: order, cycle: remaining.map(function (r) { return r.name; }) };
    done[remaining[idx].name] = true;
    order.push(remaining[idx]);
    remaining.splice(idx, 1);
  }
  return { order: order, cycle: null };
}

/**
 * Validate a model. Returns { errors: [...], warnings: [...], info } where each
 * issue is { kind, id, field, message, unknown? } and info holds the analysed
 * pieces needed by generate().
 */
function validate(model) {
  var errors = [];
  var warnings = [];
  function err(kind, id, field, message, extra) {
    var e = { kind: kind, id: id, field: field, message: message };
    if (extra) for (var k in extra) e[k] = extra[k];
    errors.push(e);
  }

  model = model || {};
  var comps = model.compartments || [];
  var params = model.parameters || [];
  var vars = model.variables || [];
  var flows = model.flows || [];
  var time = model.time || {};
  var popName = (model.population && model.population.name) || 'N';
  var byId = compartmentById(model);

  if (!comps.length) err('model', null, 'compartments', 'Add at least one compartment');

  // --- Names -------------------------------------------------------------
  var owners = {}; // name -> description of what defines it
  function claim(name, kind, id, field, what) {
    if (!IDENT_RE.test(name || '')) {
      err(kind, id, field, what + ' name "' + (name || '') + '" is not a valid R name (start with a letter; letters, digits, "_" or "." only)');
      return;
    }
    if (has(R_RESERVED, name) || has(INTERNAL, name) || has(R_FUNCS, name) || has(R_CONSTS, name)) {
      err(kind, id, field, '"' + name + '" is reserved in R or by the generated script; choose another name');
      return;
    }
    if (owners[name]) {
      err(kind, id, field, '"' + name + '" is already used by ' + owners[name]);
      return;
    }
    owners[name] = what.toLowerCase() + ' ' + name;
  }

  claim(popName, 'model', null, 'population', 'Population total');
  comps.forEach(function (c) { claim(c.name, 'compartment', c.id, 'name', 'Compartment'); });
  comps.forEach(function (c) {
    if (IDENT_RE.test(c.name || '') && owners['d' + c.name] === undefined) owners['d' + c.name] = 'the derivative of ' + c.name;
  });
  params.forEach(function (p, i) { claim(p.name, 'parameter', i, 'name', 'Parameter'); });
  vars.forEach(function (v, i) { claim(v.name, 'variable', i, 'name', 'Variable'); });
  flows.forEach(function (f) {
    if (f.output) claim(f.output, 'flow', f.id, 'output', 'Flow output');
  });

  var compNames = comps.map(function (c) { return c.name; });
  var paramNames = params.map(function (p) { return p.name; });
  var varNames = vars.map(function (v) { return v.name; });
  var inN = comps.filter(function (c) { return c.inN !== false; });

  function checkExpr(src, allowed, kind, id, field, label) {
    var a = analyze(src);
    if (!a.ok) { err(kind, id, field, label + ': ' + a.error); return a; }
    a.funcs.forEach(function (fn) {
      if (!has(R_FUNCS, fn)) err(kind, id, field, label + ': function ' + fn + '() is not supported (allowed: ' + R_FUNCS.join(', ') + ')');
    });
    a.idents.forEach(function (name) {
      if (has(R_CONSTS, name)) return;
      if (has(allowed, name)) return;
      var why = has(compNames, name) || has(varNames, name) || name === popName || name === 't'
        ? '"' + name + '" cannot be used here'
        : 'unknown name "' + name + '"';
      err(kind, id, field, label + ': ' + why, why.indexOf('unknown') === 0 ? { unknown: name } : null);
    });
    return a;
  }

  // --- Time ----------------------------------------------------------------
  var timeA = {};
  ['start', 'stop', 'step'].forEach(function (k) {
    timeA[k] = checkExpr(time[k], [], 'model', null, 'time.' + k, 'Time ' + k);
  });

  // --- Solver options --------------------------------------------------------
  var hmaxA = null;
  if (model.solver && model.solver.hmax != null && String(model.solver.hmax).trim() !== '') {
    hmaxA = checkExpr(model.solver.hmax, [], 'model', null, 'solver.hmax', 'Maximum solver step');
  }

  // --- Parameters (constants: numbers / arithmetic only) ------------------
  var paramA = params.map(function (p, i) {
    return checkExpr(p.value, [], 'parameter', i, 'value', 'Parameter ' + (p.name || '#' + (i + 1)) + ' value');
  });

  // --- Initial conditions ----------------------------------------------------
  var initA = comps.map(function (c) {
    return checkExpr(c.init, paramNames.concat(compNames), 'compartment', c.id, 'init', 'Initial value of ' + c.name);
  });
  var initTopo = stableTopo(comps.map(function (c, i) {
    return { name: c.name, index: i, deps: initA[i].ok ? initA[i].idents.filter(function (n) { return has(compNames, n); }) : [] };
  }));
  if (initTopo.cycle) err('model', null, 'init', 'Initial values depend on each other in a cycle: ' + initTopo.cycle.join(', '));

  // --- Derived variables -----------------------------------------------------
  var dynamicNames = compNames.concat(paramNames, varNames, ['t']);
  if (inN.length) dynamicNames.push(popName);
  var varA = vars.map(function (v, i) {
    return checkExpr(v.expr, dynamicNames, 'variable', i, 'expr', 'Variable ' + (v.name || '#' + (i + 1)));
  });
  var varTopo = stableTopo(vars.map(function (v, i) {
    return { name: v.name, index: i, deps: varA[i].ok ? varA[i].idents.filter(function (n) { return has(varNames, n); }) : [] };
  }));
  if (varTopo.cycle) err('model', null, 'variables', 'Variables depend on each other in a cycle: ' + varTopo.cycle.join(', '));

  // --- Flows -----------------------------------------------------------------
  var flowA = flows.map(function (f) {
    var title = flowTitle(f, byId);
    if (f.from != null && !byId[f.from]) err('flow', f.id, 'from', title + ': source compartment no longer exists');
    if (f.to != null && !byId[f.to]) err('flow', f.id, 'to', title + ': target compartment no longer exists');
    if (f.from == null && f.to == null) err('flow', f.id, 'to', title + ': a flow needs at least one compartment');
    if (f.from != null && f.from === f.to) err('flow', f.id, 'to', title + ': source and target are the same');
    var kind = f.kind === 'absolute' ? 'absolute' : 'per_capita';
    if (kind === 'per_capita' && f.from == null) {
      err('flow', f.id, 'kind', title + ': a per-capita rate needs a source compartment (use "absolute" for inflows such as births)');
    }
    return checkExpr(f.rate, dynamicNames, 'flow', f.id, 'rate', title + ' rate');
  });

  // --- Warnings --------------------------------------------------------------
  comps.forEach(function (c) {
    var touched = flows.some(function (f) { return f.from === c.id || f.to === c.id; });
    if (!touched && comps.length > 1) warnings.push({ kind: 'compartment', id: c.id, field: 'name', message: 'Compartment ' + c.name + ' has no flows' });
  });
  var used = {};
  flowA.concat(varA, initA).forEach(function (a) { (a.idents || []).forEach(function (n) { used[n] = true; }); });
  params.forEach(function (p, i) {
    if (p.name && !used[p.name]) warnings.push({ kind: 'parameter', id: i, field: 'name', message: 'Parameter ' + p.name + ' is not used' });
  });

  return {
    errors: errors,
    warnings: warnings,
    info: {
      popName: popName, byId: byId, timeA: timeA, hmaxA: hmaxA, paramA: paramA, initA: initA, initOrder: initTopo.order,
      varA: varA, varOrder: varTopo.order, flowA: flowA, inN: inN,
    },
  };
}

// ---------------------------------------------------------------------------
// Code generation
// ---------------------------------------------------------------------------

function rString(s) {
  return '"' + String(s).replace(/\\/g, '\\\\').replace(/"/g, '\\"') + '"';
}

function oneLine(s) { return String(s || '').replace(/\s+/g, ' ').trim(); }

function padRight(s, n) { while (s.length < n) s += ' '; return s; }

function withComments(lines, indent) {
  // lines: [{code, comment}] -> aligned comment column
  var width = 0;
  lines.forEach(function (l) { if (l.comment) width = Math.max(width, l.code.length); });
  return lines.map(function (l) {
    return indent + (l.comment ? padRight(l.code, width) + '  # ' + l.comment : l.code);
  });
}

function functionName(name) {
  var base = String(name || 'model').replace(/[^A-Za-z0-9_]+/g, '_').replace(/^_+|_+$/g, '');
  if (!base) base = 'model';
  if (!/^[A-Za-z]/.test(base)) base = 'm_' + base;
  return base + '_model';
}

function lineColour(hex) {
  var m = /^#?([0-9a-fA-F]{6})$/.exec(hex || '');
  if (!m) return '#333333';
  var n = parseInt(m[1], 16);
  var r = (n >> 16) & 255, g = (n >> 8) & 255, b = n & 255;
  var lum = (0.2126 * r + 0.7152 * g + 0.0722 * b) / 255;
  return lum > 0.85 ? '#333333' : '#' + m[1].toUpperCase();
}

/** Build the R term for a flow (without sign). */
function flowTerm(f, a, byId) {
  if (f.kind === 'absolute') return a;
  return a + ' * ' + byId[f.from].name;
}

/**
 * Generate R code. Returns { ok, code, errors, warnings }.
 * code is null when the model has errors.
 */
function generate(model) {
  var v = validate(model);
  if (v.errors.length) return { ok: false, code: null, errors: v.errors, warnings: v.warnings };
  var I = v.info;
  var comps = model.compartments;
  var params = model.parameters || [];
  var vars = model.variables || [];
  var flows = model.flows || [];
  var time = model.time;
  var unit = oneLine(time.unit) || 'time';
  var fname = functionName(model.id && /^[A-Za-z][A-Za-z0-9_]*$/.test(model.id) ? model.id : model.name);
  var title = oneLine(model.name) || 'Compartmental model';
  var L = [];

  // Per-flow R terms
  var terms = flows.map(function (f, i) {
    var a = I.flowA[i];
    if (f.kind === 'absolute') return a.lowPrecedence ? '(' + a.code + ')' : a.code;
    return flowTerm(f, paren(a), I.byId);
  });

  // --- Header ---------------------------------------------------------------
  L.push('## ' + title);
  if (model.description) {
    String(model.description).split(/\r?\n/).forEach(function (line) { L.push('## ' + line.replace(/\s+$/, '')); });
  }
  if (model.reference) L.push('## Reference: ' + oneLine(model.reference));
  if (model.referenceUrl) L.push('## Source: ' + oneLine(model.referenceUrl));
  L.push('##');
  L.push('## Generated from a compartmental model diagram (model format v' + FORMAT_VERSION + ').');
  L.push('## Regenerating from the same diagram always gives this exact file.');
  L.push('##');
  L.push('## Compartments: ' + comps.map(function (c) { return c.name; }).join(', '));
  L.push('## Flows (time unit: ' + unit + '):');
  var flowRows = flows.map(function (f, i) {
    var a = f.from != null ? I.byId[f.from].name : '(outside)';
    var b = f.to != null ? I.byId[f.to].name : '(outside)';
    return { left: a + ' -> ' + b, right: terms[i] };
  });
  var fw = 0;
  flowRows.forEach(function (r) { fw = Math.max(fw, r.left.length); });
  flowRows.forEach(function (r) { L.push('##   ' + padRight(r.left, fw) + '  : ' + r.right); });
  if (!flows.length) L.push('##   (none)');
  L.push('');
  L.push('library(deSolve)');
  L.push('');

  // --- Time -----------------------------------------------------------------
  L.push('# ---- Time ------------------------------------------------------------');
  L.push('times <- seq(' + I.timeA.start.code + ', ' + I.timeA.stop.code + ', by = ' + I.timeA.step.code + ')  # ' + unit);
  L.push('');

  // --- Parameters -------------------------------------------------------------
  L.push('# ---- Parameters ------------------------------------------------------');
  if (!params.length) {
    L.push('parameters <- c()');
  } else {
    L.push('parameters <- c(');
    var plines = params.map(function (p, i) {
      return { code: p.name + ' = ' + I.paramA[i].code + (i < params.length - 1 ? ',' : ''), comment: oneLine(p.description) };
    });
    withComments(plines, '  ').forEach(function (l) { L.push(l); });
    L.push(')');
  }
  L.push('');

  // --- Initial conditions -----------------------------------------------------
  L.push('# ---- Initial conditions ----------------------------------------------');
  var initNeedsEnv = I.initA.some(function (a) { return a.idents.length > 0; });
  var stateVec = 'c(' + comps.map(function (c) { return c.name + ' = ' + c.name; }).join(', ') + ')';
  if (!initNeedsEnv) {
    var ilines = comps.map(function (c, i) {
      return { code: c.name + ' = ' + I.initA[i].code + (i < comps.length - 1 ? ',' : ''), comment: oneLine(c.description) };
    });
    L.push('state <- c(');
    withComments(ilines, '  ').forEach(function (l) { L.push(l); });
    L.push(')');
  } else {
    L.push('state <- with(as.list(parameters), {');
    var olines = I.initOrder.map(function (it) {
      var c = comps[it.index];
      return { code: c.name + ' <- ' + I.initA[it.index].code, comment: oneLine(c.description) };
    });
    withComments(olines, '  ').forEach(function (l) { L.push(l); });
    L.push('  ' + stateVec);
    L.push('})');
  }
  L.push('');

  // --- Model function -----------------------------------------------------------
  L.push('# ---- Model equations -------------------------------------------------');
  L.push(fname + ' <- function(t, state, parameters) {');
  L.push('  with(as.list(c(state, parameters)), {');
  var body = [];
  if (I.inN.length) {
    body.push('# Total population');
    body.push(I.popName + ' <- ' + I.inN.map(function (c) { return c.name; }).join(' + '));
    body.push('');
  }
  if (vars.length) {
    body.push('# Derived quantities');
    withComments(I.varOrder.map(function (it) {
      var vv = vars[it.index];
      return { code: vv.name + ' <- ' + I.varA[it.index].code, comment: oneLine(vv.description) };
    }), '').forEach(function (l) { body.push(l); });
    body.push('');
  }
  var flowOutputs = flows.map(function (f, i) { return { f: f, i: i }; }).filter(function (x) { return x.f.output; });
  if (flowOutputs.length) {
    body.push('# Flows recorded as outputs');
    flowOutputs.forEach(function (x) { body.push(x.f.output + ' <- ' + terms[x.i]); });
    body.push('');
  }
  body.push('# Rates of change');
  comps.forEach(function (c) {
    var expr = '';
    flows.forEach(function (f, i) {
      var sign = 0;
      if (f.from === c.id) sign = -1;
      else if (f.to === c.id) sign = 1;
      if (!sign) return;
      if (!expr) expr = (sign < 0 ? '-' : '') + terms[i];
      else expr += (sign < 0 ? ' - ' : ' + ') + terms[i];
    });
    body.push('d' + c.name + ' <- ' + (expr || '0'));
  });
  body.push('');
  var outNames = vars.filter(function (vv) { return vv.output; }).map(function (vv) { return vv.name; })
    .concat(flowOutputs.map(function (x) { return x.f.output; }));
  var derivs = 'c(' + comps.map(function (c) { return 'd' + c.name; }).join(', ') + ')';
  if (outNames.length) {
    body.push('list(' + derivs + ', c(' + outNames.map(function (n) { return n + ' = ' + n; }).join(', ') + '))');
  } else {
    body.push('list(' + derivs + ')');
  }
  body.forEach(function (l) { L.push(l ? '    ' + l : ''); });
  L.push('  })');
  L.push('}');
  L.push('');

  // --- Solve ------------------------------------------------------------------------
  L.push('# ---- Solve -----------------------------------------------------------');
  L.push('out <- ode(y = state, times = times, func = ' + fname + ', parms = parameters'
    + (I.hmaxA ? ', hmax = ' + I.hmaxA.code : '') + ')');
  L.push('out <- as.data.frame(out)');
  L.push('');

  // --- Plot -----------------------------------------------------------------------
  L.push('# ---- Plot ------------------------------------------------------------');
  L.push('compartment_cols <- c(' + comps.map(function (c) { return c.name + ' = ' + rString(lineColour(c.color)); }).join(', ') + ')');
  L.push('matplot(out$time, out[, names(compartment_cols)], type = "l", lty = 1, lwd = 2,');
  L.push('        col = compartment_cols, xlab = ' + rString('Time (' + unit + ')') + ', ylab = "Number of individuals",');
  L.push('        main = ' + rString(title) + ')');
  L.push('legend("topright", legend = names(compartment_cols), col = compartment_cols,');
  L.push('       lty = 1, lwd = 2, bty = "n")');
  if (outNames.length) {
    L.push('');
    L.push('output_names <- c(' + outNames.map(rString).join(', ') + ')');
    L.push('matplot(out$time, out[, output_names, drop = FALSE], type = "l", lty = 1, lwd = 2,');
    L.push('        col = seq_along(output_names), xlab = ' + rString('Time (' + unit + ')') + ', ylab = "Value",');
    L.push('        main = ' + rString(title + ': outputs') + ')');
    L.push('legend("topright", legend = output_names, col = seq_along(output_names),');
    L.push('       lty = 1, lwd = 2, bty = "n")');
  }
  L.push('');

  return { ok: true, code: L.join('\n'), errors: [], warnings: v.warnings };
}

// ---------------------------------------------------------------------------
// Model utilities
// ---------------------------------------------------------------------------

/** Rename an identifier everywhere it is referenced (expressions only). Mutates model. */
function renameReferences(model, from, to) {
  if (!from || from === to) return;
  (model.compartments || []).forEach(function (c) { c.init = renameInExpr(c.init, from, to); });
  (model.parameters || []).forEach(function (p) { p.value = renameInExpr(p.value, from, to); });
  (model.variables || []).forEach(function (v) { v.expr = renameInExpr(v.expr, from, to); });
  (model.flows || []).forEach(function (f) { f.rate = renameInExpr(f.rate, from, to); });
}

/** Canonical JSON for saving (fixed key order => stable diffs). */
function serialize(model) {
  function pick(o, keys) {
    var r = {};
    keys.forEach(function (k) { if (o[k] !== undefined) r[k] = o[k]; });
    return r;
  }
  var m = {
    version: FORMAT_VERSION,
    id: model.id,
    category: model.category,
    name: model.name || '',
    description: model.description || '',
    reference: model.reference || '',
    referenceUrl: model.referenceUrl || '',
    time: pick(model.time || {}, ['start', 'stop', 'step', 'unit']),
    population: pick(model.population || { name: 'N' }, ['name']),
    solver: model.solver && model.solver.hmax ? pick(model.solver, ['hmax']) : undefined,
    compartments: (model.compartments || []).map(function (c) {
      return pick(c, ['id', 'name', 'description', 'init', 'inN', 'color', 'x', 'y']);
    }),
    parameters: (model.parameters || []).map(function (p) { return pick(p, ['name', 'symbol', 'value', 'description']); }),
    variables: (model.variables || []).map(function (v) { return pick(v, ['name', 'symbol', 'expr', 'description', 'output']); }),
    flows: (model.flows || []).map(function (f) {
      return pick(f, ['id', 'from', 'to', 'kind', 'rate', 'output', 'curve', 'ox', 'oy']);
    }),
  };
  return JSON.stringify(m, null, 2) + '\n';
}

export {
  FORMAT_VERSION,
  R_FUNCS,
  tokenize,
  analyze,
  validate,
  generate,
  autoSymbol,
  prettySegments,
  renameInExpr,
  renameReferences,
  serialize,
  functionName,
};
