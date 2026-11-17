import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { generate, validate, serialize, analyze, renameReferences, autoSymbol } from '../src/lib/generator.js';
import { LIBRARY, BY_ID } from '../src/lib/library.js';

const root = path.join(path.dirname(fileURLToPath(import.meta.url)), '..');
const clone = (o) => JSON.parse(JSON.stringify(o));
// Deep copy with every object's keys in reverse order.
const reverseKeys = (o) => (Array.isArray(o) ? o.map(reverseKeys)
  : o && typeof o === 'object' ? Object.fromEntries(Object.entries(o).reverse().map(([k, v]) => [k, reverseKeys(v)])) : o);

test('every library model is valid and matches its golden R file (run `npm run examples` to update)', () => {
  for (const m of LIBRARY) {
    const res = generate(m);
    assert.ok(res.ok, `${m.id}: ${res.errors.map((e) => e.message).join('; ')}`);
    const golden = fs.readFileSync(path.join(root, 'examples', m.id + '.R'), 'utf8');
    assert.equal(res.code, golden, `${m.id}.R is out of date`);
  }
});

test('generation is deterministic: repeated runs, JSON round-trip and key order do not matter', () => {
  for (const m of LIBRARY) {
    const a = generate(m).code;
    assert.equal(generate(m).code, a);
    assert.equal(generate(JSON.parse(serialize(m))).code, a);
    assert.equal(generate(reverseKeys(m)).code, a, m.id);
  }
});

test('diagram layout does not affect the code', () => {
  const m = clone(BY_ID.seir);
  const before = generate(m).code;
  m.compartments.forEach((c) => { c.x += 123; c.y -= 45; });
  m.flows.forEach((f) => { f.curve = 77; });
  assert.equal(generate(m).code, before);
});

test('flows become correctly signed terms (per-capita x source, absolute as-is)', () => {
  const code = generate(BY_ID.hepE).code;
  assert.match(code, /dS <- mui \* P - muo \* S - lam \* S \+ omega \* R\n/);
  assert.match(code, /dE <- lam \* S - muo \* E - gamma \* E\n/);
  assert.match(code, /list\(c\(dS, dE, dI, dR\), c\(inc = inc\)\)/);
  const seachr = generate(BY_ID.seachr).code;
  assert.match(seachr, /dA <- nu \* \(1 - gamma\) \* E - tauA \* A \+ nu \* \(1 - gamma \/ 2\) \* ER\n/);
  assert.match(seachr, /dD <- tauH \* omega \* H\n/);
});

test('rates with low-precedence operators are parenthesised', () => {
  const m = clone(BY_ID.sir);
  m.flows[1].rate = 'gamma + 0.01';
  assert.match(generate(m).code, /dI <- lam \* S - \(gamma \+ 0\.01\) \* I/);
  m.flows[1].rate = 'gamma * (1 - 0.1)';
  assert.match(generate(m).code, /dI <- lam \* S - gamma \* \(1 - 0\.1\) \* I/);
});

test('variables are emitted in dependency order', () => {
  const m = clone(BY_ID.seir_vital);
  m.variables.reverse(); // lam listed before beta
  const code = generate(m).code;
  assert.ok(code.indexOf('beta <- ') < code.indexOf('lam <- '));
});

test('initial values may depend on parameters and each other', () => {
  const m = clone(BY_ID.sir);
  m.parameters.push({ name: 'N0', value: '1000', description: '' });
  m.compartments[0].init = 'N0 - I';
  const code = generate(m).code;
  assert.match(code, /state <- with\(as\.list\(parameters\), \{\n {2}I <- 10[^\n]*\n {2}S <- N0 - I[^\n]*\n {2}R <- 0/);
});

function errorsOf(m) { return validate(m).errors.map((e) => e.message); }

test('validation: unknown names, bad syntax, reserved names, duplicates', () => {
  let m = clone(BY_ID.sir);
  m.variables[0].expr = 'betta * I / N';
  assert.ok(validate(m).errors.some((e) => e.unknown === 'betta'));

  m = clone(BY_ID.sir);
  m.flows[0].rate = 'beta * (I';
  assert.match(errorsOf(m).join(), /Unclosed/);

  m = clone(BY_ID.sir);
  m.compartments[2].name = 'F';
  assert.match(errorsOf(m).join(), /reserved/);

  m = clone(BY_ID.sir);
  m.parameters[1].name = 'beta';
  assert.match(errorsOf(m).join(), /already used/);

  m = clone(BY_ID.sir);
  m.parameters[0].value = 'gamma * 3';
  assert.match(errorsOf(m).join(), /unknown name "gamma"|cannot be used/);

  m = clone(BY_ID.sir);
  m.variables[0].expr = 'foo(I)';
  assert.match(errorsOf(m).join(), /not supported/);
});

test('validation: per-capita inflow, self-loops and cycles', () => {
  let m = clone(BY_ID.sir);
  m.flows.push({ id: 'fx', from: null, to: 'c1', kind: 'per_capita', rate: 'beta', output: '', curve: 0 });
  assert.match(errorsOf(m).join(), /per-capita rate needs a source/);

  m = clone(BY_ID.sir);
  m.flows[0].to = m.flows[0].from;
  assert.match(errorsOf(m).join(), /same/);

  m = clone(BY_ID.sir);
  m.variables.push({ name: 'a', expr: 'b', description: '' }, { name: 'b', expr: 'a', description: '' });
  assert.match(errorsOf(m).join(), /cycle/);
  assert.equal(generate(m).code, null);
});

test('renameReferences rewrites whole identifiers only', () => {
  const m = clone(BY_ID.sir);
  m.variables[0].expr = 'beta * I / N + beta2';
  renameReferences(m, 'beta', 'b');
  assert.equal(m.variables[0].expr, 'b * I / N + beta2');
});

test('expression analysis and display symbols', () => {
  assert.equal(analyze('nu*(1-gamma)').code, 'nu * (1 - gamma)');
  assert.equal(analyze('-a^2').lowPrecedence, true);
  assert.equal(analyze('a/b*c').lowPrecedence, false);
  assert.equal(analyze('(t>=t0)*c0').code, '(t >= t0) * c0');
  assert.equal(autoSymbol('gamma'), 'γ');
  assert.equal(autoSymbol('tau_A'), 'τ_A');
  assert.equal(autoSymbol('muo'), 'μ_o');
  assert.equal(autoSymbol('R0'), 'R_0');
  assert.equal(autoSymbol('lam'), 'λ');
});
