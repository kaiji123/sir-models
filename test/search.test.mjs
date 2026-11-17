import test from 'node:test';
import assert from 'node:assert/strict';
import { indexModel, search, parseQuery } from '../src/lib/search.js';
import { LIBRARY } from '../src/lib/library.js';

const docs = LIBRARY.map(indexModel);
const ids = (q, opts) => search(docs, q, opts).map((h) => h.model.id);

test('empty query returns the whole library in order', () => {
  assert.deepEqual(ids(''), LIBRARY.map((m) => m.id));
});

test('name matches rank first; all words must match', () => {
  assert.equal(ids('seir')[0], 'seir');
  assert.ok(ids('malaria').includes('ross_macdonald'));
  assert.deepEqual(ids('malaria zzzz'), []);
});

test('synonyms, prefixes and accents', () => {
  assert.ok(ids('tb').includes('tb'));
  assert.ok(ids('vacc').includes('leaky_vaccine'));
  assert.ok(ids('Ross–Macdonald').includes('ross_macdonald'));
});

test('operators: has:, cat:, -exclude, "phrase"', () => {
  assert.ok(ids('has:W').every((id) => LIBRARY.find((m) => m.id === id).compartments.some((c) => c.name === 'W')));
  assert.ok(ids('has:W').length >= 2);
  assert.ok(ids('cat:interv').every((id) => LIBRARY.find((m) => m.id === id).category === 'Interventions'));
  assert.ok(!ids('vaccin -leaky').includes('leaky_vaccine'));
  assert.ok(ids('"waning immunity"').length > 0);
  assert.deepEqual(parseQuery('-a "b c" has:E').map((t) => t.kind), ['word', 'phrase', 'has']);
});

test('category filter', () => {
  const r = ids('', { category: 'Basic' });
  assert.ok(r.length > 0 && r.every((id) => LIBRARY.find((m) => m.id === id).category === 'Basic'));
});
