// Library search: pure functions shared by the editor picker, the /library/ page and the tests.
//
// Query syntax
//   words            all words must match (prefix match, case/accent-insensitive)
//   "exact phrase"   phrase must appear
//   -word            exclude models matching word
//   has:E            model has a compartment named E (exact, case-sensitive)
//   cat:vector       category starts with "vector"
// Common abbreviations are expanded (flu -> influenza, tb -> tuberculosis, ...).

const SYNONYMS = {
  flu: ['influenza'],
  influenza: ['flu'],
  covid: ['covid', 'sars-cov-2', 'coronavirus'],
  corona: ['covid', 'coronavirus'],
  tb: ['tuberculosis'],
  tuberculosis: ['tb'],
  std: ['sexually', 'sti'],
  sti: ['sexually', 'sti'],
  mosquito: ['vector', 'mosquito', 'mosquitoes'],
  vaccine: ['vaccin'],
  vaccination: ['vaccin'],
  hiv: ['hiv', 'aids'],
  hep: ['hepatitis'],
  hbv: ['hepatitis b'],
  hcv: ['hepatitis c'],
  rsv: ['rsv', 'respiratory syncytial'],
  amr: ['resistance', 'resistant'],
  zoonotic: ['zoonos', 'zoonotic', 'spillover'],
  livestock: ['livestock', 'cattle', 'pig', 'poultry', 'animal'],
  water: ['water', 'environment'],
  waterborne: ['water'],
  ntd: ['neglected', 'helminth', 'ntd'],
  withinhost: ['within-host', 'within host'],
};

const FIELD_WEIGHTS = { name: 10, id: 8, compartments: 6, category: 4, keywords: 4, description: 2, reference: 2, details: 1 };

export function fold(s) {
  return String(s || '').normalize('NFD').replace(/[̀-ͯ]/g, '').toLowerCase();
}

/** Build a search document for a model (cheap; can be cached). */
export function indexModel(m) {
  return {
    model: m,
    compNames: new Set(m.compartments.map((c) => c.name)),
    fields: {
      name: fold(m.name),
      id: fold(m.id),
      compartments: fold(m.compartments.map((c) => c.name + ' ' + c.description).join(' ')),
      category: fold(m.category),
      keywords: fold((m.keywords || []).join(' ')),
      description: fold(m.description),
      reference: fold(m.reference),
      details: fold(
        m.parameters.map((p) => p.name + ' ' + p.description).join(' ') + ' '
        + m.variables.map((v) => v.name + ' ' + v.description).join(' '),
      ),
    },
  };
}

export function parseQuery(q) {
  const terms = [];
  const re = /(-?)(?:"([^"]*)"|(\S+))/g;
  let m;
  while ((m = re.exec(String(q || '')))) {
    const neg = m[1] === '-';
    const raw = m[2] != null ? m[2] : m[3];
    if (!raw) continue;
    const lower = raw.toLowerCase();
    if (lower.startsWith('has:') && raw.length > 4) terms.push({ kind: 'has', value: raw.slice(4), neg });
    else if (lower.startsWith('cat:') && raw.length > 4) terms.push({ kind: 'cat', value: fold(raw.slice(4)), neg });
    else terms.push({ kind: m[2] != null ? 'phrase' : 'word', value: fold(raw), neg });
  }
  return terms;
}

function wordScore(doc, word, phrase) {
  const alts = phrase ? [word] : [word, ...(SYNONYMS[word] || [])];
  let best = 0;
  for (const [field, text] of Object.entries(doc.fields)) {
    for (const alt of alts) {
      if (!alt) continue;
      let s = 0;
      if (phrase) {
        if (text.includes(alt)) s = FIELD_WEIGHTS[field];
      } else {
        // prefix of a word (start of string or after a non-alphanumeric character)
        const re = new RegExp('(^|[^a-z0-9])' + alt.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'));
        if (re.test(text)) s = FIELD_WEIGHTS[field] * (alt === word ? 1 : 0.8);
        else if (alt.length >= 4 && text.includes(alt)) s = FIELD_WEIGHTS[field] * 0.5;
      }
      if (s > best) best = s;
    }
  }
  return best;
}

/**
 * Search indexed docs. Returns [{ model, score }] sorted by score (desc) then
 * library order. An empty query returns everything in library order.
 */
export function search(docs, query, { category = null } = {}) {
  const terms = parseQuery(query);
  const out = [];
  docs.forEach((doc, order) => {
    if (category && doc.model.category !== category) return;
    let score = 0;
    for (const t of terms) {
      let s;
      if (t.kind === 'has') s = doc.compNames.has(t.value) ? 5 : 0;
      else if (t.kind === 'cat') s = doc.fields.category.startsWith(t.value) ? 4 : 0;
      else s = wordScore(doc, t.value, t.kind === 'phrase');
      if (t.neg) { if (s > 0) return; } else if (s === 0) return; else score += s;
    }
    out.push({ model: doc.model, score, order });
  });
  out.sort((a, b) => b.score - a.score || a.order - b.order);
  return out.map(({ model, score }) => ({ model, score }));
}
