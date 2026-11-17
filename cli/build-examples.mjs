#!/usr/bin/env node
// Writes examples/<id>.json and examples/<id>.R for every model in the library.
// These files are the golden snapshots checked by the tests.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { generate, serialize } from '../src/lib/generator.js';
import { LIBRARY } from '../src/lib/library.js';

const dir = path.join(path.dirname(fileURLToPath(import.meta.url)), '..', 'examples');
fs.mkdirSync(dir, { recursive: true });
let failed = 0;
for (const m of LIBRARY) {
  const res = generate(m);
  if (!res.ok) {
    failed++;
    console.error(`${m.id}: ${res.errors.map((e) => e.message).join('; ')}`);
    continue;
  }
  res.warnings.forEach((w) => console.warn(`${m.id}: warning: ${w.message}`));
  fs.writeFileSync(path.join(dir, m.id + '.json'), serialize(m));
  fs.writeFileSync(path.join(dir, m.id + '.R'), res.code);
}
console.log(`${LIBRARY.length - failed}/${LIBRARY.length} models written to examples/`);
process.exit(failed ? 1 : 0);
