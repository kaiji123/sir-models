#!/usr/bin/env node
// Usage: node cli/generate.mjs model.json [out.R]
// Prints errors and exits with status 1 if the model is invalid.
import fs from 'node:fs';
import { generate } from '../src/lib/generator.js';

const [input, output] = process.argv.slice(2);
if (!input) {
  console.error('Usage: node cli/generate.mjs <model.json> [out.R]');
  process.exit(2);
}
const res = generate(JSON.parse(fs.readFileSync(input, 'utf8')));
res.warnings.forEach((w) => console.error('warning: ' + w.message));
if (!res.ok) {
  res.errors.forEach((e) => console.error('error: ' + e.message));
  process.exit(1);
}
if (output) fs.writeFileSync(output, res.code);
else process.stdout.write(res.code);
