// Runs the generated code in R. Skipped when Rscript or deSolve is unavailable.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { generate } from '../src/lib/generator.js';
import { LIBRARY } from '../src/lib/library.js';

const here = path.dirname(fileURLToPath(import.meta.url));
const probe = spawnSync('Rscript', ['-e', 'cat(requireNamespace("deSolve", quietly = TRUE))'], { encoding: 'utf8' });
const skip = probe.status !== 0 || !probe.stdout.includes('TRUE') ? 'Rscript with deSolve not available' : false;

function tmpDir() { return fs.mkdtempSync(path.join(os.tmpdir(), 'sir-r-')); }

test('every library model runs in R and stays finite and non-negative', { skip, timeout: 300000 }, () => {
  const dir = tmpDir();
  LIBRARY.forEach((m) => fs.writeFileSync(path.join(dir, m.id + '.R'), generate(m).code));
  const r = spawnSync('Rscript', [path.join(here, 'r', 'run_all.R'), dir], { encoding: 'utf8' });
  assert.equal(r.status, 0, r.stdout + r.stderr);
});
