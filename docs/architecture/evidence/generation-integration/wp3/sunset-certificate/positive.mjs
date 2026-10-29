import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import crypto from 'node:crypto';
import {spawnSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';

const [historical, current, frozen, output] = process.argv.slice(2);
if (!historical || !current || !frozen || !output) {
  throw new Error('usage: node positive.mjs <historical-capture> <current-capture> <frozen-proof> <output.json>');
}
const scripts = path.dirname(fileURLToPath(import.meta.url));
const historicalFrames = path.resolve(historical, '../../../target/sunset-certificate/frames/source43');
const currentFrames = path.resolve(current, '../../../target/sunset-certificate/frames/current');
const hash = (file) => crypto.createHash('sha256').update(fs.readFileSync(file)).digest('hex');
const tests = [
  ['round3-controls.mjs', [historical, current, historicalFrames, currentFrames], 'correction-controls.json'],
  ['scalar-differential.mjs', [], 'scalar-differential.json'],
  ['round2-controls.mjs', [historical, current, historicalFrames, currentFrames], 'chain-controls.json'],
  ['inputs.mjs', [historical, current], 'inputs.json'],
  ['model.mjs', [current], 'model.json'],
  ['join.mjs', [historical, current], 'join.json'],
  ['frame.mjs', [historical, current], 'frame.json'],
  ['internal-frame.mjs', [historicalFrames, currentFrames, historical, current], 'internal-frame.json'],
  ['internal-mutations.mjs', [historicalFrames, currentFrames, historical, current],
    'internal-mutations.json'],
  ['mutations.mjs', [historical, current], 'mutations.json'],
];
const scratch = fs.mkdtempSync(path.join(os.tmpdir(), 'sunset-positive-'));
const results = [];
try {
  for (const [script, args, name] of tests) {
    const artifact = path.join(scratch, name);
    const run = spawnSync(process.execPath, [path.join(scripts, script), ...args, artifact],
      {encoding: 'utf8'});
    if (run.status !== 0) throw new Error(`${script} failed: ${run.stderr.slice(0, 800)}`);
    const actual = hash(artifact);
    const expected = hash(path.join(frozen, name));
    if (actual !== expected) throw new Error(`${script}: fresh result differs from retained proof`);
    results.push({script, proof: name, sha256: actual, exit: run.status});
  }
} finally {
  fs.rmSync(scratch, {recursive: true, force: true});
}
fs.writeFileSync(output, JSON.stringify({results}, null, 2) + '\n');
console.log(`${results.length} fresh positive proof replays match retained outputs`);
