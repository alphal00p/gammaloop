import './lossless-json.mjs';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import {spawnSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';

const [historical, current, historicalCapture, currentCapture, output] = process.argv.slice(2);
if (!historical || !current || !historicalCapture || !currentCapture || !output) {
  throw new Error('usage: node internal-mutations.mjs <historical-frames> <current-frames> <historical-capture> <current-capture> <output.json>');
}
const scripts = path.dirname(fileURLToPath(import.meta.url));
const bothRuns = (root, id, stage, change) => {
  for (const run of [0, 1]) {
    const file = path.join(root, `run-${run}`, `source-${id}-${stage}.json`);
    const value = JSON.parse(fs.readFileSync(file, 'utf8'));
    change(value);
    fs.writeFileSync(file, JSON.stringify(value, null, 2) + '\n');
  }
};
const bothFrames = (root, id, before, after) => {
  for (const run of [0, 1]) {
    const file = path.join(root, `run-${run}`, `source-${id}.txt`);
    const value = fs.readFileSync(file, 'utf8');
    if (!value.includes(before)) throw new Error(`missing frame mutation anchor: ${before}`);
    fs.writeFileSync(file, value.replace(before, after));
  }
};
const cases = [
  ['signed comparison edge flip', (root) => bothFrames(root, 4,
    'edge_map\tcomparison\t0\t0\t1', 'edge_map\tcomparison\t0\t0\t0'),
  'complete comparison edge map'],
  ['signed affine dependent coefficient', (root) => bothRuns(root, 4, 'comparison', (data) => {
    data.loop_momentum_basis.edge_signatures['2'].loops[1] = 1;
  }), 'signed affine edge 2'],
  ['antiparticle owner', (root) => bothRuns(root, 4, 'comparison', (data) => {
    data.edges[0][1].particle = 27;
  }), 'provisional comparison oriented edge differs'],
  ['ghost local momentum sign', (root) => bothRuns(root, 1, 'comparison', (data) => {
    data.vertices[0].numerator = data.vertices[0].numerator.replace(/^-1\*/, '');
  }), 'aggregate differs from complete local product'],
  ['massive local coefficient', (root) => bothRuns(root, 3, 'comparison', (data) => {
    data.edges[0][1].numerator = data.edges[0][1].numerator.replace('UFO::{}::MT',
      '2*UFO::{}::MT');
  }), 'aggregate differs from complete local product'],
  ['comparison spin dimension', (root) => bothRuns(root, 4, 'comparison', (data) => {
    data.edges[0][1].num_interface[0] = data.edges[0][1].num_interface[0]
      .replace('bis(4', 'bis(5');
  }), 'complete logical ports edge 0'],
  ['comparison aggregate factor', (root) => bothRuns(root, 4, 'comparison', (data) => {
    data.numerator += '*2';
  }), 'aggregate differs from complete local product'],
  ['final edge owner', (root) => bothFrames(root, 4,
    'final_edge_map\t0\t0', 'final_edge_map\t0\t1'), 'final edge owner map'],
];
const results = [];
for (const [name, change, expected] of cases) {
  const scratch = fs.mkdtempSync(path.join(os.tmpdir(), 'sunset-internal-mutation-'));
  try {
    const clone = path.join(scratch, 'frames');
    fs.cpSync(current, clone, {recursive: true});
    change(clone);
    const run = spawnSync(process.execPath, [path.join(scripts, 'internal-frame.mjs'),
      historical, clone, historicalCapture, currentCapture, path.join(scratch, 'proof.json')],
    {encoding: 'utf8'});
    const detail = `${run.stdout}\n${run.stderr}`;
    if (run.status === 0 || !detail.includes(expected)) {
      throw new Error(`${name}: missing semantic rejection: ${detail.slice(0, 800)}`);
    }
    results.push({name, expected, exit: run.status});
  } finally {
    fs.rmSync(scratch, {recursive: true, force: true});
  }
}
fs.writeFileSync(output, JSON.stringify({shaGuards: false, results}, null, 2) + '\n');
console.log(`${results.length} internal-frame mutations rejected without SHA guards`);
