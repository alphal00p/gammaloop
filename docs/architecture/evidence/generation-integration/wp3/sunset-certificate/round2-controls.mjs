import './lossless-json.mjs';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import {spawnSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';

const scripts = path.dirname(fileURLToPath(import.meta.url));
const [historical, current, historicalFrames, currentFrames, output] = process.argv.slice(2);
if (!output) throw new Error('usage: round2-controls <historical> <current> <historical-frames> <current-frames> <output>');
const walk = (root) => fs.readdirSync(root, {recursive: true, withFileTypes: true})
  .filter((entry) => entry.isFile()).map((entry) => path.join(entry.parentPath, entry.name));
const edit = (file, change) => {
  const value = JSON.parse(fs.readFileSync(file, 'utf8'));
  change(value);
  fs.writeFileSync(file, JSON.stringify(value, null, 2) + '\n');
};
function allDiagrams(capture, frames, change) {
  for (const file of walk(capture).filter((file) => /native-(source|final)\.json$/.test(file))) {
    edit(file, (diagrams) => diagrams.forEach(change));
  }
  for (const file of walk(frames).filter((file) => file.endsWith('.json'))) edit(file, change);
}
const cases = [
  ['positive-member-order', true, (capture) => {
    for (const file of walk(capture).filter((file) => file.endsWith('native-groups.json'))) {
      edit(file, (value) => value.groups.forEach((group) => group.members.reverse()));
    }
  }],
  ['positive-dot-vertex-relabel', true, (capture) => {
    for (const file of walk(capture).filter((file) => file.endsWith('.dot'))) {
      fs.writeFileSync(file, fs.readFileSync(file, 'utf8').replace(/\bv(\d+)\b/g,
        (_, id) => `v${Number(id) + 17}`));
    }
  }],
  ['mirrored-request-seed', false, (capture) => {
    for (const file of walk(capture).filter((file) => file.endsWith('grouped-definition.json'))) {
      edit(file, (value) => {value.grouping.up_to_scalar.numerical_sample_seed = 999;});
    }
  }],
  ['mirrored-native-cli-ghost-coupling', false, (capture, frames) => {
    for (const file of walk(capture).filter((file) => file.endsWith('.dot'))) {
      fs.writeFileSync(file, fs.readFileSync(file, 'utf8').replaceAll('GC_10', 'GC_12'));
    }
    allDiagrams(capture, frames, (diagram) => {
      if (diagram.name !== 'GL1') return;
      diagram.numerator = diagram.numerator.replaceAll('GC_10', 'GC_12');
      for (const vertex of diagram.vertices) vertex.numerator = vertex.numerator.replaceAll('GC_10', 'GC_12');
    });
  }],
  ['mirrored-native-cli-gluon-coefficient', false, (capture, frames) => {
    for (const file of walk(capture).filter((file) => /\/GL0\.dot$/.test(file))) {
      fs.writeFileSync(file, fs.readFileSync(file, 'utf8').replace(/(\bv\d+ \[int_id="[^"]+" num=")([^"]+)(")/g,
        '$1' + '2*($2)' + '$3'));
    }
    allDiagrams(capture, frames, (diagram) => {
      if (diagram.name !== 'GL0') return;
      diagram.numerator = `4*(${diagram.numerator})`;
      for (const vertex of diagram.vertices) vertex.numerator = `2*(${vertex.numerator})`;
    });
  }],
  ['mirrored-native-only-foreign-head', false, (capture, frames) => {
    allDiagrams(capture, frames, (diagram) => {
      diagram.overall_factor = diagram.overall_factor.replace(/feynkit_generator_factor::(?:\{\}::)?AutG/g, 'Bogus::AutG');
    });
  }],
];
const results = [];
for (const [name, positive, mutate] of cases) {
  const scratch = fs.mkdtempSync(path.join(os.tmpdir(), 'sunset-chain-control-'));
  try {
    const captures = [];
    const frames = [];
    for (const [index, [capture, frame]] of [[historical, historicalFrames], [current, currentFrames]].entries()) {
      const c = path.join(scratch, String(index), 'capture');
      const f = path.join(scratch, String(index), 'frames');
      fs.cpSync(capture, c, {recursive: true, filter: (file) => !['state', 'drawings'].includes(path.basename(file))});
      fs.cpSync(frame, f, {recursive: true});
      mutate(c, f);
      captures.push(c);
      frames.push(f);
    }
    const checks = [
      ['inputs.mjs', captures],
      ['join.mjs', captures],
      ['frame.mjs', captures],
      ['internal-frame.mjs', [...frames, ...captures]],
    ].map(([script, args]) => {
      const run = spawnSync(process.execPath, [path.join(scripts, script), ...args,
        path.join(scratch, `${script}.json`)], {encoding: 'utf8', env: process.env});
      return {script, exit: run.status, stderr: run.stderr.slice(-1400).replaceAll(scratch, '<fixture>')};
    });
    const pass = checks.every((check) => check.exit === 0);
    if (pass !== positive) throw new Error(`${name}: unexpected semantic result ${JSON.stringify(checks)}`);
    if (name.includes('ghost-coupling') || name.includes('gluon-coefficient')) {
      if (!checks.some((check) => check.stderr.includes('complete live local rule differs') ||
          check.stderr.includes('declared tensor/layout differs'))) {
        throw new Error(`${name}: did not reach production local-rule boundary`);
      }
    }
    results.push({name, positive, mirroredBothProducers: true, shaGuards: false, checks});
  } finally {
    fs.rmSync(scratch, {recursive: true, force: true});
  }
}
fs.writeFileSync(output, JSON.stringify({results}, null, 2) + '\n');
console.log(`${results.filter((test) => test.positive).length} positive transforms; ` +
  `${results.filter((test) => !test.positive).length} new complete-pipeline negatives`);
