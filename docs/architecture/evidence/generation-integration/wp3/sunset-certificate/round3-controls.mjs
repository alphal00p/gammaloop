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
const card = (capture, change) => {
  for (const file of walk(capture).filter((file) => file.endsWith('-parameters.json'))) edit(file, change);
};
const dots = (capture, change) => {
  for (const file of walk(capture).filter((file) => file.endsWith('.dot'))) {
    fs.writeFileSync(file, change(fs.readFileSync(file, 'utf8')));
  }
};
const declaration = (name, value) => '    ' + name + ' = ' + JSON.stringify(JSON.stringify(value)) + ';\n';
const afterNum = (source, value) => source.replace('    num = "1";', '    num = "1";\n' + value);
const cases = [
  ['card-extra-entry', false, (c) => card(c, (x) => {x.unexpected = [1, 0];})],
  ['card-internal-ZERO', false, (c) => card(c, (x) => {x.ZERO = [0, 0];})],
  ['card-other-real', false, (c) => card(c, (x) => {x.Gf[0] = 0.9;})],
  ['card-other-imaginary', false, (c) => card(c, (x) => {x.Gf[1] = 0.1;})],
  ['card-single-sided', false, (c, f, route) => {if (route === 1) card(c, (x) => {x.aS[0] = 0.9;});}],
  ['positive-card-key-order', true, (c) => {
    for (const file of walk(c).filter((file) => file.endsWith('-parameters.json'))) {
      const x = JSON.parse(fs.readFileSync(file));
      fs.writeFileSync(file, JSON.stringify(Object.fromEntries(Object.entries(x).reverse())));
    }
  }],
  ['global-projector-value', false, (c) => dots(c, (s) => s.replace('projector = "1"', 'projector = "2"'))],
  ['global-projector-interface', false, (c) => dots(c, (s) => afterNum(s, declaration('projector_interface', ['spenso::mink(5)'])))],
  ['global-malformed-interface', false, (c) => dots(c, (s) => afterNum(s, declaration('numerator_prefactor_interface', null)))],
  ['global-interface-single-sided', false, (c, f, route) => {
    if (route === 1) dots(c, (s) => afterNum(s, declaration('projector_interface', ['spenso::mink(5)'])));
  }],
  ['global-tsv-single-sided', false, (c, f, route) => {
    if (route !== 1) return;
    for (const file of walk(c).filter((file) => file.endsWith('-graphs.tsv'))) {
      const lines = fs.readFileSync(file, 'utf8').trimEnd().split('\n');
      fs.writeFileSync(file, lines.map((line, i) => {if (!i) return line; const x = line.split('\t'); x[6] = '2'; return x.join('\t');}).join('\n') + '\n');
    }
  }],
  ['positive-global-scalar-and-empty-ports', true, (c) => {
    dots(c, (s) => afterNum(s, declaration('numerator_prefactor_interface', []) +
      declaration('projector_interface', []))
      .replace('num = "1"', 'num = "1+0"').replace('projector = "1"', 'projector = "1*1"'));
    for (const file of walk(c).filter((file) => file.endsWith('-graphs.tsv'))) {
      const lines = fs.readFileSync(file, 'utf8').trimEnd().split('\n');
      fs.writeFileSync(file, lines.map((line, i) => {if (!i) return line; const x = line.split('\t'); x[6] = '(1+0)'; return x.join('\t');}).join('\n') + '\n');
    }
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
      mutate(c, f, index);
      captures.push(c);
      frames.push(f);
    }
    const checks = [
      ['inputs.mjs', captures],
      ['model.mjs', [captures[0]]],
      ['model.mjs', [captures[1]]],
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
