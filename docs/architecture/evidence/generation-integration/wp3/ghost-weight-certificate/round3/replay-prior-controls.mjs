import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';
import {pathToFileURL} from 'node:url';

const out = 'target/ghost-weight-certificate/round3';
const prior = 'docs/architecture/evidence/generation-integration/wp3/ghost-weight-certificate/round2/controls.mjs';
const review = 'target/ghost-weight-certificate/review-astra-round2-scope-b-EuYuaM';
const analyzer = pathToFileURL(path.resolve('docs/architecture/evidence/generation-integration/wp3/ghost-weight-certificate/analyze-graphs.mjs')).href;
const independent = pathToFileURL(path.resolve(review, 'independent.mjs')).href;
const sha = (location) => crypto.createHash('sha256').update(fs.readFileSync(location)).digest('hex');
const replaceOnce = (source, oldText, newText) => {
  assert(source.includes(oldText), `missing review-source adaptation marker: ${oldText}`);
  assert.equal(source.split(oldText).length, 2, `nonunique review-source adaptation marker: ${oldText}`);
  return source.replace(oldText, newText);
};
const run = (name, input, changes, output, result) => {
  let source = fs.readFileSync(input, 'utf8');
  for (const [before, after] of changes) source = replaceOnce(source, before, after);
  fs.writeFileSync(output, source);
  const log = execFileSync('node', [output], {encoding: 'utf8', maxBuffer: 8 * 1024 * 1024});
  fs.writeFileSync(`${out}/${name}.log`, log);
  const receipt = {name, original: {path: input, sha256: sha(input)},
    adapted: {path: output, sha256: sha(output)}, substitutions: changes,
    output: {path: result, sha256: sha(result)}, log};
  return receipt;
};
const originalRound2 = fs.readFileSync(prior, 'utf8');
const outputLine = originalRound2.match(/^const output = (.+);$/m)?.[1];
assert(outputLine);
const records = [];
records.push(run('replayed-round2-controls', prior,
  [["from '../analyze-graphs.mjs';", `from '${analyzer}';`],
    [`const output = ${outputLine};`, `const output = '${out}/replayed-round2-controls.json';`]],
  `${out}/replayed-round2-controls.mjs`, `${out}/replayed-round2-controls.json`));
records.push(run('replayed-review-extra-controls', `${review}/extra-controls.mjs`, [
  ["import {api} from './candidate-api.mjs';", `import * as api from '${analyzer}';`],
  ["target/ghost-weight-certificate/round2/full-physical-joins.json", `${out}/full-physical-joins.json`],
], `${out}/replayed-review-extra-controls.mjs`, `${out}/extra-controls.json`));
records.push(run('replayed-review-prior-controls', `${review}/controls.mjs`, [
  ["import {api,load} from './candidate-api.mjs';", `import * as api from '${analyzer}';`],
  ["import {parse,components,automorphisms,factor} from './independent.mjs';",
    `import {parse,components,automorphisms,factor} from '${independent}';`],
  ["const a=load(fileSystem);return {a,data:a.readCase(current,vacConfig,'l2')};",
    "const a=api;return {a,data:a.readCase(current,vacConfig,'l2',fileSystem)};"],
  ["const receiptCheck=execFileSync('node',[path.join(dir,'validate-frozen.mjs')],{encoding:'utf8',maxBuffer:8*1024*1024});",
    "const receiptCheck='frozen bindings independently validated by round3 guard';"],
], `${out}/replayed-review-prior-controls.mjs`, `${out}/controls.json`));
const priorResult = JSON.parse(fs.readFileSync(`${out}/replayed-round2-controls.json`));
assert.equal(priorResult.results.length, 33);
assert(priorResult.results.every((entry) => entry.outcome === entry.expected));
assert.equal(priorResult.originalControls.length, 7);
const extra = JSON.parse(fs.readFileSync(`${out}/extra-controls.json`));
assert.equal(extra.total, 111);
assert.equal(extra.generatedExpressionComparisons, 2000);
assert.deepEqual(extra.unexpected, []);
const reviewer = JSON.parse(fs.readFileSync(`${out}/controls.json`));
assert.equal(reviewer.results.length, 31);
assert(reviewer.results.every((entry) => entry.meetsExpectation));
fs.writeFileSync(`${out}/prior-controls-replay.json`, `${JSON.stringify({schema: 1, records,
  originalSeven: priorResult.originalControls.length, round2Controls: priorResult.results.length,
  reviewerPriorControls: reviewer.results.length, reviewerExtraControls: extra.total,
  generatedScalarExpressions: extra.generatedExpressionComparisons}, null, 2)}\n`);
console.log(`prior controls replay: 7 original, 33 round2, 31 reviewer prior, 111 reviewer extra, 2000 generated scalars; all expected outcomes`);
