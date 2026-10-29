import fs from 'node:fs';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';

const root = 'docs/architecture/evidence/generation-integration/wp3/ghost-weight-certificate';
const oldReceipt = JSON.parse(fs.readFileSync(`${root}/receipt.json`));
const oldSummary = JSON.parse(fs.readFileSync(`${root}/physical-summary.json`));
const newSummary = JSON.parse(fs.readFileSync(`${root}/round2/physical-summary.json`));
const oldFull = JSON.parse(fs.readFileSync('target/ghost-weight-certificate/full-physical-joins.json'));
const newFull = JSON.parse(fs.readFileSync('target/ghost-weight-certificate/round2/full-physical-joins.json'));
const sha = (bytes) => crypto.createHash('sha256').update(bytes).digest('hex');
const hashFile = (location) => sha(fs.readFileSync(location));
const frozen = {
  receipt: '0e3679012f92c565f7d9895fd20a56a61c5243102bdd0003b433e8f3e96d19b7',
  summary: '48709016e63adbb9c0c755ac99b9287b611d438375adbc557f0c3d0e873edd6f',
  full: '1562ae44ef80613286f24b6784ea2abbd1d430844cb962ba9aa09662d8078bd6',
};
assert.equal(hashFile(`${root}/receipt.json`), frozen.receipt);
assert.equal(hashFile(`${root}/physical-summary.json`), frozen.summary);
assert.equal(hashFile('target/ghost-weight-certificate/full-physical-joins.json'), frozen.full);
let boundFiles = 0;
function validateBindings(value) {
  if (Array.isArray(value)) return value.forEach(validateBindings);
  if (!value || typeof value !== 'object') return;
  if ('path' in value && 'sha256' in value) {
    const bytes = value.path.endsWith('/ghost-weight-certificate/analyze-graphs.mjs')
      ? execFileSync('jj', ['file', 'show', '-r', 'generation-ghost-weight-certificate-round1',
        `${root}/analyze-graphs.mjs`], {maxBuffer: 8 * 1024 * 1024})
      : fs.readFileSync(value.path);
    assert.equal(sha(bytes), value.sha256, `changed frozen binding: ${value.path}`);
    boundFiles++;
  }
  Object.values(value).forEach(validateBindings);
}
validateBindings(oldReceipt);
assert.equal(boundFiles, 199, 'round1 receipt binding count');
assert.equal(hashFile(newSummary.fullWitnessPath), newSummary.fullWitnessSha256);
assert.deepEqual(Object.keys(oldFull.cases).sort(), Object.keys(newFull.cases).sort());
let graphRecords = 0;
let joins = 0;
let cuts = 0;
let thresholds = 0;
let properSubgraphThresholds = 0;
for (const [key, oldCase] of Object.entries(oldFull.cases)) {
  const next = newFull.cases[key];
  for (const field of ['historicalDefinition', 'currentDefinition',
    'historicalWitnessSha256', 'currentWitnessSha256', 'joins']) {
    assert.deepEqual(next[field], oldCase[field], `${key}/${field} changed`);
  }
  for (const route of ['historical', 'current']) {
    assert.equal(next[route].length, oldCase[route].length, `${key}/${route} graph count`);
    for (let index = 0; index < oldCase[route].length; index++) {
      const previous = oldCase[route][index];
      const corrected = next[route][index];
      for (const field of Object.keys(previous)) {
        assert.deepEqual(corrected[field], previous[field], `${key}/${route}/${index}/${field} changed`);
      }
      graphRecords++;
      cuts += corrected.cuts;
      thresholds += corrected.thresholds;
      properSubgraphThresholds += corrected.properSubgraphThresholds;
    }
  }
  joins += oldCase.joins.length;
  for (const [field, previous] of Object.entries(oldSummary.cases[key])) {
    assert.deepEqual(newSummary.cases[key][field], previous, `${key} summary ${field} changed`);
  }
}
assert.deepEqual({graphRecords, joins, cuts, thresholds, properSubgraphThresholds},
  {graphRecords: 8318, joins: 4159, cuts: 26368, thresholds: 73964, properSubgraphThresholds: 47716});
const controls = JSON.parse(fs.readFileSync('target/ghost-weight-certificate/round2/controls.json'));
assert.equal(controls.results.length, 33);
assert(controls.results.every((entry) => entry.expected === entry.outcome));
const report = {schema: 1, oldReceiptSha256: frozen.receipt, oldSummarySha256: frozen.summary,
  oldFullSha256: frozen.full, oldAnalyzerFromRef: 'generation-ghost-weight-certificate-round1',
  originalBoundFiles: boundFiles, unchangedOriginalRawAndDerivedFields: true,
  graphRecords, joins, cuts, thresholds, properSubgraphThresholds,
  correctedFullSha256: hashFile(newSummary.fullWitnessPath),
  correctedSummarySha256: hashFile(`${root}/round2/physical-summary.json`),
  controlsSha256: hashFile('target/ghost-weight-certificate/round2/controls.json')};
fs.writeFileSync('target/ghost-weight-certificate/round2/raw-guard.json', `${JSON.stringify(report, null, 2)}\n`);
let correctedBoundFiles = 0;
if (!process.argv.includes('--raw-only')) {
  const correction = JSON.parse(fs.readFileSync(`${root}/round2/correction.json`));
  const validateCorrectionBindings = (value) => {
    if (Array.isArray(value)) return value.forEach(validateCorrectionBindings);
    if (!value || typeof value !== 'object') return;
    if ('path' in value && 'sha256' in value) {
      assert.equal(hashFile(value.path), value.sha256, `changed corrected binding: ${value.path}`);
      correctedBoundFiles++;
    }
    Object.values(value).forEach(validateCorrectionBindings);
  };
  validateCorrectionBindings(correction);
}
console.log(`round2 guard: ${boundFiles} frozen and ${correctedBoundFiles} corrected bindings; ${graphRecords} unchanged records; ${joins} joins; ${cuts} cuts; ${thresholds} thresholds (${properSubgraphThresholds} proper subgraphs)`);
