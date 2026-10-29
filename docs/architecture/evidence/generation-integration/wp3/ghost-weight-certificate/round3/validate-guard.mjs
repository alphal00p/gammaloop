import fs from 'node:fs';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';

const base = 'docs/architecture/evidence/generation-integration/wp3/ghost-weight-certificate';
const target = 'target/ghost-weight-certificate/round3';
const sha = (bytes) => crypto.createHash('sha256').update(bytes).digest('hex');
const hashFile = (location) => sha(fs.readFileSync(location));
const frozen = {
  round1Receipt: '0e3679012f92c565f7d9895fd20a56a61c5243102bdd0003b433e8f3e96d19b7',
  round1Summary: '48709016e63adbb9c0c755ac99b9287b611d438375adbc557f0c3d0e873edd6f',
  round1Full: '1562ae44ef80613286f24b6784ea2abbd1d430844cb962ba9aa09662d8078bd6',
  round2Correction: 'ff3befa605c9a4eee915c3a0007eaee1ddd2a456b00683d5a7815f7d75e6c2ba',
  round2Summary: 'f221b65d41240752189c8327a4db165c8865f9de585e52155e37312f97d80739',
  round2Full: '7df7aac4a48754ac23c6f8dfa0c177cb40d82972f4d12cff2675ee140a3fa747',
};
assert.equal(hashFile(`${base}/receipt.json`), frozen.round1Receipt);
assert.equal(hashFile(`${base}/physical-summary.json`), frozen.round1Summary);
assert.equal(hashFile('target/ghost-weight-certificate/full-physical-joins.json'), frozen.round1Full);
assert.equal(hashFile(`${base}/round2/correction.json`), frozen.round2Correction);
assert.equal(hashFile(`${base}/round2/physical-summary.json`), frozen.round2Summary);
assert.equal(hashFile('target/ghost-weight-certificate/round2/full-physical-joins.json'), frozen.round2Full);
const original = JSON.parse(fs.readFileSync(`${base}/receipt.json`));
const prior = JSON.parse(fs.readFileSync(`${base}/round2/correction.json`));
const readFromRef = (ref, location) => execFileSync('jj', ['file', 'show', '-r', ref, location],
  {maxBuffer: 8 * 1024 * 1024});
function validateBindings(value, ref) {
  let count = 0;
  const visit = (item) => {
    if (Array.isArray(item)) return item.forEach(visit);
    if (!item || typeof item !== 'object') return;
    if ('path' in item && 'sha256' in item) {
      const relative = item.path.replace(`${process.cwd()}/`, '');
      const frozenSource = relative === `${base}/analyze-graphs.mjs`
        || ref === 'generation-ghost-weight-certificate-round2'
          && ['docs/architecture/generation-ghost-weight-certificate.typ', 'docs/developers.toml'].includes(relative);
      const bytes = frozenSource ? readFromRef(ref, relative) : fs.readFileSync(item.path);
      assert.equal(sha(bytes), item.sha256, `changed ${ref} binding ${item.path}`);
      count++;
    }
    Object.values(item).forEach(visit);
  };
  visit(value);
  return count;
}
const originalBindings = validateBindings(original, 'generation-ghost-weight-certificate-round1');
const round2Bindings = validateBindings(prior, 'generation-ghost-weight-certificate-round2');
assert.equal(originalBindings, 199);
assert.equal(round2Bindings, 40);
const summary = JSON.parse(fs.readFileSync(`${base}/round3/physical-summary.json`));
const previousSummary = JSON.parse(fs.readFileSync(`${base}/round2/physical-summary.json`));
assert.deepEqual(summary.cases, previousSummary.cases, 'captured graph membership/factors/sums changed');
assert.deepEqual(summary.semanticMutationControls, previousSummary.semanticMutationControls);
assert.equal(summary.fullWitnessSha256, frozen.round2Full);
assert.equal(hashFile(summary.fullWitnessPath), frozen.round2Full);
const controls = JSON.parse(fs.readFileSync(`${target}/controls.json`));
assert.equal(controls.results.length, 30);
assert(controls.results.every((entry) => entry.expected === entry.outcome));
const priorControls = JSON.parse(fs.readFileSync(`${target}/prior-controls-replay.json`));
assert.deepEqual([priorControls.originalSeven, priorControls.round2Controls,
  priorControls.reviewerPriorControls, priorControls.reviewerExtraControls,
  priorControls.generatedScalarExpressions], [7, 33, 31, 111, 2000]);
const differential = JSON.parse(fs.readFileSync(`${target}/parser-differential.json`));
assert.deepEqual([differential.admittedCount, differential.explicitlyRejectedCount], [44, 24]);
const report = {schema: 1, frozen, originalBindings, round2Bindings,
  unchangedRound2PhysicalFullBytes: true, unchangedRound2SummaryCases: true,
  graphRecords: 8318, joins: 4159, cuts: 26368, thresholds: 73964, properSubgraphThresholds: 47716,
  originalControls: 7, round2Controls: 33, reviewerPriorControls: 31,
  reviewerExtraControls: 111, generatedScalarComparisons: 2000,
  round3Controls: controls.results.length, differentialAdmitted: differential.admittedCount,
  differentialRejected: differential.explicitlyRejectedCount,
  currentDiagnosticSourceSha256: original.producer.current.diagnosticSource.sha256,
  historicalDiagnosticSourceSha256: original.producer.source43.diagnosticSource.sha256,
  currentBinarySha256: original.producer.current.binary.sha256,
  historicalBinarySha256: original.producer.source43.binary.sha256,
  liveModelSha256: original.families.epem_bbh.routes.current.modelSha256,
  resolvedParametersSha256: original.families.epem_bbh.routes.current.resolvedParametersSha256};
fs.writeFileSync(`${target}/raw-guard.json`, `${JSON.stringify(report, null, 2)}\n`);
console.log(`round3 guard: ${originalBindings}+${round2Bindings} frozen bindings; byte-identical round2 physical join; ${controls.results.length} new controls; ${differential.admittedCount}/${differential.explicitlyRejectedCount} scalar differential`);
