import fs from 'node:fs';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';

const root = 'docs/architecture/evidence/generation-integration/wp3';
const task = 'target/ghost-weight-certificate/migration';
const approval = 'generation-ghost-weight-certificate-approved';
const sha = (bytes) => crypto.createHash('sha256').update(bytes).digest('hex');
const hash = (file) => sha(fs.readFileSync(file));
const frozen = [
  {label: 'round1', receipt: `${root}/ghost-weight-certificate/receipt.json`,
    sha256: '0e3679012f92c565f7d9895fd20a56a61c5243102bdd0003b433e8f3e96d19b7',
    ref: 'generation-ghost-weight-certificate-round1', expected: 199},
  {label: 'round2', receipt: `${root}/ghost-weight-certificate/round2/correction.json`,
    sha256: 'ff3befa605c9a4eee915c3a0007eaee1ddd2a456b00683d5a7815f7d75e6c2ba',
    ref: 'generation-ghost-weight-certificate-round2', expected: 40},
  {label: 'approved', receipt: `${root}/ghost-weight-certificate/round3/correction.json`,
    sha256: 'ceaa4f77ae6d479ee32e6e96da586638f5e4c44938afd717ec4db3d760dc5749',
    ref: approval, expected: 44},
];
const refFile = (ref, file) => execFileSync('jj', ['file', 'show', '-r', ref, file],
  {maxBuffer: 8 * 1024 * 1024});
const counts = {};
for (const item of frozen) {
  assert.equal(hash(item.receipt), item.sha256, `${item.label} receipt moved`);
  const data = JSON.parse(fs.readFileSync(item.receipt));
  let bound = 0;
  const visit = (value) => {
    if (Array.isArray(value)) return value.forEach(visit);
    if (!value || typeof value !== 'object') return;
    if ('path' in value && 'sha256' in value) {
      const relative = value.path.replace(`${process.cwd()}/`, '');
      const sourceInRef = relative === 'tests/tests/test_feyngen.rs' ? approval
        : relative === 'docs/developers.toml' ? item.ref
          : relative.endsWith('/ghost-weight-certificate/analyze-graphs.mjs')
            && item.label !== 'approved' ? item.ref
              : relative === 'docs/architecture/generation-ghost-weight-certificate.typ'
                && item.label === 'round2' ? item.ref : null;
      const bytes = sourceInRef ? refFile(sourceInRef, relative) : fs.readFileSync(value.path);
      assert.equal(sha(bytes), value.sha256, `${item.label} binding changed: ${value.path}`);
      bound++;
    }
    Object.values(value).forEach(visit);
  };
  visit(data);
  assert.equal(bound, item.expected, `${item.label} binding count`);
  counts[item.label] = bound;
}
const oldFull = 'target/ghost-weight-certificate/full-physical-joins.json';
const round2Full = 'target/ghost-weight-certificate/round2/full-physical-joins.json';
const round3Full = 'target/ghost-weight-certificate/round3/full-physical-joins.json';
assert.equal(hash(oldFull), '1562ae44ef80613286f24b6784ea2abbd1d430844cb962ba9aa09662d8078bd6');
assert.equal(hash(round2Full), '7df7aac4a48754ac23c6f8dfa0c177cb40d82972f4d12cff2675ee140a3fa747');
assert.equal(hash(round3Full), hash(round2Full));
const reviewRoot = 'target/ghost-weight-certificate/review-astra-round3-combined-5LJwn8';
const review = JSON.parse(fs.readFileSync(`${reviewRoot}/review.json`));
assert.equal(review.verdict, 'APPROVE_CERTIFICATE');
let reviewedArtifacts = 0;
for (const [name, expected] of Object.entries(review.artifacts)) {
  assert.equal(hash(`${reviewRoot}/${name}`), expected, `approved review artifact changed: ${name}`);
  reviewedArtifacts++;
}
execFileSync('node', [`${root}/ghost-weight-migration/check-five-row-diff.mjs`]);
execFileSync('node', [`${root}/ghost-weight-migration/verify-test-logs.mjs`]);
execFileSync('node', [`${root}/ghost-weight-migration/verify-final-capture.mjs`]);
const diff = JSON.parse(fs.readFileSync(`${task}/five-row-diff.json`));
const capture = JSON.parse(fs.readFileSync(`${task}/final-capture-comparison.json`));
assert.equal(diff.changedLineCount, 5);
assert.deepEqual(capture.counts, {dots: 4239, graphWitnesses: 15, processDefinitions: 15, live: 63});
assert.equal(capture.oldByteTreeSha256IgnoringSettingsFolder,
  capture.finalByteTreeSha256IgnoringSettingsFolder);
const oldBinary = `${task}/old-literals-test_feyngen`;
const newBinary = `${task}/new-literals-test_feyngen`;
assert.equal(hash(oldBinary), 'fc4cccfde6411b2036b20a3c0e58189618f667e4b696aba50da5faed5c66d948');
assert.equal(hash(newBinary), '811bc93ca99e97497cf7b0397fed9418761886e8b449b3b2602b480e62631da3');
const result = {schema: 1, approvedCertificateCommit: review.candidate.commit,
  reviewPath: `${reviewRoot}/review.json`, reviewSha256: hash(`${reviewRoot}/review.json`),
  reviewedArtifacts, frozenBindings: counts, frozenFullJoinSha256: [hash(oldFull), hash(round2Full), hash(round3Full)],
  oldBinary: {path: oldBinary, sha256: hash(oldBinary)},
  newBinary: {path: newBinary, sha256: hash(newBinary)},
  oldSourceSha256: diff.originalSourceSha256,
  newSourceSha256: diff.migratedSourceSha256,
  sourceModelSha256: hash('assets/models/json/sm/sm.json'),
  cargoLockSha256: hash('Cargo.lock'), flakeLockSha256: hash('flake.lock'),
  finalCapture: capture,
  exactFiveLineDiff: true,
  proofSourceAndRawInputsStable: true};
fs.writeFileSync(`${task}/proof-guard.json`, `${JSON.stringify(result, null, 2)}\n`);
console.log(`migration proof guard: ${counts.round1}+${counts.round2}+${counts.approved} frozen bindings; ${reviewedArtifacts} approved review artifacts; five lines only; ${capture.totalSelectedFiles} final captures equal`);
