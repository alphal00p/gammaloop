import fs from 'node:fs';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';

const base = 'docs/architecture/evidence/generation-integration/wp3/ghost-weight-certificate/round3';
const guard = execFileSync('node', [`${base}/validate-guard.mjs`], {encoding: 'utf8',
  maxBuffer: 8 * 1024 * 1024}).trim();
assert(guard.includes('199+40 frozen bindings'));
const receipt = JSON.parse(fs.readFileSync(`${base}/correction.json`));
assert.equal(receipt.status, 'corrected five-row scalar checker; NEW independent Astra correction review required');
assert.equal(receipt.checks.fullJoinByteEqualToRound2, true);
assert.deepEqual(receipt.checks.frozenBindings, [199, 40]);
assert.deepEqual(receipt.checks.differential, [44, 24]);
const sha = (location) => crypto.createHash('sha256').update(fs.readFileSync(location)).digest('hex');
let bindings = 0;
function validate(value) {
  if (Array.isArray(value)) return value.forEach(validate);
  if (!value || typeof value !== 'object') return;
  if ('path' in value && 'sha256' in value) {
    assert.equal(sha(value.path), value.sha256, `changed round3 receipt binding ${value.path}`);
    bindings++;
  }
  Object.values(value).forEach(validate);
}
validate(receipt);
const summary = JSON.parse(fs.readFileSync(`${base}/physical-summary.json`));
assert.equal(sha(summary.fullWitnessPath), summary.fullWitnessSha256);
assert.equal(receipt.corrected.fullJoin.sha256, summary.fullWitnessSha256);
const original = JSON.parse(fs.readFileSync('docs/architecture/evidence/generation-integration/wp3/ghost-weight-certificate/source44-rows.json'));
assert.equal(original.cases.length, 5);
const testLines = fs.readFileSync('tests/tests/test_feyngen.rs', 'utf8').split('\n');
for (const row of original.cases) {
  const source = testLines[row.current.line - 1];
  assert(source.includes(row.exact_process_fragment), `${row.id} current command changed`);
  assert(source.includes(`@"${row.current.expected}"`), `${row.id} current literal changed`);
  const owner = testLines.slice(0, row.current.line).reverse()
    .find((line) => /^\s*fn test_/.test(line));
  assert(owner?.includes(`${row.qualified_test.split('::')[1]}(`), `${row.id} current function changed`);
}
console.log(`${guard}; ${bindings} round3 bindings; five unchanged source44 row identities`);
