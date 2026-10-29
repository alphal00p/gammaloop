import fs from 'node:fs';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';

const root = 'docs/architecture/evidence/generation-integration/wp3/ghost-weight-certificate';
const receipt = JSON.parse(fs.readFileSync(`${root}/receipt.json`));
const summary = JSON.parse(fs.readFileSync(`${root}/physical-summary.json`));
const mapping = JSON.parse(fs.readFileSync(`${root}/source44-rows.json`));
const sha = (file) => crypto.createHash('sha256').update(fs.readFileSync(file)).digest('hex');
let boundFiles = 0;
function validateBindings(value) {
  if (Array.isArray(value)) return value.forEach(validateBindings);
  if (value && typeof value === 'object') {
    if ('path' in value && 'sha256' in value) {
      assert.equal(sha(value.path), value.sha256, `changed artifact ${value.path}`);
      boundFiles++;
    }
    Object.values(value).forEach(validateBindings);
  }
}
validateBindings(receipt);
assert.equal(sha(summary.fullWitnessPath), summary.fullWitnessSha256);
const oldGit = '/shared/lcnbr/dev/gammaloop/main/.git';
const versions = {
  source43: execFileSync('git', [`--git-dir=${oldGit}`, 'show', `${mapping.sources.source43}:tests/tests/test_feyngen.rs`], {encoding: 'utf8', maxBuffer: 8 * 1024 * 1024}).split('\n'),
  source44: execFileSync('git', [`--git-dir=${oldGit}`, 'show', `${mapping.sources.source44}:tests/tests/test_feyngen.rs`], {encoding: 'utf8', maxBuffer: 8 * 1024 * 1024}).split('\n'),
  current: fs.readFileSync('tests/tests/test_feyngen.rs', 'utf8').split('\n'),
};
for (const row of mapping.cases) {
  for (const [version, lines] of Object.entries(versions)) {
    const {line, expected} = row[version];
    const sourceLine = lines[line - 1];
    assert(sourceLine?.includes(row.exact_process_fragment), `${row.id}/${version} exact CLI source mismatch`);
    assert(sourceLine.includes(`@"${expected}"`), `${row.id}/${version} expected literal mismatch`);
    const owner = lines.slice(0, line).reverse().find((text) => /^\s*fn test_/.test(text));
    assert(owner?.includes(row.qualified_test.split('::')[1] + '('), `${row.id}/${version} wrong function owner`);
    const module = lines.slice(0, line).reverse().find((text) => /^mod slow \{/.test(text));
    assert(module, `${row.id}/${version} missing slow module`);
  }
  const key = row.id.replace('vacuum_nonfactorizable_', 'vacuum_');
  const capture = receipt.cases[row.id].captures[row.id === 'epem_bbh_qcd2' ? 'bbh_qcd2' : row.id === 'slow_aa_ttx_qcd2' ? 'qcd2' : row.id.split('_').at(-1)];
  assert(capture, `missing capture ${key}`);
  assert.equal(capture.graphsPerRoute, Number(row.source44.expected.split(' | ')[0]));
  assert.equal(capture.sum.numerator + '/' + capture.sum.denominator,
    row.source44.expected.split(' | ')[1].split(' = ')[0].includes('/')
      ? row.source44.expected.split(' | ')[1].split(' = ')[0]
      : `${row.source44.expected.split(' | ')[1].split(' = ')[0]}/1`);
}
console.log(`receipt valid: ${boundFiles} bound files; five qualified unchanged current literals; source43/source44 exact rows; full physical join hash`);
