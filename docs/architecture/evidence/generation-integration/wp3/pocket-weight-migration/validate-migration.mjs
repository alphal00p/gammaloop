import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import { execFileSync } from 'node:child_process';

const root = process.cwd();
const evidence = 'docs/architecture/evidence/generation-integration/wp3/pocket-weight-migration/migration.json';
const receipt = JSON.parse(fs.readFileSync(evidence));
const hash = (value) => crypto.createHash('sha256').update(value).digest('hex');
function fileHash(file, expected) {
  if (hash(fs.readFileSync(file)) !== expected) throw new Error(`SHA-256 mismatch: ${file}`);
}
function sourceAt(revision) {
  return execFileSync('jj', ['file', 'show', '-r', revision, 'tests/tests/test_feyngen.rs'], {
    encoding: 'utf8',
    maxBuffer: 16 * 1024 * 1024,
  });
}
function commands(log) {
  return log.split('\n').map((line) => line.match(/Running command (generate xs .*)/)?.[1]).filter(Boolean);
}

const approvedSource = sourceAt(receipt.source.approved_certificate_bookmark);
const helperSource = sourceAt(receipt.source.helper_bookmark);
const currentSource = fs.readFileSync('tests/tests/test_feyngen.rs', 'utf8');
if (hash(approvedSource) !== receipt.source.approved_test_source_sha256 ||
    hash(helperSource) !== receipt.source.helper_test_source_sha256 ||
    hash(currentSource) !== receipt.source.migration_test_source_sha256) {
  throw new Error('Source revision binding changed');
}
const named = 'fn test_generate_aa_ttx_cross_section_slow_filter()';
const begin = helperSource.indexOf(named);
const end = helperSource.indexOf('\n    }', begin);
if (begin < 0 || end < 0) throw new Error('Named failing test is missing');
const old = '@"874 | -266 = -266"';
const next = '@"874 | -186 = -186"';
const body = helperSource.slice(begin, end);
if (body.split(old).length !== 2 || currentSource !== helperSource.slice(0, begin) + body.replace(old, next) + helperSource.slice(end)) {
  throw new Error('Migration differs from the one named snapshot replacement');
}

const oldLog = fs.readFileSync(receipt.controls.approved_old_fail.log, 'utf8');
const newLog = fs.readFileSync(receipt.controls.migration_new_pass.log, 'utf8');
fileHash(receipt.controls.approved_old_fail.binary, receipt.controls.approved_old_fail.binary_sha256);
fileHash(receipt.controls.approved_old_fail.log, receipt.controls.approved_old_fail.log_sha256);
fileHash(receipt.controls.helper_cleanup.binary, receipt.controls.helper_cleanup.binary_sha256);
fileHash(receipt.controls.helper_cleanup.log, receipt.controls.helper_cleanup.log_sha256);
fileHash(receipt.controls.migration_new_pass.binary, receipt.controls.migration_new_pass.binary_sha256);
fileHash(receipt.controls.migration_new_pass.log, receipt.controls.migration_new_pass.log_sha256);
if (commands(oldLog).length !== 3 || JSON.stringify(commands(oldLog)) !== JSON.stringify(commands(newLog)) ||
    JSON.stringify(commands(oldLog)) !== JSON.stringify(receipt.inputs.generate_commands)) {
  throw new Error('Old and new generation commands differ');
}
if (!oldLog.includes('874 | -266 = -266') || !oldLog.includes('874 | -186 = -186') ||
    !oldLog.includes('test failing::test_generate_aa_ttx_cross_section_slow_filter ... FAILED') ||
    !newLog.includes('test failing::test_generate_aa_ttx_cross_section_slow_filter ... ok') ||
    !newLog.includes('test result: ok. 1 passed; 0 failed;')) {
  throw new Error('Old-fail/new-pass result is missing');
}
const counts = [...newLog.matchAll(/Generated cross-section (\d+) graphs/g)].map((match) => Number(match[1]));
if (JSON.stringify(counts) !== JSON.stringify([4, 40, 874])) throw new Error('QCD graph count matrix changed');

for (const approval of Object.values(receipt.approvals)) {
  fileHash(approval.review, approval.sha256);
  if (JSON.parse(fs.readFileSync(approval.review)).verdict !== 'APPROVE_CERTIFICATE') {
    throw new Error('Certificate approval is missing');
  }
}
const frozen = receipt.approved_evidence_frozen;
fileHash(frozen.scope_b_raw_bindings_path, frozen.scope_b_raw_bindings_sha256);
const rawBindings = JSON.parse(fs.readFileSync(frozen.scope_b_raw_bindings_path));
if (Object.keys(rawBindings).length !== frozen.scope_b_raw_binding_count) {
  throw new Error('Scope B raw binding inventory changed');
}
for (const [file, binding] of Object.entries(rawBindings)) {
  const bytes = fs.readFileSync(file);
  if (bytes.length !== binding.bytes || hash(bytes) !== binding.sha256) {
    throw new Error(`Scope B raw input changed: ${file}`);
  }
}
fileHash('docs/architecture/evidence/generation-integration/wp3/pocket-weight-diagnosis/diagnosis.json', frozen.diagnosis_json_sha256);
fileHash('docs/architecture/evidence/generation-integration/wp3/pocket-weight-diagnosis/live-import-provenance.json', frozen.live_import_json_sha256);
fileHash('docs/architecture/evidence/generation-integration/wp3/pocket-weight-diagnosis/graph-join-index.json', frozen.graph_join_index_sha256);
fileHash('tests/artifacts/pocket_weight_diagnosis_round2/filtered-graph-witness.tsv', frozen.filtered_tsv_sha256);
fileHash('tests/artifacts/pocket_weight_diagnosis_round2/unfiltered-graph-witness.tsv', frozen.unfiltered_tsv_sha256);
fileHash(path.join('/shared/lcnbr/dev/gammaloop/generation-pocket-weight-historical-source43',
  'target/pocket-weight-diagnosis/binaries/source43-pristine-test_feyngen'), frozen.historical_pristine_binary_sha256);
for (const [file, expected] of [
  ['assets/models/json/sm/sm.json', receipt.inputs.source_model_sha256],
  ['assets/models/json/sm/restrict_default.json', receipt.inputs.default_restriction_source_sha256],
  ['Cargo.lock', receipt.inputs.cargo_lock_sha256],
  ['flake.lock', receipt.inputs.flake_lock_sha256],
  ['nix-ci.nix', receipt.inputs.nix_ci_sha256],
]) fileHash(file, expected);

const original = 'tests/artifacts/pocket_weight_diagnosis_round2';
const helper = 'target/pocket-weight-diagnosis/migration/helper-state/pocket_weight_diagnosis_round2';
for (const label of ['filtered', 'unfiltered']) {
  for (const suffix of ['graph-witness.tsv', 'process-definition.json']) {
    const file = `${label}-${suffix}`;
    if (!fs.readFileSync(path.join(original, file)).equals(fs.readFileSync(path.join(helper, file)))) {
      throw new Error(`Helper cleanup changed ${file}`);
    }
  }
  const fragment = `${label}-dot/processes/cross_sections/aa_ttx/default`;
  const before = path.join(original, fragment);
  const after = path.join(helper, fragment);
  const names = fs.readdirSync(before).sort();
  if (JSON.stringify(names) !== JSON.stringify(fs.readdirSync(after).sort()) ||
      names.length !== (label === 'filtered' ? 874 : 1106)) {
    throw new Error(`Helper cleanup changed ${label} DOT inventory`);
  }
  for (const name of names) {
    if (!fs.readFileSync(path.join(before, name)).equals(fs.readFileSync(path.join(after, name)))) {
      throw new Error(`Helper cleanup changed ${label} DOT ${name}`);
    }
  }
}
for (const stage of ['post-import', 'filtered-after-generation', 'unfiltered-after-generation']) {
  for (const kind of ['model', 'parameters']) {
    const file = `${stage}-live-${kind}.json`;
    if (!fs.readFileSync(path.join(original, file)).equals(fs.readFileSync(path.join(helper, file)))) {
      throw new Error(`Helper cleanup changed ${file}`);
    }
  }
}
console.log('single-row source diff, old-fail/new-pass, frozen certificates and helper capture verified');
