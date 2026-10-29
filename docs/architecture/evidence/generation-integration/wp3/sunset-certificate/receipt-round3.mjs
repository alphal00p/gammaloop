import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';

const [destination] = process.argv.slice(2);
if (!destination) throw new Error('usage: receipt-round3.mjs <receipt-round3.json>');
const scripts = path.dirname(fileURLToPath(import.meta.url));
const workspace = path.resolve(scripts, '../../../../../..');
const root = path.join(workspace, 'target/sunset-certificate/round3');
const read = (file) => JSON.parse(fs.readFileSync(file, 'utf8'));
const sha256 = (bytes) => crypto.createHash('sha256').update(bytes).digest('hex');
const artifact = (file) => ({path: path.resolve(file), bytes: fs.statSync(file).size,
  sha256: sha256(fs.readFileSync(file))});
const jj = (...args) => execFileSync('jj', ['--ignore-working-copy', ...args],
  {cwd: workspace, encoding: 'utf8'}).trim();

const prior = read(path.join(scripts, 'receipt-round2.json'));
const priorBindings = [];
function bindings(value) {
  if (!value || typeof value !== 'object') return;
  if (value.path && value.sha256) priorBindings.push(value);
  else Object.values(value).forEach(bindings);
}
bindings(prior);
const retained = priorBindings.map((entry) => {
  const live = artifact(entry.path);
  if (live.sha256 === entry.sha256) return {path: entry.path, sha256: entry.sha256,
    retainedAt: entry.path};
  const relative = path.relative(workspace, entry.path);
  assert(!relative.startsWith('..'), `foreign round2 artifact changed: ${entry.path}`);
  const archives = [path.join(workspace, 'target/sunset-certificate/round2/round1-original-files', relative),
    path.join(root, 'round2-original-files', relative)];
  const archived = archives.filter(fs.existsSync).map(artifact)
    .find((candidate) => candidate.sha256 === entry.sha256);
  assert(archived, `round2 artifact lost: ${relative}`);
  assert.equal(archived.sha256, entry.sha256, `round2 artifact lost: ${relative}`);
  return {path: entry.path, sha256: entry.sha256, retainedAt: archived.path};
});
const protectedFiles = ['Cargo.lock', 'flake.lock', 'nix-ci.nix'].map((name) => {
  const current = artifact(path.join(workspace, name));
  const frozen = execFileSync('jj', ['file', 'show', '-r', '6732cfbb', name], {cwd: workspace});
  assert.equal(current.sha256, sha256(frozen), `${name}: protected bytes changed`);
  return current;
});
const refs = Object.fromEntries(['generation-sunset-certificate',
  'generation-sunset-certificate-round1', 'generation-sunset-certificate-round2',
  'generation-sunset-historical-source43'].map((name) =>
  [name, jj('log', '-r', name, '--no-graph', '-T', 'commit_id')]));
assert.equal(refs['generation-sunset-certificate'], '6732cfbb59099e23c315665440dbdae67a3944d4');
assert.equal(refs['generation-sunset-certificate-round1'], refs['generation-sunset-certificate']);
assert.equal(refs['generation-sunset-certificate-round2'], '6aa6d311432d64b68e9208d3db0326aafbb080a3');
assert.equal(refs['generation-sunset-historical-source43'],
  '74fcd8e57e9fe256d4bc34b5d8ad93c97fc28848');
const controls = read(path.join(root, 'correction-controls.json')).results;
assert.equal(controls.length, 12);
assert.equal(controls.filter((test) => test.positive).length, 2);
assert(controls.every((test) => test.checks.every((check) => check.exit !== null)));
const sourceClosure = read(path.join(root, 'source-closure.json'));
const probe = artifact(path.join(workspace,
  'crates/feynkit-generator/examples/sunset_live_rule_probe.rs'));
assert.equal(probe.sha256, artifact(path.join(root, 'sunset-live-rule-probe.rs')).sha256);
const binary = artifact(path.join(root, 'binaries/sunset-live-rule-probe'));
assert.equal(binary.sha256,
  artifact(path.join(workspace, 'target/debug/examples/sunset_live_rule_probe')).sha256);
const changed = jj('diff', '--summary')
  .split('\n').filter(Boolean);
const nonReceiptPaths = changed.map((line) => line.slice(2))
  .filter((name) => !name.endsWith('/receipt-round3.json'));
const nonReceiptStats = jj('diff', '--stat', '--', ...nonReceiptPaths);
const sourceScripts = fs.readdirSync(scripts).filter((name) => name.endsWith('.mjs'))
  .sort().map((name) => artifact(path.join(scripts, name)));
const artifacts = ['inputs.json', 'model.json', 'join.json', 'frame.json',
  'internal-frame.json', 'mutations.json', 'internal-mutations.json',
  'scalar-differential.json', 'chain-controls.json', 'correction-controls.json',
  'positive.json', 'source-closure.json', 'title-audit.json',
  'a-round2-old.log', 'a-round2.log', 'b-round2.log', 'check-proof-final.log',
  'build-final.log', 'clippy.log', 'rustfmt.log', 'docs-check.log',
  'a-round2-old.mjs', 'a-round2.mjs', 'b-round2.mjs',
  'check-proof-final.sh', 'docs-check.sh', 'source-closure.mjs',
  'note.html', 'note.pdf'].map((name) => artifact(path.join(root, name)));
const result = {
  schema: 3, status: 'corrected bounded candidate; STOP for fresh independent Astra review',
  owner: {name: 'Lucien Huber', email: 'im@lcnbr.ch'},
  workspace, changeId: 'wmuropzp', bookmark: 'generation-sunset-certificate-round3',
  parent: refs['generation-sunset-certificate-round2'], preservedReferences: refs,
  scope: 'Only two grouped two-loop vacuum sunset rows, paired raw controls and certificate repairs.',
  originalReceipt: artifact(path.join(scripts, 'receipt.json')),
  round2Receipt: artifact(path.join(scripts, 'receipt-round2.json')),
  round2BindingsRetained: retained,
  round2Accounting: {bindings: retained.length,
    unchanged: retained.filter((entry) => entry.retainedAt === entry.path).length,
    archived: retained.filter((entry) => entry.retainedAt !== entry.path).length},
  protectedFiles, sourceScripts, changed, nonReceiptStats,
  sourceClosure: artifact(path.join(root, 'source-closure.json')),
  documentation: ['docs/architecture/generation-sunset-certificate.typ',
    'docs/developers.toml'].map((name) => artifact(path.join(workspace, name))),
  probe: {source: probe, binary, actualClosure: sourceClosure,
    features: 'Cargo default/native; current Symbolica kernel 1a1829b1f3d0a88b53b44085b3b9f34ef77945b3',
    liveObservations: 172},
  validation: {
    fresh: 'round3 proof replay, original A/B and round2 A/B controls, 12 round3 controls, 172 live owner observations, pinned fmt, targeted Clippy, Node syntax, Typst HTML/PDF and docs registry',
    reused: 'Frozen source43/current producer captures and five qualified diagnostic ELFs; prior independent math and 336 original artifact bindings',
    controls: controls.map(({name, positive}) => ({name, positive})),
    oldAcceptNewReject: ['mirrored-input-card-missing-as', 'mirrored-input-card-as',
      'mirrored-dot-global-interface', 'mirrored-tsv-global-prefactor',
      'mirrored-tsv-global-prefactor-ranked'],
    originalAssertions: 'Both exit 101 at first grouped row; pristine nonfactorizable row not reached.',
    failedAttempts: [
      'Frozen round2 checker accepted the five A counterexamples that the corrected checker rejects; exact old/new outcomes are retained in a-round2-old.log and a-round2.log.',
      'The first round3 proof command used an obsolete control-output path and stopped; check-proof.log and check-proof.sh retain that attempt. check-proof-final.sh completed.',
      'An initial probe build produced the same source with the earlier command context; build.log and build-final.log retain both runs. No failed attempt is presented as a passing validation.',
    ],
  },
  dispositions: {
    'A-R2-1/B3': 'TSV global prefactor, DOT num/projector and both global interface declarations admitted through existing native tensor owner, each bound to the unit closed native global factors within each producer.',
    'A-R2-2': 'All 15 external name/complex-value pairs checked against imported Model after default restriction; omissions, extras and changes reject; internal ZERO excluded.',
    'A-R2-3': 'IoError and JsonError module imports and typed error sources corrected; rebuilt actual probe ELF and source closure qualified.',
  },
  artifacts,
  limits: ['No expectation, grouping algorithm, physics, lock, kernel or NixCI edit.',
    'No full pristine vacuum test pass, historical-kernel replay, private cache history or general grouping claim.',
    'No migration before new independent Astra review.',
    'Archived diagnostic formatting failures remain historical failures.'],
};
fs.writeFileSync(destination, JSON.stringify(result, null, 2) + '\n');
console.log(`round3 receipt: ${retained.length} round2 bindings retained; ` +
  `${result.round2Accounting.archived} archived revisions; ${controls.length} controls`);
