import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';

const [destination] = process.argv.slice(2);
if (!destination) throw new Error('usage: receipt-round5.mjs <receipt-round5.json>');
const scripts = path.dirname(fileURLToPath(import.meta.url));
const workspace = path.resolve(scripts, '../../../../../..');
const root = path.join(workspace, 'target/sunset-certificate/round5');
const read = (file) => JSON.parse(fs.readFileSync(file, 'utf8'));
const sha256 = (file) => crypto.createHash('sha256').update(fs.readFileSync(file)).digest('hex');
const artifact = (file) => ({path: path.resolve(file), bytes: fs.statSync(file).size,
  sha256: sha256(file)});
const output = (name) => artifact(path.join(root, name));
const old = read(path.join(root, 'frozen-seven/results.json'));
const corrected = read(path.join(root, 'review-seven/results.json'));
assert.equal(old.length, 7);
assert.equal(corrected.length, 7);
for (const [index, prior] of old.entries()) {
  const now = corrected[index];
  assert.equal(prior.name, now.name);
  assert.equal(prior.accepted, true, `${prior.name}: frozen round4 accepted`);
  assert.equal(now.accepted, false, `${now.name}: corrected checker rejects`);
  assert.equal(prior.shaGuards, false);
  assert.equal(now.shaGuards, false);
  assert.equal(now.checks.find((check) => check.stage === 'join').exit, 1);
}
const controls = read(path.join(root, 'round5-controls.json')).results;
assert.equal(controls.length, 6);
assert.equal(controls.filter((test) => test.positive).length, 2);
const owners = read(path.join(root, 'local-owner-observations.json'));
assert.equal(owners.observations.length, 14);
const closure = read(path.join(root, 'source-closure.json'));
assert.equal(closure.source.sha256, sha256(path.join(workspace,
  'crates/feynkit-generator/examples/sunset_live_rule_probe.rs')));
assert.equal(closure.binaries[0].sha256, closure.binaries[1].sha256);
assert.equal(closure.features, '["default", "native"]');
assert.equal(closure.dependencies.find((dependency) => dependency.name === 'symbolica').kernel,
  '1a1829b');
assert(closure.dependencies.some((dependency) => dependency.name === 'linnet'));
const guard = read(path.join(root, 'guard-round5.json'));
assert.equal(guard.checkedFiles, 11606);
assert.equal(guard.originalBindings, 336);
assert.equal(guard.round2Bindings, 743);
assert.equal(guard.round3Bindings, 1096);
assert.equal(guard.round4Bindings, 79);
const join = read(path.join(root, 'join.json'));
const internal = read(path.join(root, 'internal-frame.json'));
for (const route of ['source43', 'current']) {
  for (const [name, raw, final, total, singleton] of [
    ['all', 6, 5, '-43/24', '-19/24'],
    ['nonfactorizable', 5, 4, '-23/12', '-11/12'],
  ]) {
    const result = join.native[route][name];
    assert.equal(result.sourceCount, raw);
    assert.equal(result.finalCount, final);
    assert.equal(result.rawTotal, total);
    assert.equal(result.groupedTotal, total);
    assert.equal(result.singletonConstant, singleton);
  }
  assert(join.liveRules[route].result.includes('passed: 20 diagrams'));
  assert(internal.liveRules[route].result.includes('passed: 66 diagrams'));
}
const sourceScripts = fs.readdirSync(scripts).filter((name) => name.endsWith('.mjs'))
  .sort().map((name) => artifact(path.join(scripts, name)));
const proof = ['inputs.json', 'model.json', 'join.json', 'frame.json',
  'internal-frame.json', 'scalar-differential.json', 'mutations.json',
  'internal-mutations.json', 'chain-controls.json', 'correction-controls.json',
  'round4-controls.json', 'round5-controls.json', 'positive.json',
  'review-seven/results.json', 'frozen-seven/results.json',
  'local-owner-observations.json', 'source-closure.json',
  'archive-round4.json', 'guard-round5.json', 'title-audit.json'].map(output);
const commands = ['check-proof.sh', 'docs-check.sh', 'source-closure.mjs',
  'archive-round4.mjs', 'guard-round5.mjs', 'local-owner-observations.mjs',
  'frozen-reproduce-round4.mjs'].map(output);
const logs = ['build.log', 'check-proof.log', 'clippy.log', 'rustfmt.log',
  'review-seven.log', 'frozen-seven.log', 'round5-controls.log',
  'source-closure.log', 'archive-round4.log', 'guard-round5.log',
  'docs-check.log'].map(output);
const changed = execFileSync('jj', ['diff', '--name-only'],
  {cwd: workspace, encoding: 'utf8'}).trim().split('\n').filter(Boolean);
const nonReceiptPaths = changed.filter((name) => !name.endsWith('/receipt-round5.json'));
const nonReceiptStats = execFileSync('jj', ['diff', '--stat', '--', ...nonReceiptPaths],
  {cwd: workspace, encoding: 'utf8'}).trim();
const result = {
  schema: 5, status: 'bounded correction candidate; STOP for new independent Astra review',
  owner: {name: 'Lucien Huber', email: 'im@lcnbr.ch'},
  workspace, changeId: 'zrlnsssqvwtpwvtltplxxxxpyyrzylnx',
  bookmark: 'generation-sunset-certificate-round5',
  frozenRound4: '19d2ef715d05bec8b1c4ec09c7583cb710346e9f',
  priorReceipts: ['receipt.json', 'receipt-round2.json', 'receipt-round3.json',
    'receipt-round4.json'].map((name) => artifact(path.join(scripts, name))),
  review: [artifact(path.join(workspace, 'target/sunset-review-combined-round4-j3q756/review.json')),
    artifact(path.join(workspace, 'target/sunset-review-combined-round4-j3q756/reproduce-round4.mjs'))],
  guard: output('guard-round5.json'), preservedReferences: guard.previousRefs,
  changed, nonReceiptStats, sourceScripts,
  probe: {source: closure.source, manifest: closure.manifest, binary: closure.binaries[1],
    oldSource: output('sunset-live-rule-probe-round4.rs'),
    oldBinary: output('binaries/sunset-live-rule-probe-round4'),
    sourceClosure: output('source-closure.json'), features: closure.features,
    currentKernel: '1a1829b1f3d0a88b53b44085b3b9f34ef77945b3',
    observations: {sourceFinal: 40, comparisonRepresentative: 132, total: 172}},
  proof, commands, logs,
  documentation: ['docs/architecture/generation-sunset-certificate.typ',
    'docs/developers.toml'].map((name) => artifact(path.join(workspace, name))),
  renders: [output('note.html'), output('note.pdf')],
  corrections: {
    parsedOwners: 'Linnet parsed HedgePair source/sink and node-owner IDs bind the actual DOT edges to native half_edge_order. Display port labels are retained only as labels.',
    parsedInterfaces: 'The same Rust NumeratorInterface::parse handles each effective global/local declaration exactly once; absent metadata is distinct from present valid arrays, and present malformed/null declarations reject.',
    oldAcceptNewReject: old.map((item) => item.name),
    positiveControls: controls.filter((test) => test.positive).map((test) => test.name),
    extraNegativeControls: controls.filter((test) => !test.positive).map((test) => test.name),
    directOwnerObservations: owners.observations.length,
  },
  validation: {
    fresh: 'Seven exact SHA-free old/new full-pipeline differentials, fourteen direct production parser observations, ten proof replays, 22 outer/eight internal mutations, round2/round3/round4/round5 controls, 172 live-rule observations, pinned fmt/Clippy, Node syntax, scoped Typst/title/registry, complete prior-file guard.',
    reused: 'Frozen original source43/current producer captures and five ELFs, prior independent math and original assertion controls. No generation producer rerun or historical-kernel replay.',
    attempts: ['The first direct parser-observation invocation lacked SYMBOLICA_LICENSE and could not parse its stdout. The licensed run completed; the diagnostic source and seven-case proof were unchanged.',
      'The initial observation expectation treated an unmodified historical extra-escape specimen as mutated. The corrected fourteen-file inventory records each route separately.'],
    originalAssertions: 'Both unchanged test functions exit 101 at the first grouped label assertion; the pristine nonfactorizable assertion is not reached.',
  },
  limits: ['No grouped expectation, snapshot, production grouping, lock, kernel, NixCI or main integration edit.',
    'Source43 capture replay on the current kernel is not historical-kernel execution.',
    'No whole-graph numerator expansion, hidden HEP/Dirac preparation, private cache history, general physics or whole pristine-test pass claim.',
    'No migration before a new independent Astra correction review.'],
};
fs.writeFileSync(destination, JSON.stringify(result, null, 2) + '\n');
console.log(`round5 receipt: ${old.length} old accepts/new rejects, ` +
  `${result.probe.observations.total} live-rule observations, ${guard.checkedFiles} guarded hashes`);
