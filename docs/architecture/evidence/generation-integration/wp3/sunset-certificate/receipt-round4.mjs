import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';

const [destination] = process.argv.slice(2);
if (!destination) throw new Error('usage: receipt-round4.mjs <receipt-round4.json>');
const scripts = path.dirname(fileURLToPath(import.meta.url));
const workspace = path.resolve(scripts, '../../../../../..');
const root = path.join(workspace, 'target/sunset-certificate/round4');
const read = (file) => JSON.parse(fs.readFileSync(file, 'utf8'));
const sha256 = (bytes) => crypto.createHash('sha256').update(bytes).digest('hex');
const artifact = (file) => ({path: path.resolve(file), bytes: fs.statSync(file).size,
  sha256: sha256(fs.readFileSync(file))});
const output = (name) => artifact(path.join(root, name));
const old = ['frozen-review-heldout/results.json', 'frozen-review-reproduce/results.json']
  .flatMap((name) => read(path.join(root, name)));
const corrected = ['review-heldout/results.json', 'review-reproduce/results.json']
  .flatMap((name) => read(path.join(root, name)));
assert.equal(old.length, 8);
assert.equal(corrected.length, 8);
for (const [index, prior] of old.entries()) {
  const now = corrected[index];
  assert.equal(prior.name, now.name);
  assert.equal(prior.accepted, true, `${prior.name}: frozen round3 accepted`);
  assert.equal(now.accepted, false, `${now.name}: corrected checker rejects`);
  assert.equal(prior.shaGuards, false);
  assert.equal(now.shaGuards, false);
  assert.equal(now.checks.find((check) => check.stage === 'join').exit, 1);
}
const extra = read(path.join(root, 'round4-controls.json')).results;
assert.equal(extra.length, 8);
assert.equal(extra.filter((test) => test.positive).length, 4);
const observations = read(path.join(root, 'dot-observations.json'));
assert.equal(observations.observations.length, 16);
const closure = read(path.join(root, 'source-closure.json'));
assert.equal(closure.source.sha256,
  artifact(path.join(workspace, 'crates/feynkit-generator/examples/sunset_live_rule_probe.rs')).sha256);
assert.equal(closure.binaries[0].sha256, closure.binaries[1].sha256);
assert.equal(closure.features, '["default", "native"]');
assert.equal(closure.dependencies.find((dependency) => dependency.name === 'symbolica').kernel,
  '1a1829b');
assert(closure.dependencies.some((dependency) => dependency.name === 'linnet'));
const guard = read(path.join(root, 'guard-round4.json'));
assert.equal(guard.checkedFiles, 10289);
assert.equal(guard.originalBindings, 336);
assert.equal(guard.round2Bindings, 743);
assert.equal(guard.round3Bindings, 1096);
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
const scriptsInventory = fs.readdirSync(scripts).filter((name) => name.endsWith('.mjs'))
  .sort().map((name) => artifact(path.join(scripts, name)));
const proof = ['inputs.json', 'model.json', 'join.json', 'frame.json', 'internal-frame.json',
  'scalar-differential.json', 'mutations.json', 'internal-mutations.json',
  'chain-controls.json', 'correction-controls.json', 'round4-controls.json', 'positive.json',
  'review-heldout/results.json', 'review-reproduce/results.json',
  'frozen-review-heldout/results.json', 'frozen-review-reproduce/results.json',
  'dot-observations.json', 'source-closure.json', 'archive-round3.json',
  'guard-round4.json', 'title-audit.json'].map(output);
const commands = ['check-proof.sh', 'docs-check.sh', 'archive-round3.mjs',
  'guard-round4.mjs', 'dot-observations.mjs', 'source-closure.mjs',
  'review-heldout.mjs', 'review-reproduce.mjs',
  'frozen-review-heldout.mjs', 'frozen-review-reproduce.mjs'].map(output);
const logs = ['build.log', 'build-final.log', 'rust-final.log', 'rustfmt-apply.log',
  'check-proof.log', 'check-proof-final.log', 'review-heldout.log', 'review-reproduce.log',
  'frozen-review-heldout.log', 'frozen-review-reproduce.log', 'source-closure.log',
  'guard-round4.log', 'docs-check.log'].map(output);
const changed = execFileSync('jj', ['diff', '--name-only'],
  {cwd: workspace, encoding: 'utf8'}).trim().split('\n').filter(Boolean);
const nonReceiptPaths = changed.filter((name) => !name.endsWith('/receipt-round4.json'));
const nonReceiptStats = execFileSync('jj', ['diff', '--stat', '--', ...nonReceiptPaths],
  {cwd: workspace, encoding: 'utf8'}).trim();
const result = {
  schema: 4, status: 'bounded correction candidate; STOP for new independent Astra review',
  owner: {name: 'Lucien Huber', email: 'im@lcnbr.ch'},
  workspace, changeId: 'zppnyxunmykstkproknukopwuqyxpzyo',
  bookmark: 'generation-sunset-certificate-round4',
  frozenRound3: '30d4906436d9f62db22b621cd717c91f6a845d2e',
  priorReceipts: ['receipt.json', 'receipt-round2.json', 'receipt-round3.json']
    .map((name) => artifact(path.join(scripts, name))),
  guard: output('guard-round4.json'), preservedReferences: guard.previousRefs,
  changed, nonReceiptStats, sourceScripts: scriptsInventory,
  probe: {source: closure.source, manifest: closure.manifest, binary: closure.binaries[1],
    oldSource: output('sunset-live-rule-probe-round3.rs'),
    oldBinary: output('binaries/sunset-live-rule-probe-round3'),
    sourceClosure: output('source-closure.json'), features: closure.features,
    currentKernel: '1a1829b1f3d0a88b53b44085b3b9f34ef77945b3',
    observations: {sourceFinal: 40, comparisonRepresentative: 132, total: 172}},
  proof, commands, logs,
  documentation: ['docs/architecture/generation-sunset-certificate.typ',
    'docs/developers.toml'].map((name) => artifact(path.join(workspace, name))),
  renders: [output('note.html'), output('note.pdf')],
  corrections: {
    effectiveDOT: 'Linnet DotGraph::from_file supplies complete effective graph, vertex and edge attributes. Unit closed native admission binds the effective global num, projector and interfaces to TSV/native owners.',
    round3Fixture: 'Four round3 controls now place their same declarations inside digraph; old source is archived. The positive and ten negative intentions are unchanged.',
    oldAcceptNewReject: old.map((item) => item.name),
    newControls: extra.map(({name, positive}) => ({name, positive})),
    parserObservations: 16,
  },
  validation: {
    fresh: 'Ten proof replays, 22 outer and eight internal mutations, round2/round3/round4 controls, eight reviewer old/new differentials, sixteen production parser observations, 172 live-rule observations, pinned fmt/targeted Clippy, Node syntax, scoped Typst/title/registry and complete previous-file guard.',
    reused: 'Frozen source43/current producer captures and original five ELFs; no producer rebuild. Prior independent physical math and original assertion exits reused after guard.',
    attempts: [
      'Initial diagnostic build used node_data and omitted the DotGraph type annotation; the first build.log and command output retain those failures. The final build and Clippy pass are recorded separately.',
      'Initial complete proof stopped on a round3 control that placed declarations before digraph; check-proof.log retains failure. The same declarations now occur within the graph.',
      'First full semantic replay reached all controls but stopped at pinned rustfmt; rustfmt.log and rustfmt-apply.log retain the failure and correction. The final source was rebuilt and all checks rerun.',
      'First archive guard missed the intentionally saved old probe ELF; the corrected archive manifest and successful complete guard are retained.',
    ],
    originalAssertions: 'Both unchanged test functions exit 101 at the first grouped label assertion and never reach the pristine nonfactorizable assertion.',
  },
  limits: ['No grouped expectation, snapshot, production grouping, lock, kernel, NixCI or main integration edit.',
    'Current-kernel replay of source43 captures is not historical-kernel execution.',
    'No whole-graph numerator expansion, private cache history, general physics or whole pristine-test pass claim.',
    'No migration before a new independent Astra correction review.'],
};
fs.writeFileSync(destination, JSON.stringify(result, null, 2) + '\n');
console.log(`round4 receipt: ${old.length} old accepts/new rejects, ` +
  `${result.probe.observations.total} live-rule observations, ${guard.checkedFiles} guarded hashes`);
