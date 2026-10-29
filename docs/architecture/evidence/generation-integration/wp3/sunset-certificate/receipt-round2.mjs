import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';

const [destination] = process.argv.slice(2);
if (!destination) throw new Error('usage: receipt-round2.mjs <receipt-round2.json>');
const scripts = path.dirname(fileURLToPath(import.meta.url));
const workspace = path.resolve(scripts, '../../../../../..');
const root = path.join(workspace, 'target/sunset-certificate/round2');
const read = (file) => JSON.parse(fs.readFileSync(file, 'utf8'));
const hash = (bytes) => crypto.createHash('sha256').update(bytes).digest('hex');
const artifact = (file) => ({path: path.resolve(file), bytes: fs.statSync(file).size,
  sha256: hash(fs.readFileSync(file))});
const oldReceipt = path.join(scripts, 'receipt.json');
const original = read(oldReceipt);
const originalBindings = [];
function bindings(value) {
  if (!value || typeof value !== 'object') return;
  if (value.path && value.sha256) originalBindings.push(value);
  else Object.values(value).forEach(bindings);
}
bindings(original);
assert.equal(originalBindings.length, 336);
const retained = originalBindings.map((entry) => {
  const live = artifact(entry.path);
  if (live.sha256 === entry.sha256) return {original: entry, retained: live, location: 'unchanged'};
  const relative = path.relative(workspace, entry.path);
  assert(!relative.startsWith('..'), `historical artifact changed: ${relative}`);
  const archive = artifact(path.join(root, 'round1-original-files', relative));
  assert.equal(archive.sha256, entry.sha256, `original artifact lost: ${relative}`);
  return {original: entry, retained: archive, location: 'original source archived before revision'};
});
assert.equal(retained.filter((entry) => entry.location !== 'unchanged').length, 9);
const protectedPaths = ['Cargo.lock', 'flake.lock', 'nix-ci.nix'];
const protectedFiles = protectedPaths.map((file) => {
  const expected = execFileSync('jj', ['file', 'show', '-r', '6732cfbb', file], {cwd: workspace});
  const actual = artifact(path.join(workspace, file));
  assert.equal(actual.sha256, hash(expected), `${file}: protected bytes changed`);
  return actual;
});
const jj = (...args) => execFileSync('jj', ['--ignore-working-copy', ...args],
  {cwd: workspace, encoding: 'utf8'}).trim();
const refs = Object.fromEntries(['generation-sunset-certificate',
  'generation-sunset-certificate-round1', 'generation-sunset-historical-source43',
  'generation-sunset-certificate-resumed-partial'].map((name) =>
  [name, jj('log', '-r', name, '--no-graph', '-T', 'commit_id')]));
assert.equal(refs['generation-sunset-certificate'], '6732cfbb59099e23c315665440dbdae67a3944d4');
assert.equal(refs['generation-sunset-certificate-round1'], refs['generation-sunset-certificate']);
assert.equal(refs['generation-sunset-historical-source43'], '74fcd8e57e9fe256d4bc34b5d8ad93c97fc28848');
assert.equal(refs['generation-sunset-certificate-resumed-partial'], '5a8c974714077056f85f0d251ddede2b818c49b5');
const resultAt = (file) => read(path.join(root, file));
const allPass = (checks) => checks.every((check) => check.exit === 0);
const bOld = resultAt('old-reviewer-b/independent-controls.json').results;
const bNew = resultAt('repaired/b-review-replay/independent-controls.json').results;
const comparisons = bOld.map((old) => {
  const fresh = bNew.find((entry) => entry.name === old.name);
  assert(fresh, `missing unchanged B control ${old.name}`);
  const positive = ['baseline', 'positive_metric_symmetric_swap',
    'positive_scalar_representative_equivalence'].includes(old.name);
  assert.equal(allPass(fresh.checks), positive, `B control result ${old.name}`);
  return {reviewer: 'B', name: old.name, positive,
    layoutProofLimit: old.name === 'comparison_interface_permutation',
    originalAccepted: allPass(old.checks), correctedAccepted: allPass(fresh.checks)};
});
const aOld = resultAt('old-reviewer-a/results.json');
const aNew = resultAt('repaired/a-provenance-last/results.json');
for (const old of aOld) {
  const fresh = aNew.find((entry) => entry.name === old.name);
  const positive = old.name === 'unchanged baseline';
  assert(fresh);
  assert.equal(allPass(fresh.results), positive, `A control ${old.name}`);
  assert(allPass(old.results));
  comparisons.push({reviewer: 'A', name: old.name, positive,
    originalAccepted: true, correctedAccepted: allPass(fresh.results)});
}
const ownerOld = resultAt('old-model-owner.log');
const ownerNew = resultAt('reviewer-a-model.log');
assert(allPass(ownerOld.results));
assert(!allPass(ownerNew.results));
comparisons.push({reviewer: 'A', name: 'mirrored wrong native/observed model owner',
  positive: false, originalAccepted: true, correctedAccepted: false});
const scalar = resultAt('scalar-differential.json');
assert.equal(scalar.accepted, 40);
assert.equal(scalar.explicitlyRejected, 17);
assert.equal(resultAt('mutations.json').results.length, 22);
assert.equal(resultAt('internal-mutations.json').results.length, 8);
assert.equal(resultAt('positive.json').results.length, 9);
const chain = resultAt('chain-controls.json').results;
assert.equal(chain.filter((entry) => entry.positive).length, 2);
assert.equal(chain.filter((entry) => !entry.positive).length, 4);
const joined = resultAt('join.json');
for (const route of ['source43', 'current']) {
  assert.equal(joined.native[route].all.groupedTotal, '-43/24');
  assert.equal(joined.native[route].nonfactorizable.groupedTotal, '-23/12');
  assert(joined.liveRules[route].result.includes('passed: 20 diagrams;'));
  assert(resultAt('internal-frame.json').liveRules[route].result.includes('passed: 66 diagrams;'));
}
const proofNames = ['inputs', 'model', 'join', 'frame', 'internal-frame', 'scalar-differential',
  'mutations', 'internal-mutations', 'chain-controls', 'positive', 'title-audit'];
const proofArtifacts = Object.fromEntries(proofNames.map((name) =>
  [name, artifact(path.join(root, `${name}.json`))]));
const reviewFiles = ['old-reviewer-b/independent-controls.json',
  'repaired/b-review-replay/independent-controls.json', 'old-reviewer-a/results.json',
  'repaired/a-provenance-last/results.json', 'old-model-owner.log', 'reviewer-a-model.log'];
const supplementalBinary = artifact(path.join(root, 'binaries/sunset-live-rule-probe'));
assert.equal(supplementalBinary.sha256,
  artifact(path.join(workspace, 'target/debug/examples/sunset_live_rule_probe')).sha256);
assert.equal(artifact(path.join(root, 'sunset-live-rule-probe.rs')).sha256,
  artifact(path.join(workspace, 'crates/feynkit-generator/examples/sunset_live_rule_probe.rs')).sha256);
const result = {
  schema: 2, status: 'candidate; STOP for new independent Astra A/B review before migration',
  owner: {name: 'Lucien Huber', email: 'im@lcnbr.ch', registry: 'lucien-huber'},
  workspace, changeId: 'vtvkkxko', bookmark: 'generation-sunset-certificate-round2',
  originalCandidate: refs['generation-sunset-certificate'], preservedReferences: refs,
  physicalParent: original.physicalParent, kernels: original.kernels,
  scope: 'Exactly the two grouped two-loop vacuum rows and their paired raw controls. No expectation or production algorithm edits.',
  originalReceipt: artifact(oldReceipt), originalBindingsRetained: retained,
  originalArtifactAccounting: {total: 336, unchangedPaths: 327, archivedOriginalSources: 9},
  protectedFiles, historicalTrackedGuard: artifact(path.join(root, 'historical-guard.json')),
  sourceScripts: fs.readdirSync(scripts).filter((name) => name.endsWith('.mjs')).sort()
    .map((name) => artifact(path.join(scripts, name))),
  supplementalOwner: {
    source: artifact(path.join(workspace, 'crates/feynkit-generator/examples/sunset_live_rule_probe.rs')),
    sourceArchive: artifact(path.join(root, 'sunset-live-rule-probe.rs')),
    manifest: artifact(path.join(workspace, 'crates/feynkit-generator/Cargo.toml')),
    binary: supplementalBinary,
    profile: 'dev; Cargo defaults; normal shared cache; no CARGO_TARGET_DIR override',
    modelOwner: 'Model::from_json -> canonical Model::fingerprint before diagram ID resolution',
    localOwner: 'FeynmanDiagramRulesExt::apply_feynman_rules -> Generator::construct_numerators',
    tensorOwner: 'SymbolicTensor::permuted and same_value; explicit named-slot bijections retained for inferred interfaces only',
    observations: {nativeSourceFinal: 40, sourceComparisonRepresentative: 132},
    presentationMap: 'Source43 Q tensor-label:q and edge/vertex index-label display tags -> current corresponding head; identity components/ports. No other alias.',
    aggregate: 'Complete factorized product of all fragments, with local sums undistributed; unit prefactor/projector and closed interface checked separately',
  },
  freshValidation: {
    exactProofReplays: 9, originalNegatives: 30,
    pinnedScalar: {acceptedMatches: 40, unsupportedRejected: 17, oracle: scalar.oracle,
      source: scalar.oracleSource},
    reviewComparisons: comparisons, shaGuardsInSemanticControls: false,
    additionalChainControls: chain.map(({name, positive}) => ({name, positive})),
    javascript: 'Pinned Node v24.19.0 --check all proof scripts; final-newline/trailing-whitespace/tab checks. No configured JavaScript formatter was found.',
    rust: 'Pinned cargo fmt --all -- --check and cargo clippy --locked -p feynkit-generator --example sunset_live_rule_probe -- -D warnings pass; explicit module-scope imports reviewed.',
    rustPolicyDeparture: 'The standalone probe uses a private typed thiserror result at main rather than adding eyre to feynkit-generator dev dependencies: adding that dependency would require a prohibited Cargo.lock edit. Existing typed errors are preserved for IO/model/diagram/generation failures.',
    documentation: 'Fresh Typst HTML/PDF, single source-title/canonical-title audit and alphal00p-docs-builder check; not the full five-product documentation suite.',
  },
  reusedValidation: {
    producerRuns: 'Original qualified five ELFs/captures; A freshly replayed diagnostic/frame producers and both original assertion failures. No producer rebuild for round2.',
    originalAssertionExit: 101, pristineNonfactorizableReached: false,
    archivedFormatting: 'A found rustfmt/import-policy failures in frozen diagnostic sources/patches. Old restored-worktree fmt pass does not cover them; archives were not rewritten.',
    independentMath: 'B independent-math evidence was checked as a reference; no copied-code claim of new independence.',
  },
  dispositions: {
    A1: 'CLI/native local expressions, coefficients, complete interfaces, graph indices, group IDs and master flags bound within each producer; final payload equals selected real source.',
    A2: 'Complete live model owner recomputes fingerprint before native/frame ID resolution; every diagram model checked; unchanged wrong-model negative rejects.',
    A3: 'Complete four resolved definitions and relevant post-import settings checked against intended policies and each other; unilateral and mirrored seed mutations reject.',
    B1: 'Bounded exact grammar, qualified identities, head arguments/occurrences, full consumption, pinned differential, undefined forms rejected.',
    B2: 'Unique physical raw source occurrences partition all members exactly once, including singleton coverage and actual master membership; all Group occurrences match.',
    B3: 'Full production local-rule replay plus factorized aggregate proof closes model/local/aggregate/selected-native/CLI chain. Mirrored ghost/gluon coefficients reach local owner and reject.',
    B4: 'Live endpoint slots, carrier/antiparticle and oriented map checked; named logical-port permutations retained and proved by tensor permutation owner; unsupported declared-order change rejects.',
  },
  proofArtifacts,
  reviewerHarnessSources: ['old-reviewer-a.mjs', 'old-reviewer-b.mjs', 'old-model-owner.mjs',
    'repaired/b-controls.mjs', 'repaired/a-model-owner.mjs',
    ...fs.readdirSync(path.join(root, 'original-checker')).map((file) => `original-checker/${file}`)]
    .map((file) => artifact(path.join(root, file))),
  reviewerArtifacts: reviewFiles.map((file) => artifact(path.join(root, file))),
  commands: ['run-checks.sh', 'verify-final.sh', 'docs-check.sh'].map((file) => artifact(path.join(root, file))),
  logs: ['verification-final.log', 'clippy-final.log', 'rustfmt-final.log', 'docs-check-final.log',
    'live-rule-build-clean.log', 'old-reviewer-a.log', 'old-reviewer-b.log',
    'live-rule-pair.log'].map((file) => artifact(path.join(root, file))),
  documentation: ['docs/architecture/generation-sunset-certificate.typ', 'docs/developers.toml']
    .map((file) => artifact(path.join(workspace, file))),
  renders: ['note.html', 'note.pdf'].map((file) => artifact(path.join(root, file))),
  limits: ['No expectation migration; new independent A/B review remains required.',
    'Finite exported production 4D/SU(3) preparation only; no graph expansion or hidden HEP/Dirac preparation.',
    'Unobserved private tensor-cache/collision histories are unclaimed.',
    'Unsupported declared-layout permutations are rejected without calling them necessarily physics-changing.',
    'All original raw captures, binaries, logs, proofs and refs preserved; nine revised source originals additionally archived.',
    'No main/other-workspace implementation edits, fetch/push/PR/upload/publication/CI enabling or lock/kernel change.'],
};
fs.writeFileSync(destination, JSON.stringify(result, null, 2) + '\n');
console.log('round2 receipt: 336 original identities retained; 30 original negatives; all reviewer escapes rejected; 172 live-rule observations');
