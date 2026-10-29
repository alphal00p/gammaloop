import assert from 'node:assert/strict';
import crypto from 'node:crypto';
import fs from 'node:fs';
import path from 'node:path';

const workspace = process.cwd();
const migration = 'target/sunset-certificate/migration';
const correction = `${migration}/full-command-correction`;
const output = 'docs/architecture/evidence/generation-integration/wp3/sunset-certificate/migration-receipt.json';
const sha256 = (file) => crypto.createHash('sha256').update(fs.readFileSync(file)).digest('hex');
const artifact = (file) => ({
  path: path.resolve(file),
  bytes: fs.statSync(file).size,
  sha256: sha256(file),
});
const exit = (name) => Number(fs.readFileSync(`${migration}/${name}.exit`, 'utf8').match(/\d+/)[0]);
const correctedFull = JSON.parse(fs.readFileSync(`${correction}/commands-results.json`, 'utf8'));
assert.equal(correctedFull.results.length, 2);
for (const result of correctedFull.results) {
  assert.equal(result.exit, 101, `${result.name}: full-function exit`);
  assert.equal(result.selected, 1, `${result.name}: selected test count`);
  assert.equal(result.failed, 1, `${result.name}: failed test count`);
  assert.equal(result.passed, 0, `${result.name}: passed test count`);
  assert(!result.command.includes('--ignored'), `${result.name}: zero-selection --ignored flag`);
  const log = fs.readFileSync(result.log.path, 'utf8');
  assert.match(log, /running 1 test/);
  assert.match(log, /test result: FAILED\. 0 passed; 1 failed;/);
  assert.doesNotMatch(log, /running 0 tests/);
}
const [oldFull, finalFull] = correctedFull.results;
assert.equal(oldFull.name, 'old-full-corrected');
assert.equal(finalFull.name, 'final-full-corrected');
assert.match(fs.readFileSync(oldFull.log.path, 'utf8'), /Snapshot: vacuum_amplitude_generation-2/);
assert.match(fs.readFileSync(finalFull.log.path, 'utf8'), /failed in line 2095/);
const testSource = 'tests/tests/test_feyngen.rs';
const oldSource = `${migration}/test_feyngen-round5.rs`;
const changes = [
  {
    cell: 'all grouped two-loop',
    request: '{} > {} | g ghg t u d QED==0 [{2}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling',
    old: '5 | -19/24+Group(0,1,-1/2)+Group(5,1,-1/2) = -43/24',
    migrated: '5 | -19/24+Group(2,1,-1/2)+Group(4,1,-1/2) = -43/24',
    oldResult: {exit: exit('all_old'), executed: 1, status: 'genuine assertion failure'},
    migratedResult: {exit: exit('all_new'), executed: 1, status: 'pass'},
  },
  {
    cell: 'nonfactorizable grouped two-loop',
    request: '{} > {} | g ghg t u d QED==0 [{2}] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1 1',
    old: '4 | 1/12+Group(0,1,-1/2)+Group(4,1,-1/2) = -11/12',
    migrated: '4 | -11/12+Group(2,1,-1/2)+Group(4,1,-1/2) = -23/12',
    oldResult: {exit: exit('nonfactorizable_old'), executed: 1, status: 'genuine assertion failure'},
    migratedResult: {exit: exit('nonfactorizable_new'), executed: 1, status: 'pass'},
  },
];
assert.deepEqual(changes.map((entry) => [entry.oldResult.exit, entry.migratedResult.exit]), [[101, 0], [101, 0]]);
let expected = fs.readFileSync(oldSource, 'utf8');
for (const change of changes) {
  assert.equal(expected.split(change.old).length, 2, `old literal is not unique: ${change.cell}`);
  expected = expected.replace(change.old, change.migrated);
}
assert.equal(fs.readFileSync(testSource, 'utf8'), expected, 'test source is not an exact two-literal migration');

const logs = [
  'old-full.log', 'final-full.log', 'all_old.log', 'all_new.log',
  'nonfactorizable_old.log', 'nonfactorizable_new.log', 'final-build.log',
  'harness-build.log', 'harness-rustfmt.log', 'harness-static.log',
  'certificate-replay.log', 'guard.log', 'docs-final.log',
].map((name) => artifact(`${migration}/${name}`));
const frozenReceipts = [
  'receipt.json', 'receipt-round2.json', 'receipt-round3.json',
  'receipt-round4.json', 'receipt-round5.json',
].map((name) => artifact(`docs/architecture/evidence/generation-integration/wp3/sunset-certificate/${name}`));
const receipt = {
  schema: 2,
  status: 'full-command receipt correction candidate; STOP for new independent focused Astra review',
  owner: {name: 'Lucien Huber', email: 'im@lcnbr.ch'},
  workspace,
  changeId: 'wonnsqrpmkoxooolznvrkotwxnsurluo',
  bookmark: 'generation-sunset-migration',
  partialReference: 'generation-sunset-migration-partial-resume',
  authorization: {
    approvedCertificate: '9ede0c12fd0bd70458d10eeaaf636cb51c39543a',
    review: artifact('target/sunset-review-combined-round5-NGQS7n/review.json'),
    commands: artifact('target/sunset-review-combined-round5-NGQS7n/commands.json'),
    verdict: 'APPROVE_CERTIFICATE',
    source43: '74fcd8e57e9fe256d4bc34b5d8ad93c97fc28848',
    source43Kernel: 'c78cffc20f318df1fce50539a587cbe8ffc717c0',
    currentPhysicalBase: '59b57cead76e576fa9f23e3996f5614de939fe35',
    currentKernel: '1a1829b1f3d0a88b53b44085b3b9f34ef77945b3',
    migrationRound1: 'a202e83bc94954d4d1f900b4518febce7338b2f1',
    focusedReview: artifact('target/sunset-migration-astra-review-M9ve1U/review.json'),
    focusedReviewCommands: artifact('target/sunset-migration-astra-review-M9ve1U/commands.json'),
    finding: 'M1: --ignored selected zero tests in both recorded full-function command strings',
  },
  source: {
    old: artifact(oldSource),
    migrated: artifact(testSource),
    exactTwoLiteralDelta: true,
    changes,
    untouched: {
      ungroupedSource44Rows: 'all five remain at their approved-parent literals; not imported here',
      otherGroupedRows: 'unchanged',
      requestsModelsCardsWorkersSeedSamplesFiltersCountsProduction: 'unchanged',
      physicalCode: 'unchanged',
    },
  },
  builds: {
    profile: 'ci-optim',
    package: 'gammaloop-integration-tests@0.1.0',
    features: ['default'],
    normalSharedCache: true,
    oldProducer: artifact(`${migration}/binaries/current-diagnostic-test_feyngen-round5`),
    finalProducer: artifact(`${migration}/binaries/final-test_feyngen`),
    scopedProducer: artifact(`${migration}/binaries/scoped-test_feyngen`),
    finalBinding: artifact(`${migration}/final-build-binding.json`),
    scopedBinding: artifact(`${migration}/harness-build-binding.json`),
    scopedSource: artifact(`${migration}/harness-scoped.rs`),
    commands: {
      final: 'cargo test --locked -p gammaloop-integration-tests --test test_feyngen --profile ci-optim --no-run --message-format=json-render-diagnostics',
      scoped: 'cargo test --locked -p gammaloop-integration-tests --test sunset_migration_scoped --profile ci-optim --no-run --message-format=json-render-diagnostics',
      environment: 'nix develop; privately preserved or loaded SYMBOLICA_LICENSE; normal shared target/cache',
    },
  },
  validation: {
    recording: 'INSTA_UPDATE=no',
    commands: {
      oldFull: oldFull.command,
      finalFull: finalFull.command,
      scoped: 'run-scoped.sh invokes scoped-test_feyngen migration_scope::{all_old,all_new,nonfactorizable_old,nonfactorizable_new} --exact --nocapture',
    },
    fullCommandCorrection: {
      results: correctedFull.results,
      binding: artifact(`${correction}/commands-results.json`),
      runner: artifact(`${correction}/run-full.mjs`),
      runnerLog: artifact(`${correction}/run-full.log`),
      focusedGuard: artifact(`${correction}/guard.json`),
      focusedGuardScript: artifact(`${correction}/guard.mjs`),
      invariant: 'each full command must select exactly one test, fail exactly one test, exit 101, and omit --ignored',
    },
    negativeControl: {
      qualification: 'round-1 --ignored commands selected zero tests and receive no pass or assertion-failure credit',
      preservedRound1: 'generation-sunset-migration-round1 at a202e83bc94954d4d1f900b4518febce7338b2f1',
      reviewOldLog: artifact('target/sunset-migration-astra-review-M9ve1U/old-full.log'),
      reviewFinalStdout: artifact('target/sunset-migration-astra-review-M9ve1U/final-full-receipt-command.stdout'),
      reviewResults: artifact('target/sunset-migration-astra-review-M9ve1U/full-test-results.json'),
    },
    oldFull: {
      exit: oldFull.exit,
      selected: oldFull.selected,
      failed: oldFull.failed,
      reached: 'original all-grouped two-loop cell',
      status: 'expected assertion failure',
    },
    finalFull: {
      exit: finalFull.exit,
      selected: finalFull.selected,
      failed: finalFull.failed,
      reached: 'migrated all-grouped two-loop cell passed; stopped at unchanged three-loop grouped snapshot on test_feyngen.rs:2095',
      limitation: 'did not reach either nonfactorizable two-loop row and is not a whole-test pass',
    },
    scoped: {
      sourceMode: 'verbatim include of final test_feyngen.rs plus four narrow tests using existing feyngen_str',
      setup: 'sm.json import; 16 feyngen workers; original helper generation command and policies',
      interpretation: 'the second-row scoped pass proves its exact assertion only, not pristine full-function reachability',
    },
    documentation: [
      artifact(`${correction}/note.pdf`),
      artifact(`${correction}/note.html`),
      artifact(`${correction}/title-audit.json`),
      artifact(`${correction}/docs.log`),
    ],
    logs,
  },
  certifiedBinding: {
    members: [{sourceOccurrence: 2, ratio: 1}, {sourceOccurrence: 4, ratio: 1}],
    master: 'GL4',
    weights: ['-1/2', '-1/2'],
    all: {raw: 6, final: 5, total: '-43/24'},
    nonfactorizable: {raw: 5, final: 4, total: '-23/12'},
    numeratorExpansion: 'none; graph numerator remains factorized',
    proofReplay: artifact(`${migration}/certificate-replay.json`),
    guard: artifact(`${migration}/guard.json`),
    reuseQualification: 'certificate math and captures reused only after exact final-source/hash/producer and preservation binding',
    historicalQualification: 'current-kernel replay of historical captures is not historical-kernel execution',
  },
  checks: {
    fresh: [
      'corrected old/final exact original named test executions with recording disabled; one selected failure and exit 101 each',
      'four scoped old-fail/new-pass executions',
      'scoped rustfmt check',
      'scoped Cargo check/Clippy compilation (two unrelated vendored Clarabel warnings)',
      'certificate replay and frozen-artifact guard',
    ],
    reused: 'approved certificate proof and retained source43/current captures after hash verification',
  },
  preserved: {
    references: {
      'generation-sunset-certificate': '6732cfbb59099e23c315665440dbdae67a3944d4',
      'generation-sunset-certificate-round1': '6732cfbb59099e23c315665440dbdae67a3944d4',
      'generation-sunset-certificate-round2': '6aa6d311432d64b68e9208d3db0326aafbb080a3',
      'generation-sunset-certificate-round3': '30d4906436d9f62db22b621cd717c91f6a845d2e',
      'generation-sunset-certificate-resumed-partial': '5a8c974714077056f85f0d251ddede2b818c49b5',
      'generation-sunset-certificate-round4': '19d2ef715d05bec8b1c4ec09c7583cb710346e9f',
      'generation-sunset-certificate-approved': '9ede0c12fd0bd70458d10eeaaf636cb51c39543a',
      'generation-sunset-certificate-round5': '9ede0c12fd0bd70458d10eeaaf636cb51c39543a',
      'generation-sunset-historical-source43': '74fcd8e57e9fe256d4bc34b5d8ad93c97fc28848',
      'generation-sunset-migration-partial-resume': '6653cdf7c9e697e8fe6d538e8abb1ba4a40f636b',
      'generation-sunset-migration-round1': 'a202e83bc94954d4d1f900b4518febce7338b2f1',
    },
    receipts: frozenReceipts,
    protectedFiles: ['Cargo.lock', 'flake.lock', 'nix-ci.nix'],
    artifacts: 'all certificate captures, proof outputs, probe ELFs, review directories, migration sources, ELFs, commands and logs retained',
  },
};
fs.writeFileSync(output, `${JSON.stringify(receipt, null, 2)}\n`);
console.log(`wrote ${output}`);
