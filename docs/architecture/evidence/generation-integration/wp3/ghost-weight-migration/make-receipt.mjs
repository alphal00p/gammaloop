import fs from 'node:fs';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';

const evidence = 'docs/architecture/evidence/generation-integration/wp3/ghost-weight-migration';
const target = 'target/ghost-weight-certificate/migration';
const sha = (location) => crypto.createHash('sha256').update(fs.readFileSync(location)).digest('hex');
const file = (location) => ({path: location, sha256: sha(location)});
const source = JSON.parse(fs.readFileSync(`${target}/five-row-diff.json`));
const testResults = JSON.parse(fs.readFileSync(`${target}/test-results.json`));
const capture = JSON.parse(fs.readFileSync(`${target}/final-capture-comparison.json`));
const guard = JSON.parse(fs.readFileSync(`${target}/proof-guard.json`));
const harness = JSON.parse(fs.readFileSync(`${target}/vacuum-harness-source.json`));
const reviewPath = 'target/ghost-weight-certificate/review-astra-round3-combined-5LJwn8/review.json';
const review = JSON.parse(fs.readFileSync(reviewPath));
assert.equal(review.verdict, 'APPROVE_CERTIFICATE');
assert.equal(source.changedLineCount, 5);
assert.equal(capture.totalSelectedFiles, 4332);
assert.deepEqual(guard.frozenBindings, {round1: 199, round2: 40, approved: 44});
const rows = source.rows.map((row) => ({...row,
  result: testResults.cases[row.id] ?? testResults.cases[row.id.replace('vacuum_nonfactorizable_', 'vacuum_')]}));
assert(rows.every((row) => row.result));
const caseNames = [
  ['old-epem', 'slow::test_generate_epem_a_qqh', 'old-literals-test_feyngen'],
  ['new-epem', 'slow::test_generate_epem_a_qqh', 'new-literals-test_feyngen'],
  ['old-aa', 'slow::test_generate_aa_ttx_cross_section', 'old-literals-test_feyngen'],
  ['new-aa', 'slow::test_generate_aa_ttx_cross_section', 'new-literals-test_feyngen'],
  ...[2, 3, 4].flatMap((loop) => [
    [`vacuum-old-l${loop}`, `scoped_vacuum_l${loop}`, 'vacuum-old-assertions'],
    [`vacuum-new-l${loop}`, `scoped_vacuum_l${loop}`, 'vacuum-new-assertions'],
  ]),
  ['new-original-vacuum', 'slow::test_vacuum_amplitude_generation', 'new-literals-test_feyngen'],
  ['final-diagnostic', 'diagnose_source44_ungrouped_weights', 'new-literals-test_feyngen'],
];
const commands = caseNames.map(([id, test, executable]) => ({id,
  replayCommand: `nix develop -c env -u INSTA_FORCE_PASS INSTA_UPDATE=no TESTS_GAMMALOOP_STATE_PATH=<isolated-${id}-state> ${target}/${executable} ${id === 'final-diagnostic' ? '--ignored ' : ''}--exact ${test} --nocapture --test-threads=1`,
  actualTestArgv: [...(id === 'final-diagnostic' ? ['--ignored'] : []), '--exact', test,
    '--nocapture', '--test-threads=1'],
  executable: file(`${target}/${executable}`),
  log: file(`${target}/${id}.log`),
  environment: {INSTA_UPDATE: 'no', INSTA_FORCE_PASS: 'unset',
    TESTS_GAMMALOOP_STATE_PATH: `${target}/${id}-state`}}));
const receipt = {
  schema: 1,
  status: 'prepared five-literal migration; NEW independent Astra migration/test review required',
  author: 'Lucien Huber <im@lcnbr.ch>',
  parent: '121ac233d5ab014f4fadefccce702041bbbdf487',
  jjChange: 'qouoqyzw',
  bookmarks: {approvedCertificate: 'generation-ghost-weight-certificate-approved',
    certificate: 'generation-ghost-weight-certificate', migration: 'generation-ghost-weight-migration',
    approvedNamedAaMigration: 'generation-pocket-weight-diagnosis-approved-migration',
    priorCertificate: 'generation-pocket-weight-diagnosis-approved-certificate'},
  approval: file(reviewPath),
  receiptSources: {generator: file(`${evidence}/make-receipt.mjs`),
    validator: file(`${evidence}/validate-receipt.mjs`)},
  sourceAndFeatureIdentity: {
    approvedSourceRef: 'generation-ghost-weight-certificate-approved',
    approvedTestSourceSha256: source.originalSourceSha256,
    migratedTestSource: file('tests/tests/test_feyngen.rs'),
    frozenPhysicsParent: '59b57cead76e576fa9f23e3996f5614de939fe35',
    source43Historical: '74fcd8e57e9fe256d4bc34b5d8ad93c97fc28848',
    currentKernel: '1a1829b1f3d0a88b53b44085b3b9f34ef77945b3',
    source43Kernel: 'c78cffc20f318df1fce50539a587cbe8ffc717c0',
    cargoLock: file('Cargo.lock'), flakeLock: file('flake.lock'),
    model: file('assets/models/json/sm/sm.json'),
    cargoFeatures: 'Cargo defaults; ci-optim profile; no feature override',
    finalBuildCommand: 'nix develop -c cargo test --locked -p gammaloop-integration-tests --test test_feyngen --profile ci-optim --no-run -j 2',
    normalCargoProducer: '/common/build-cache/lcnbr/9c/5e5e2f7a370b0f/ci-optim/deps/test_feyngen-9581d43dfd98ca33',
    oldImmutableBinary: file(`${target}/old-literals-test_feyngen`),
    newImmutableBinary: file(`${target}/new-literals-test_feyngen`),
    finalCargoBuildLog: file(`${target}/final-cargo-build.log`),
    finalTargetedClippyLog: file(`${target}/final-clippy.log`),
    pinnedSourceChecksLog: file(`${target}/source-checks.log`),
  },
  literalMigration: {rows, diffSource: file(`${evidence}/check-five-row-diff.mjs`),
    exactDiff: file(`${target}/five-row-diff.json`), negativeInventory: source.negativeInventory},
  vacuumHarness: {
    attribution: 'task-local observational assertion harness, generated from exact old/new owner source and using existing scalar generation helpers; original grouped prefix unexecuted',
    generator: file(`${evidence}/make-vacuum-harness.mjs`),
    sourceBinding: file(`${target}/vacuum-harness-source.json`),
    oldSource: file(`${target}/vacuum-old-assertions.rs`),
    newSource: file(`${target}/vacuum-new-assertions.rs`),
    oldBinary: file(`${target}/vacuum-old-assertions`),
    newBinary: file(`${target}/vacuum-new-assertions`),
    buildScript: file(`${target}/build-vacuum-harness.sh`),
    successfulBuildLog: file(`${target}/vacuum-harness-build-attempt3.log`),
    sourceOwnerHelperHashes: harness.helperRanges,
    originalSetup: harness.old.setup,
    exactCommands: harness.old.assertions.map((row) => row.command),
    originalWorkerChoice: 16,
    currentPinnedSymbolicaRlib: file('/common/build-cache/lcnbr/9c/5e5e2f7a370b0f/ci-optim/deps/libsymbolica-219a1a53d9363ef7.rlib'),
  },
  testExecution: {
    verificationSource: file(`${evidence}/verify-test-logs.mjs`),
    structuredResults: file(`${target}/test-results.json`),
    commands,
    results: testResults,
    oldTargetFailures: 5,
    newTargetPasses: 5,
    completeNamedCrossSectionTestPasses: 0,
    originalVacuumWholeTestPasses: 0,
    laterUnchangedGroupedFailuresRetained: ['generate_epem_a_qqh-7',
      'generate_aa_ttx_cross_section-4', 'vacuum_amplitude_generation-2'],
    recordAndUpdateModes: {INSTA_UPDATE: 'no', INSTA_FORCE_PASS: 'unset', INSTA_RECORD: 'unset'},
  },
  physicalPreservation: {
    freshFinalBinaryCapture: file(`${target}/final-diagnostic.log`),
    comparisonSource: file(`${evidence}/verify-final-capture.mjs`),
    comparison: file(`${target}/final-capture-comparison.json`),
    counts: capture.counts,
    finalVsApprovedTreeSha256: capture.finalByteTreeSha256IgnoringSettingsFolder,
    proofGuardSource: file(`${evidence}/validate-proof.mjs`),
    proofGuard: file(`${target}/proof-guard.json`),
    approvedReviewArtifactsRehashed: guard.reviewedArtifacts,
    frozenBindingsRehashed: guard.frozenBindings,
    round1JoinSha256: guard.frozenFullJoinSha256[0],
    round2AndRound3JoinSha256: guard.frozenFullJoinSha256[1],
    modelAndResolvedParameters: {liveModelSha256: '7e250b3ff20c98277ee19af863cc030c7a5fd47fea979931b28aba744f8a9de3',
      resolvedCardSha256: '232964f2b9783482f8ccb0e8812eb6e34a0d2dd6d5e8a8905a5b98377ae7fe1d'},
  },
  documentation: {
    note: file('docs/architecture/generation-ghost-weight-migration.typ'),
    registry: file('docs/developers.toml'),
    pdf: file(`${target}/migration.pdf`), html: file(`${target}/migration.html`),
    scopedAudit: file(`${target}/migration-doc-audit.json`),
    scopedCheckLog: file(`${target}/migration-doc-check.log`),
    registryCheckLog: file(`${target}/docs-builder-check.log`),
  },
  retainedFailedAttempts: [
    {stage: 'direct scratch rustfmt outside Nix', result: 'rustfmt unavailable outside development shell; no source changed'},
    {stage: 'first scoped-harness Nix build', log: file(`${target}/vacuum-harness-build.log`),
      result: 'rustfmt import ordering only; corrected in scratch generator'},
    {stage: 'second scoped-harness Nix build', log: file(`${target}/vacuum-harness-build-attempt2.log`),
      result: 'scratch rustc lacked Cargo macro environment and one extern path typo; corrected in task-owned build script'},
  ],
  limitations: [
    'The named slow cross-section tests pass the migrated cells but stop at later unchanged grouped assertions; no whole-test pass is claimed.',
    'The original vacuum function stops at an earlier unchanged grouped assertion; isolated 16-worker harnesses validate only the three original ungrouped cells.',
    'Whole-file rustfmt still reports pre-existing formatting differences outside these five rustfmt-skipped snapshot lines; the scoped harness formatting and final targeted Clippy pass are recorded separately.',
    'Grouped sunset rows, other grouped references and source46 combined color-cycle validation remain outside this migration.',
    'No graph numerator expansion, broad snapshot rerecording, production fix, main integration or generation gate approval is implied.',
  ],
};
assert.equal(receipt.sourceAndFeatureIdentity.approvedTestSourceSha256, harness.oldOwnerSha256);
assert.equal(receipt.sourceAndFeatureIdentity.migratedTestSource.sha256, harness.newOwnerSha256);
fs.writeFileSync(`${evidence}/receipt.json`, `${JSON.stringify(receipt, null, 2)}\n`);
