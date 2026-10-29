import fs from 'node:fs';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';

const base = 'docs/architecture/evidence/generation-integration/wp3/ghost-weight-certificate';
const target = 'target/ghost-weight-certificate/round3';
const aRoot = 'target/ghost-weight-certificate/review-astra-round2-scope-a-uZNf9j';
const bRoot = 'target/ghost-weight-certificate/review-astra-round2-scope-b-EuYuaM';
const file = (location) => ({path: location, sha256: crypto.createHash('sha256')
  .update(fs.readFileSync(location)).digest('hex')});
const round1 = JSON.parse(fs.readFileSync(`${base}/receipt.json`));
const round2 = JSON.parse(fs.readFileSync(`${base}/round2/correction.json`));
const a = JSON.parse(fs.readFileSync(`${aRoot}/review.json`));
const b = JSON.parse(fs.readFileSync(`${bRoot}/review.json`));
const repro = JSON.parse(fs.readFileSync(`${bRoot}/finding-repro.json`));
const summary = JSON.parse(fs.readFileSync(`${base}/round3/physical-summary.json`));
const guard = JSON.parse(fs.readFileSync(`${target}/raw-guard.json`));
const differential = JSON.parse(fs.readFileSync(`${target}/parser-differential.json`));
const controls = JSON.parse(fs.readFileSync(`${target}/controls.json`));
const priorControls = JSON.parse(fs.readFileSync(`${target}/prior-controls-replay.json`));
const accessibility = JSON.parse(fs.readFileSync(`${target}/semantic-table.json`));
assert.equal(a.verdict, 'APPROVE_CERTIFICATE');
assert.equal(b.verdict, 'REQUEST_CHANGES');
assert.equal(repro.mirroredJoin.accepted, true);
assert.equal(repro.duplicateInitial.localChecks, 'accepted');
assert.equal(guard.originalBindings, 199);
assert.equal(guard.round2Bindings, 40);
assert.deepEqual([differential.admittedCount, differential.explicitlyRejectedCount], [44, 24]);
assert.equal(controls.results.length, 30);
assert(controls.results.every((entry) => entry.expected === entry.outcome));
assert.deepEqual([priorControls.originalSeven, priorControls.round2Controls,
  priorControls.reviewerPriorControls, priorControls.reviewerExtraControls,
  priorControls.generatedScalarExpressions], [7, 33, 31, 111, 2000]);
assert.deepEqual(accessibility, {scoped: true, errors: [], table_count: 1,
  th_count: 6, table_markup_byte_equal_to_round2: true});
const receipt = {
  schema: 3,
  status: 'corrected five-row scalar checker; NEW independent Astra correction review required',
  change: 'tqlnmrvv',
  parent: '91b04d891860087f2bf9c363d12683f1c8dcabd3',
  author: 'Lucien Huber <im@lcnbr.ch>',
  bookmarks: {
    candidate: 'generation-ghost-weight-certificate',
    round2: 'generation-ghost-weight-certificate-round2',
    round1: 'generation-ghost-weight-certificate-round1',
    approvedMigration: 'generation-pocket-weight-diagnosis-approved-migration',
    approvedCertificate: 'generation-pocket-weight-diagnosis-approved-certificate',
  },
  scope: 'five source44 ungrouped rows only; no expectation or production edit',
  reviews: {
    scopeAApprovedRound2: file(`${aRoot}/review.json`),
    scopeBChangesRound2: file(`${bRoot}/review.json`),
    oldPrecedenceAndInitialReproducer: file(`${bRoot}/finding-repro.json`),
    oldReproducerSource: file(`${bRoot}/finding-repro.mjs`),
    oldExtraControlResults: file(`${bRoot}/extra-controls.json`),
    pinnedParserTokenLog: file(`${bRoot}/precedence.log`),
  },
  frozen: {
    round1Receipt: file(`${base}/receipt.json`),
    round1Summary: file(`${base}/physical-summary.json`),
    round1FullJoin: file('target/ghost-weight-certificate/full-physical-joins.json'),
    round2Receipt: file(`${base}/round2/correction.json`),
    round2Summary: file(`${base}/round2/physical-summary.json`),
    round2FullJoin: file('target/ghost-weight-certificate/round2/full-physical-joins.json'),
  },
  corrected: {
    analyzer: file(`${base}/analyze-graphs.mjs`),
    parserDifferentialSource: file(`${base}/round3/parser-differential.mjs`),
    receiptGeneratorSource: file(`${base}/round3/make-receipt.mjs`),
    receiptValidatorSource: file(`${base}/round3/validate-receipt.mjs`),
    parserDifferential: file(`${target}/parser-differential.json`),
    pinnedScalarOracleSource: file(`${target}/scalar-oracle.rs`),
    pinnedScalarOracleBinary: file(`${target}/scalar-oracle`),
    pinnedScalarOracleBuildLog: file(`${target}/scalar-oracle-initial-pass.log`),
    pinnedScalarOracleProbeLog: file(`${target}/scalar-oracle-probes.log`),
    priorControlsReplaySource: file(`${base}/round3/replay-prior-controls.mjs`),
    priorControlsReplay: file(`${target}/prior-controls-replay.json`),
    newControlsSource: file(`${base}/round3/controls.mjs`),
    newControls: file(`${target}/controls.json`),
    rawGuardSource: file(`${base}/round3/validate-guard.mjs`),
    rawGuard: file(`${target}/raw-guard.json`),
    fullJoin: file(`${target}/full-physical-joins.json`),
    summary: file(`${base}/round3/physical-summary.json`),
    analyzerLog: file(`${target}/analyzer-replay.log`),
    report: file('docs/architecture/generation-ghost-weight-certificate.typ'),
    registry: file('docs/developers.toml'),
    html: file(`${target}/certificate.html`),
    pdf: file(`${target}/certificate.pdf`),
    accessibility: file(`${target}/semantic-table.json`),
    docsCheckLog: file(`${target}/scoped-doc-check.log`),
    registryCheckLog: file(`${target}/docs-builder-check.log`),
    pinnedSourceChecksLog: file(`${target}/source-checks.log`),
  },
  reusedProducersAndInputs: {
    attribution: 'immutable round1/round2 receipts and independent Astra A qualifications; not rebuilt or rerun for this checker correction',
    currentDiagnosticSource: round1.producer.current.diagnosticSource,
    historicalDiagnosticSource: round1.producer.source43.diagnosticSource,
    currentElf: round1.producer.current.binary,
    historicalElf: round1.producer.source43.binary,
    currentKernel: round1.producer.current.kernel,
    historicalKernel: round1.producer.source43.kernel,
    liveModelSha256: round1.families.epem_bbh.routes.current.modelSha256,
    resolvedParametersSha256: round1.families.epem_bbh.routes.current.resolvedParametersSha256,
    scopedSixteenWorkerCurrent: round2.sourceBoundSixteenWorkerControls.current,
    scopedSixteenWorkerHistorical: round2.sourceBoundSixteenWorkerControls.historical,
  },
  findings: {
    'B1-R2': 'Unary minus now applies after power; -2^0 and (-2)^0 are distinct. Multiplication/division, signed bounded exponents, exact tokens, parentheses and full consumption are checked. Additive/chained/undefined/out-of-scope forms reject. The real GL1 mutation fails locally, on TSV-only/DOT-only/simultaneous reads for both routes, and in the mirrored route setup without hash guards.',
    'B5-R2': 'Initial left/right arrays require unique valid physical half-edges and exact sewn incidence on each side. Duplicate, missing, wrong and out-of-domain entries fail locally; reordering and prior relabel controls pass.',
    closedPrior: 'B1-R1 complete factor parsing, B2 open directed Grassmann flow, B3 selected-subgraph threshold pairing, B4 vacuum sign-head prose and A1 table headers remain closed; table markup is byte identical to approved round2.',
  },
  oldCheckerFailNewCheckerPass: {
    oldChecker: {tsvOnly: 'accepted', dotOnly: 'accepted', both: 'accepted',
      mirroredHistoricalCurrent: repro.mirroredJoin, duplicateInitialLocal: repro.duplicateInitial.localChecks},
    newChecker: {controls: controls.results, parserDifferential: {
      admitted: differential.admittedCount, rejected: differential.explicitlyRejectedCount}},
  },
  checks: {
    completeReplay: 'ten historical/current filtered and zero-color-off datasets pass; same five exact sums and ghost-only omission totals',
    fullJoinByteEqualToRound2: file(`${target}/full-physical-joins.json`).sha256
      === file('target/ghost-weight-certificate/round2/full-physical-joins.json').sha256,
    frozenBindings: [guard.originalBindings, guard.round2Bindings],
    originalControls: priorControls.originalSeven,
    round2Controls: priorControls.round2Controls,
    reviewerPriorControls: priorControls.reviewerPriorControls,
    reviewerExtraControls: priorControls.reviewerExtraControls,
    generatedScalarExpressions: priorControls.generatedScalarExpressions,
    newControls: controls.results.length,
    newRejected: controls.results.filter((entry) => entry.outcome === 'rejected').length,
    newAccepted: controls.results.filter((entry) => entry.outcome === 'accepted').length,
    differential: [differential.admittedCount, differential.explicitlyRejectedCount],
    scopedAccessibility: accessibility,
  },
  limitations: [
    'The scalar grammar is deliberately bounded and its 44-case differential is a scope certificate, not a general Symbolica parser equivalence theorem.',
    'The full pristine vacuum test and grouped prefix remain unclaimed; scoped 16-worker command evidence is reused from independent A.',
    'Zero-color-off membership is bound, without symbolic proof of every removed color-zero graph.',
    'No numerator expansion, grouped source44 row, source46 combined validation, expectation migration, production fix, main integration or gate release is implied.',
  ],
};
assert.equal(receipt.corrected.fullJoin.sha256, summary.fullWitnessSha256);
assert.equal(receipt.corrected.fullJoin.sha256, receipt.frozen.round2FullJoin.sha256);
fs.writeFileSync(`${base}/round3/correction.json`, `${JSON.stringify(receipt, null, 2)}\n`);
