import fs from 'node:fs';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';

const base = 'docs/architecture/evidence/generation-integration/wp3/ghost-weight-certificate';
const aRoot = 'target/pocket-weight-diagnosis/review-astra-ghost-scope-a-L13RwH';
const bRoot = 'target/ghost-weight-certificate/review-astra-scope-b-8wU1pd';
const target = 'target/ghost-weight-certificate/round2';
const file = (location) => ({path: location, sha256: crypto.createHash('sha256')
  .update(fs.readFileSync(location)).digest('hex')});
const a = JSON.parse(fs.readFileSync(`${aRoot}/review.json`));
const b = JSON.parse(fs.readFileSync(`${bRoot}/review.json`));
const old = JSON.parse(fs.readFileSync(`${base}/receipt.json`));
const newSummary = JSON.parse(fs.readFileSync(`${base}/round2/physical-summary.json`));
const guard = JSON.parse(fs.readFileSync(`${target}/raw-guard.json`));
const controls = JSON.parse(fs.readFileSync(`${target}/controls.json`));
const accessibility = JSON.parse(fs.readFileSync(`${target}/semantic-table.json`));
assert.equal(a.verdict, 'REQUEST_CHANGES');
assert.equal(b.verdict, 'REQUEST_CHANGES');
assert.equal(guard.originalBoundFiles, 199);
assert.equal(guard.graphRecords, 8318);
assert.equal(controls.results.length, 33);
assert(controls.results.every((entry) => entry.expected === entry.outcome));
assert.deepEqual(accessibility, {scoped: true, errors: [], table_count: 1, th_count: 6});
const certificate = {
  schema: 2,
  status: 'corrected five-row certificate; new independent Astra A/B correction reviews required',
  change: 'vswwmovz',
  parent: 'decb85d4a15c66b6c69675cfd641465dfbb70757',
  round1Bookmark: 'generation-ghost-weight-certificate-round1',
  candidateBookmark: 'generation-ghost-weight-certificate',
  approvedMigrationBookmark: 'generation-pocket-weight-diagnosis-approved-migration',
  approvedMigrationCommit: '309631912ff5693ab38df4c748e5edfa3f94cd00',
  source43: old.source.source43,
  source44: old.source.source44,
  reviewInputs: {
    scopeA: file(`${aRoot}/review.json`),
    scopeB: file(`${bRoot}/review.json`),
    scopeBOldControls: file(`${bRoot}/controls.json`),
    scopeBOldFactorRepro: file(`${bRoot}/factor-repro.json`),
    scopeBSupplement: file(`${bRoot}/supplement.json`),
    scopeBIndependentAnalyzer: file(`${bRoot}/independent.mjs`),
    scopeAOldAccessibilityFailure: file(`${aRoot}/semantic-table-reproducer.log`),
    scopeADocumentCommands: file(`${aRoot}/docs-commands.json`),
  },
  round1: {
    receipt: file(`${base}/receipt.json`),
    summary: file(`${base}/physical-summary.json`),
    fullPhysicalJoin: file('target/ghost-weight-certificate/full-physical-joins.json'),
    sourceAnalyzerSha256AtRound1: old.physicalCertificate.analyzer.sha256,
    preEditReceiptValidation: file(`${target}/round1-receipt-validation.log`),
  },
  corrected: {
    analyzer: file(`${base}/analyze-graphs.mjs`),
    source44Mapping: file(`${base}/source44-rows.json`),
    report: file('docs/architecture/generation-ghost-weight-certificate.typ'),
    registry: file('docs/developers.toml'),
    summary: file(`${base}/round2/physical-summary.json`),
    fullPhysicalJoin: file(`${target}/full-physical-joins.json`),
    analyzerReplayLog: file(`${target}/analyzer-replay.log`),
    semanticControls: file(`${target}/controls.json`),
    rawGuard: file(`${target}/raw-guard.json`),
    scopedHtml: file(`${target}/certificate.html`),
    scopedPdf: file(`${target}/certificate.pdf`),
    scopedAccessibility: file(`${target}/semantic-table.json`),
    scopedDocCheckLog: file(`${target}/scoped-doc-check.log`),
    registryTitleCheckLog: file(`${target}/docs-builder-check.log`),
    pinnedJsSyntaxAndWhitespaceLog: file(`${target}/js-checks.log`),
    accessibilityRegexAttempt: file(`${target}/semantic-table-attempt1.json`),
  },
  dispositions: {
    B1: 'Complete bounded scalar grammar consumes every TSV and DOT factor token; exact BigInt rational evaluation, structural head occurrence and argument checks, both recorded scalar values and graph-derived factor checks; all six reviewer factor mutations and full mirrored readCase corruption rejected without hash guards; valid equivalent spellings accepted.',
    B2: 'Open Grassmann chains require one directed incoming and outgoing edge at internal vertices and model/role-correct source/sink external endpoints; reversed incoming positron arrow rejected locally; valid two-edge ghost-circuit reversal accepted.',
    B3: 'Both cut sides and crossing arrays have unique known half-edges; physical and selected-subgraph thresholds require exact paired crossing incidence; proper subgraphs remain valid, missing mate rejected locally.',
    B4: 'Vacuum prose now states ExternalFermionOrderingSign(1) is retained, no external legs and no AntiFermionSpinSumSign head.',
    A1: 'Six unchanged table label texts wrapped in table.header; scoped Typst HTML contains one table and six th cells; DocumentAudit.errors empty versus retained old failure.',
  },
  preservedInputs: {
    rawGuard: guard,
    currentDiagnosticSource: old.producer.current.diagnosticSource,
    source43DiagnosticSource: old.producer.source43.diagnosticSource,
    currentBinary: old.producer.current.binary,
    source43Binary: old.producer.source43.binary,
    currentKernel: old.producer.current.kernel,
    source43Kernel: old.producer.source43.kernel,
    modelSha256: old.families.epem_bbh.routes.current.modelSha256,
    resolvedParametersSha256: old.families.epem_bbh.routes.current.resolvedParametersSha256,
    sourceModelSha256: old.producer.sourceModel.current.sha256,
  },
  sourceBoundSixteenWorkerControls: {
    attribution: 'reused independent Astra scope A review evidence, not rerun by this correction',
    current: a.vacuumAdjudication.current16.binary,
    historical: a.vacuumAdjudication.historical16.binary,
    currentScratchSource: a.vacuumAdjudication.current16.scratchSource,
    historicalScratchSource: a.vacuumAdjudication.historical16.scratchSource,
    commands: a.vacuumAdjudication.commands,
    result: a.vacuumAdjudication.result,
  },
  sourceMatchedRustChecks: {
    attribution: 'reused independent Astra scope A checks because test_feyngen.rs SHA-256 is unchanged',
    helperSourceSha256: old.producer.current.diagnosticSource.sha256,
    rustfmtReadOnly: file(`${aRoot}/rustfmt-readonly.log`),
    targetedClippy: file(`${aRoot}/clippy.log`),
    result: a.documentationAndRust.rustfmtInterpretation,
    clippyResult: a.documentationAndRust.clippy,
  },
  checks: {
    analyzer: 'complete historical/current filtered and zero-color-off replay passed all ten datasets',
    originalControls: 7,
    newControls: controls.results.length,
    newRejectingControls: controls.results.filter((entry) => entry.outcome === 'rejected').length,
    newAcceptingControls: controls.results.filter((entry) => entry.outcome === 'accepted').length,
    oldCheckerFactorEscapes: 'scope B controls.json and factor-repro.json retain accepted invalid factors and mirrored full readCase corruption',
    oldAccessibilityFailure: 'scope A semantic-table-reproducer.log reports one table and zero th cells',
    newAccessibility: accessibility,
  },
  verdict: 'Unchanged captured graph mathematics is now checked by complete scalar parsing and local Grassmann/cut invariants. Five proposed sums and ghost-only omission counterfactuals remain exact. The first observed unequal boundary remains historical snapshot versus its own producer, without a claimed first buggy production commit.',
  limitations: [
    'No expectation migration or production physics edit is authorized by this correction.',
    'The scoped 16-worker harnesses do not execute the pristine vacuum test’s grouped prefix or establish a whole-test pass.',
    'Zero-color-off membership is bound; a symbolic proof of each removed color-zero graph is not claimed.',
    'No graph numerator expansion, grouped source44 row, source46 color-cycle combined validation, main integration or gate approval is claimed.',
  ],
};
assert.equal(certificate.corrected.fullPhysicalJoin.sha256, newSummary.fullWitnessSha256);
fs.writeFileSync(`${base}/round2/correction.json`, `${JSON.stringify(certificate, null, 2)}\n`);
