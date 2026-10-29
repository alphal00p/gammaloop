import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';

const workspace = process.cwd();
const oldWorkspace = '/shared/lcnbr/dev/gammaloop/generation-pocket-weight-historical-source43';
const current = path.join(workspace, 'target/ghost-weight-certificate/current/states/ghost_weight_certificate');
const historical = path.join(oldWorkspace, 'target/ghost-weight-certificate/historical/states/ghost_weight_certificate');
const evidence = 'docs/architecture/evidence/generation-integration/wp3/ghost-weight-certificate';
const sha = (file) => crypto.createHash('sha256').update(fs.readFileSync(file)).digest('hex');
const file = (location) => ({path: location, sha256: sha(location)});
const same = (left, right, message) => assert.deepEqual(left, right, message);
const rows = JSON.parse(fs.readFileSync(path.join(evidence, 'source44-rows.json'))).cases;
const physical = JSON.parse(fs.readFileSync(path.join(evidence, 'physical-summary.json')));
const specs = [
  {id: 'epem_bbh_qcd2', family: 'epem_bbh', label: 'bbh_qcd2', kind: 'cross_sections', pristine: 'epem'},
  {id: 'slow_aa_ttx_qcd2', family: 'slow_aa_ttx', label: 'qcd2', kind: 'cross_sections', pristine: 'aa'},
  {id: 'vacuum_nonfactorizable_l2', family: 'vacuum_nonfactorizable', label: 'l2', kind: 'amplitudes'},
  {id: 'vacuum_nonfactorizable_l3', family: 'vacuum_nonfactorizable', label: 'l3', kind: 'amplitudes'},
  {id: 'vacuum_nonfactorizable_l4', family: 'vacuum_nonfactorizable', label: 'l4', kind: 'amplitudes'},
];
const families = {};
const inputCommands = (log) => {
  const groups = [];
  for (const line of fs.readFileSync(log, 'utf8').split('\n')) {
    const marker = 'Running command ';
    const offset = line.indexOf(marker);
    if (offset < 0) continue;
    const command = line.slice(offset + marker.length);
    if (command === 'import model sm.json') groups.push([command]);
    else if (command.startsWith('generate ') || command.startsWith('set global kv global.n_cores.feyngen=')) {
      assert(groups.length, 'generation before import');
      groups.at(-1).push(command);
    }
  }
  same(groups.length, 3, 'three independent CLI model imports');
  return groups;
};
const currentCommands = inputCommands(path.join(workspace, 'target/ghost-weight-certificate/current-capture.log'));
const historicalCommands = inputCommands(path.join(oldWorkspace, 'target/ghost-weight-certificate/historical-capture.log'));
same(historicalCommands, currentCommands, 'historical and current effective CLI commands differ');
for (const family of ['epem_bbh', 'slow_aa_ttx', 'vacuum_nonfactorizable']) {
  const stages = family === 'epem_bbh' ? ['dd_qcd0', 'bbh_qcd0', 'bbh_qcd1', 'bbh_qcd2', 'bbh_qcd2-no-zero-filter']
    : family === 'slow_aa_ttx' ? ['qcd0', 'qcd1', 'qcd2', 'qcd2-no-zero-filter']
      : ['l2', 'l2-no-zero-filter', 'l3', 'l3-no-zero-filter', 'l4', 'l4-no-zero-filter'];
  const routes = {};
  for (const [route, root] of [['current', current], ['source43', historical]]) {
    const folder = path.join(root, family);
    const live = {};
    for (const stage of ['post-import', 'pre-generation', ...stages.map((name) => `${name}-after-generation`)]) {
      const model = path.join(folder, `${stage}-live-model.json`);
      const parameters = path.join(folder, `${stage}-live-parameters.json`);
      const settings = path.join(folder, `${stage}-live-settings.json`);
      const settingsObject = JSON.parse(fs.readFileSync(settings));
      same(settingsObject.restriction, 'default', `${family}/${route}/${stage} restriction`);
      live[stage] = {model: file(model), resolvedParameters: file(parameters), settings: file(settings),
        settingsValue: settingsObject};
    }
    const modelHash = live['post-import'].model.sha256;
    const parametersHash = live['post-import'].resolvedParameters.sha256;
    for (const stage of Object.values(live)) {
      same(stage.model.sha256, modelHash, `${family}/${route} model changed after import`);
      same(stage.resolvedParameters.sha256, parametersHash, `${family}/${route} card changed after import`);
    }
    routes[route] = {live, modelSha256: modelHash, resolvedParametersSha256: parametersHash};
  }
  same(routes.source43.modelSha256, routes.current.modelSha256, `${family} live models differ`);
  same(routes.source43.resolvedParametersSha256, routes.current.resolvedParametersSha256, `${family} cards differ`);
  const normalizeSettings = (value) => {
    const copy = structuredClone(value);
    delete copy.cli_settings.state?.folder;
    return copy;
  };
  for (const stage of ['post-import', 'pre-generation', ...stages.map((name) => `${name}-after-generation`)]) {
    same(normalizeSettings(routes.source43.live[stage].settingsValue),
      normalizeSettings(routes.current.live[stage].settingsValue), `${family}/${stage} resolved settings differ`);
  }
  families[family] = {stages, routes, postImportAndPostGenerationModelCardStable: true,
    normalizedSettingsEqualAtEveryStage: true,
    orderedEffectiveInputCommands: currentCommands[['epem_bbh', 'slow_aa_ttx', 'vacuum_nonfactorizable'].indexOf(family)],
    source43AndCurrentEffectiveCommandsByteEqual: true,
    coreSetting: family === 'vacuum_nonfactorizable' ? {originalSlowTest: 16, diagnostic: 2} : {default: 1},
    lifecycle: 'new CLI; import model sm.json; each feyngen_str command appends --clear-existing-processes --only-diagrams; each collection has one process; save dot; capture live state after generation'};
}
const cases = {};
for (const spec of specs) {
  const sourceRow = rows.find((row) => row.id === spec.id);
  const labels = [spec.label, `${spec.label}-no-zero-filter`];
  const captures = {};
  for (const label of labels) {
    const key = `${spec.id === 'epem_bbh_qcd2' ? 'epem_bbh_qcd2' : spec.id === 'slow_aa_ttx_qcd2' ? 'slow_aa_ttx_qcd2' : spec.id.replace('nonfactorizable_', '')}/${label}`;
    const row = physical.cases[key];
    assert(row, `physical case absent ${key}`);
    const paths = Object.fromEntries([['source43', historical], ['current', current]].map(([route, root]) => {
      const folder = path.join(root, spec.family);
      const definition = path.join(folder, `${label}-process-definition.json`);
      const witness = path.join(folder, `${label}-graph-witness.tsv`);
      return [route, {definition: file(definition), graphWitness: file(witness)}];
    }));
    same(paths.source43.definition.sha256, paths.current.definition.sha256, `${key} definition differs`);
    same(paths.source43.graphWitness.sha256, paths.current.graphWitness.sha256, `${key} raw graph witness differs`);
    captures[label] = {...row, paths};
  }
  cases[spec.id] = {qualifiedTest: sourceRow.qualified_test,
    source43Expected: sourceRow.source43.expected, source44Expected: sourceRow.source44.expected,
    currentUnchangedExpected: sourceRow.current.expected,
    command: `generate ${spec.kind === 'cross_sections' ? 'xs' : 'amp'} ${sourceRow.exact_process_fragment} --clear-existing-processes --only-diagrams`,
    actualFiltered: captures[spec.label].sum, ghostSignOmitted: captures[spec.label].ghostlessSum,
    captures,
    pristine: spec.pristine ? {
      current: file(path.join(workspace, `target/ghost-weight-certificate/controls/current-${spec.pristine === 'epem' ? 'epem' : 'aa'}.log`)),
      source43: file(path.join(oldWorkspace, `target/ghost-weight-certificate/controls/historical-${spec.pristine === 'epem' ? 'epem' : 'aa'}.log`)),
      result: 'old assertion reached and failed against exact generated result; preceding cells passed',
    } : null};
}
const approvedRoot = path.join(workspace, 'tests/artifacts/pocket_weight_diagnosis_round2');
const slowRoot = path.join(current, 'slow_aa_ttx');
const reuse = {approvedCase: 'failing::test_generate_aa_ttx_cross_section_slow_filter',
  extendedCase: 'slow::test_generate_aa_ttx_cross_section',
  priorScopeAReview: file('target/pocket-weight-diagnosis/review-astra-scope-a-round2-JOlXh9/review.json'),
  priorScopeBReview: file('target/pocket-weight-diagnosis/review-astra-scope-b-YmoQie/review.json'),
  approvedMigrationReview: file('target/pocket-weight-diagnosis/review-astra-migration-YGwHhO/review.json'),
  populations: {}};
for (const [approvedLabel, slowLabel] of [['filtered', 'qcd2'], ['unfiltered', 'qcd2-no-zero-filter']]) {
  const files = {};
  for (const [kind, suffix] of [['witness', 'graph-witness.tsv'], ['definition', 'process-definition.json']]) {
    const earlier = file(path.join(approvedRoot, `${approvedLabel}-${suffix}`));
    const duplicate = file(path.join(slowRoot, `${slowLabel}-${suffix}`));
    same(earlier.sha256, duplicate.sha256, `slow photon and approved failing photon ${approvedLabel}/${kind} differ`);
    files[kind] = {earlier, duplicate, byteEqual: true};
  }
  const folder = 'processes/cross_sections/aa_ttx/default';
  const earlierDots = path.join(approvedRoot, `${approvedLabel}-dot`, folder);
  const duplicateDots = path.join(slowRoot, `${slowLabel}-dot`, folder);
  const dotNames = fs.readdirSync(earlierDots).filter((name) => name.endsWith('.dot')).sort();
  same(dotNames, fs.readdirSync(duplicateDots).filter((name) => name.endsWith('.dot')).sort(), `${approvedLabel} DOT membership`);
  const makeTreeDigest = (root) => crypto.createHash('sha256')
    .update(dotNames.map((name) => `${name}\0${sha(path.join(root, name))}\n`).join('')).digest('hex');
  const earlierTree = makeTreeDigest(earlierDots);
  same(earlierTree, makeTreeDigest(duplicateDots), `${approvedLabel} DOT bytes differ`);
  reuse.populations[approvedLabel] = {files, dotCount: dotNames.length,
    dotTreeSha256: earlierTree, earlierDotDirectory: earlierDots, duplicateDotDirectory: duplicateDots};
}
for (const [earlierName, duplicateName] of [['post-import-live-model.json', 'post-import-live-model.json'],
  ['post-import-live-parameters.json', 'post-import-live-parameters.json']]) {
  same(sha(path.join(approvedRoot, earlierName)), sha(path.join(slowRoot, duplicateName)),
    'slow photon live model/card differs from approved failing case');
}
const normalize = (settingsFile) => {
  const settings = JSON.parse(fs.readFileSync(settingsFile));
  delete settings.cli_settings.state.folder;
  return settings;
};
same(normalize(path.join(approvedRoot, 'post-import-live-settings.json')),
  normalize(path.join(slowRoot, 'post-import-live-settings.json')),
  'slow photon normalized live settings differ from approved failing case');
reuse.postImportModelCardAndNormalizedSettingsEqual = true;
const receipt = {
  schema: 1,
  status: 'review candidate; no expectation or production change in this certificate',
  source44Mapping: file(path.join(evidence, 'source44-rows.json')),
  source: {
    frozenParent: '59b57cead76e576fa9f23e3996f5614de939fe35',
    approvedMigration: '309631912ff5693ab38df4c748e5edfa3f94cd00',
    source43: '74fcd8e57e9fe256d4bc34b5d8ad93c97fc28848',
    source44: '087bda82ce5a14dab092776ea470e66278dda12d',
    source43Baseline: '81f7e3c0dcf9909472d22cc70a2bb2a79cef4c7c',
    source43DiagnosticAdaptation: 'c2127c1c3342a365c5f64c6431cfb5259bb9beac',
    historicalDiagnosticBookmark: 'generation-ghost-weight-historical-diagnostic',
    historicalPristineBookmark: 'generation-ghost-weight-historical-source43-pristine',
  },
  producer: {
    current: {
      binary: file(path.join(workspace, 'target/ghost-weight-certificate/current-test_feyngen')),
      diagnosticSource: file(path.join(workspace, 'tests/tests/test_feyngen.rs')),
      cargoLock: file(path.join(workspace, 'Cargo.lock')),
      flakeLock: file(path.join(workspace, 'flake.lock')),
      kernel: '1a1829b1f3d0a88b53b44085b3b9f34ef77945b3',
      captureLog: file(path.join(workspace, 'target/ghost-weight-certificate/current-capture.log')),
      pristineBinary: file(path.join(workspace, 'target/pocket-weight-diagnosis/migration/named-new-pass-test_feyngen')),
    },
    source43: {
      binary: file(path.join(oldWorkspace, 'target/ghost-weight-certificate/historical-test_feyngen')),
      diagnosticSource: file(path.join(oldWorkspace, 'target/ghost-weight-certificate/historical-diagnostic-source.rs')),
      cargoLock: file(path.join(oldWorkspace, 'Cargo.lock')),
      flakeLock: file(path.join(oldWorkspace, 'flake.lock')),
      kernel: 'c78cffc20f318df1fce50539a587cbe8ffc717c0',
      captureLog: file(path.join(oldWorkspace, 'target/ghost-weight-certificate/historical-capture.log')),
      pristineBinary: file(path.join(oldWorkspace, 'target/pocket-weight-diagnosis/binaries/source43-pristine-test_feyngen')),
    },
    features: 'Cargo defaults; ci-optim profile; no feature override; Nix develop with recorded locks; two build jobs',
    sourceModel: {
      current: file(path.join(workspace, 'assets/models/json/sm/sm.json')),
      source43: file(path.join(oldWorkspace, 'assets/models/json/sm/sm.json')),
    },
  },
  families, cases, approvedFailingPhotonCaseReuse: reuse,
  physicalCertificate: {
    analyzer: file(path.join(evidence, 'analyze-graphs.mjs')),
    compact: file(path.join(evidence, 'physical-summary.json')),
    full: file(path.join(workspace, 'target/ghost-weight-certificate/full-physical-joins.json')),
    method: physical.method,
    semanticMutationControls: physical.semanticMutationControls,
  },
  verdict: 'All five source43/current per-graph physical memberships, mapped cuts, factors, and weights agree. Independently derived directed AutG and closed Grassmann signs give all five current sums. Omitting only ghost-loop parity gives all five old literals. The first observed unequal boundary is each old snapshot against its own producer; no first buggy production commit is established. Independent Astra certificate review required before migration.',
  limitations: [
    'The pristine epep and slow-aa tests stop at their stale named rows; later assertions are unclaimed.',
    'The vacuum capture uses two generator cores instead of the test’s 16 and skips preceding unrelated grouped/full-mode assertions. Exact pristine vacuum16 execution and whole-test pass are unclaimed.',
    'Zero-color-off membership checks do not independently prove every removed graph has symbolic zero color.',
    'No complete graph numerator was expanded, no grouped row or source46 color-cycle combined validation was attempted.',
  ],
};
same(receipt.producer.sourceModel.source43.sha256, receipt.producer.sourceModel.current.sha256, 'source SM model bytes differ');
fs.writeFileSync(path.join(evidence, 'receipt.json'), `${JSON.stringify(receipt, null, 2)}\n`);
