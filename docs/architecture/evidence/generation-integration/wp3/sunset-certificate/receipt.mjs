import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import {fileURLToPath} from 'node:url';

const [historical, current, proofRoot, output] = process.argv.slice(2);
if (!historical || !current || !proofRoot || !output) {
  throw new Error('usage: node receipt.mjs <historical-capture> <current-capture> <proof-dir> <output.json>');
}
const scripts = path.dirname(fileURLToPath(import.meta.url));
const workspace = path.resolve(current, '../../..');
const historicalWorkspace = path.resolve(historical, '../../..');
const sha256 = (file) => crypto.createHash('sha256').update(fs.readFileSync(file)).digest('hex');
const artifact = (file) => ({path: file, bytes: fs.statSync(file).size, sha256: sha256(file)});
const allFiles = (directory) => fs.readdirSync(directory, {recursive: true, withFileTypes: true})
  .filter((entry) => entry.isFile())
  .map((entry) => path.join(entry.parentPath, entry.name)).sort();
const capture = (root) => allFiles(root).filter((file) =>
  !file.includes(`${path.sep}state${path.sep}`) &&
  !file.includes(`${path.sep}drawings${path.sep}`))
  .map(artifact);
const joined = JSON.parse(fs.readFileSync(path.join(proofRoot, 'join.json')));
const input = JSON.parse(fs.readFileSync(path.join(proofRoot, 'inputs.json')));
const ensure = (file, expected) => {
  if (sha256(file) !== expected) throw new Error(`identity mismatch: ${file}`);
  return artifact(file);
};
const source43Pristine = path.join(historicalWorkspace,
  'target/sunset-certificate/binaries/source43-pristine-test_feyngen');
const currentBinary = path.join(workspace,
  'target/sunset-certificate/binaries/current-diagnostic-test_feyngen');
const source43Binary = path.join(historicalWorkspace,
  'target/sunset-certificate/binaries/source43-diagnostic-test_feyngen');
const roots = {source43: historical, current};
const result = {
  schema: 1,
  status: 'bounded preparation certificate; new independent Astra proof review required',
  physicalParent: '59b57cead76e576fa9f23e3996f5614de939fe35',
  source43: '74fcd8e57e9fe256d4bc34b5d8ad93c97fc28848',
  source44: '087bda82ce5a14dab092776ea470e66278dda12d',
  workspace, historicalWorkspace,
  changeId: 'plkmzuwm', bookmark: 'generation-sunset-certificate',
  historicalBookmark: 'generation-sunset-historical-source43',
  sourceInputs: {
    currentModel: ensure(path.join(workspace, 'assets/models/json/sm/sm.json'),
      'a0193ae09a6f920c8fcfbe8231c6a70e820b82b573257e07c216bdcd8e57acd9'),
    historicalModel: ensure(path.join(historicalWorkspace, 'assets/models/json/sm/sm.json'),
      'a0193ae09a6f920c8fcfbe8231c6a70e820b82b573257e07c216bdcd8e57acd9'),
    currentRestriction: ensure(path.join(workspace, 'assets/models/json/sm/restrict_default.json'),
      'c683f9abb9459375865ef694edea38890ffd78a0b8f118a25db3c3de09a8a1cb'),
    historicalRestriction: ensure(path.join(historicalWorkspace,
      'assets/models/json/sm/restrict_default.json'),
    'c683f9abb9459375865ef694edea38890ffd78a0b8f118a25db3c3de09a8a1cb'),
    currentLock: artifact(path.join(workspace, 'Cargo.lock')),
    historicalLock: artifact(path.join(historicalWorkspace, 'Cargo.lock')),
    currentFlakeLock: artifact(path.join(workspace, 'flake.lock')),
    historicalFlakeLock: artifact(path.join(historicalWorkspace, 'flake.lock')),
  },
  documentation: {
    note: artifact(path.join(workspace, 'docs/architecture/generation-sunset-certificate.typ')),
    registry: artifact(path.join(workspace, 'docs/developers.toml')),
  },
  kernels: {source43: 'c78cffc20f318df1fce50539a587cbe8ffc717c0',
    current: '1a1829b1f3d0a88b53b44085b3b9f34ef77945b3'},
  compile: {features: 'Cargo defaults, no explicit feature flags', profile: 'ci-optim',
    workers: 16, cargoBuildJobs: 2},
  binaries: {
    source43Pristine: ensure(source43Pristine,
      '5e09a8dc940f56bbde9c9484f70bf76537892ebc39e094491e2981c6285d7f69'),
    source43Diagnostic: ensure(source43Binary,
      'b58dea4f73cfd7c28a2f72b2e743f1f1b8ce1e8d10158a623b5a00783e8c0358'),
    currentDiagnostic: ensure(currentBinary,
      'ada5c1c77ce8c8fe42ab2c4aad4c982f7e3f55503a217992455d3d5a83d316dd'),
    source43FrameDiagnostic: ensure(path.join(historicalWorkspace,
      'target/sunset-certificate/binaries/source43-frame-diagnostic-test_feyngen'),
    '734dd778ff0334340a97b5e33dfe6aecacbb2fe2b56f1fd34346614bb2f7f8c2'),
    currentFrameDiagnostic: ensure(path.join(workspace,
      'target/sunset-certificate/binaries/current-frame-diagnostic-test_feyngen'),
    '9e873c0af39b155bad7b9661dfd950c2ee213c6fef7d70335a07be8cbc04243a'),
  },
  diagnosticSources: {
    source43: artifact(path.join(historicalWorkspace, 'target/sunset-certificate/source43-diagnostic-source.rs')),
    current: artifact(path.join(workspace, 'target/sunset-certificate/current-diagnostic-source.rs')),
  },
  frameObservationPatches: {
    source43: artifact(path.join(historicalWorkspace,
      'target/sunset-certificate/source43-frame-observation.patch')),
    current: artifact(path.join(workspace,
      'target/sunset-certificate/current-frame-observation.patch')),
  },
  commands: {
    build: 'CARGO_BUILD_JOBS=2 nix develop -c cargo test --locked -p gammaloop-integration-tests --test test_feyngen --profile ci-optim --no-run',
    historicalPristine: 'nix develop -c target/sunset-certificate/binaries/source43-pristine-test_feyngen slow::test_vacuum_amplitude_generation --exact --nocapture',
    diagnostic: 'nix develop -c target/sunset-certificate/binaries/<producer>-diagnostic-test_feyngen diagnose_two_loop_vacuum_sunsets --exact --ignored --nocapture',
    currentAssertion: 'nix develop -c target/sunset-certificate/binaries/current-diagnostic-test_feyngen slow::test_vacuum_amplitude_generation --exact --nocapture',
    frameDiagnostic: 'nix develop -c env SUNSET_FRAME_CAPTURE=<task-owned-frame-root> target/sunset-certificate/binaries/<producer>-frame-diagnostic-test_feyngen diagnose_two_loop_vacuum_sunsets --exact --ignored --nocapture',
    replayOrder: ['inputs.mjs', 'model.mjs', 'join.mjs', 'frame.mjs', 'mutations.mjs',
      'internal-frame.mjs', 'internal-mutations.mjs', 'positive.mjs', 'receipt.mjs',
      'verify-receipt.mjs'],
  },
  exits: {source43Pristine: 101, currentAssertion: 101, source43Diagnostic: 0,
    currentDiagnostic: 0, source43FrameDiagnostic: 0, currentFrameDiagnostic: 0,
    firstFrameBuildSource43: 101, firstFrameBuildCurrent: 101,
    firstFrameBuildReason: 'serde_json was a dev-dependency; corrected to existing FeynmanDiagram::to_json without manifest changes',
    positiveControls: 7,
    rejectedSemanticMutations: 22, rejectedInternalFrameMutations: 8},
  logs: {
    source43Pristine: artifact(path.join(historicalWorkspace, 'target/sunset-certificate/source43-pristine.log')),
    source43Build: artifact(path.join(historicalWorkspace, 'target/sunset-certificate/source43-build.log')),
    source43Capture: artifact(path.join(historicalWorkspace, 'target/sunset-certificate/source43-capture.log')),
    currentBuild: artifact(path.join(workspace, 'target/sunset-certificate/current-build.log')),
    currentCapture: artifact(path.join(workspace, 'target/sunset-certificate/current-capture.log')),
    currentAssertion: artifact(path.join(workspace, 'target/sunset-certificate/current-assertion-control.log')),
    source43FirstFrameBuildFailure: artifact(path.join(historicalWorkspace,
      'target/sunset-certificate/source43-frame-build.log')),
    currentFirstFrameBuildFailure: artifact(path.join(workspace,
      'target/sunset-certificate/current-frame-build.log')),
    source43FrameBuild: artifact(path.join(historicalWorkspace,
      'target/sunset-certificate/source43-frame-build-final.log')),
    currentFrameBuild: artifact(path.join(workspace,
      'target/sunset-certificate/current-frame-build-final.log')),
    source43FrameCapture: artifact(path.join(historicalWorkspace,
      'target/sunset-certificate/source43-frame-capture.log')),
    currentFrameCapture: artifact(path.join(workspace,
      'target/sunset-certificate/current-frame-capture.log')),
  },
  validation: {
    targetedClippy: {exit: 0,
      command: 'CARGO_BUILD_JOBS=2 nix develop -c cargo clippy --locked -p gammaloop-integration-tests --test test_feyngen --profile ci-optim -- -D warnings',
      log: artifact(path.join(workspace, 'target/sunset-certificate/clippy.log'))},
    finalFormattingTypstAndDocs: {exit: 0,
      command: "nix develop -c bash -c 'cargo fmt --all -- --check && typst compile docs/architecture/generation-sunset-certificate.typ target/sunset-certificate/proof/note.pdf && cargo run --locked -p alphal00p-docs-builder -- check'",
      log: artifact(path.join(workspace, 'target/sunset-certificate/final-static-doc-check.log'))},
    pythonRuff: {
      changedPythonFiles: 0,
      firstAttempt: {exit: 127,
        command: 'nix develop -c ruff check tests/resources/python/tests/unit/test_process_generation.py',
        log: artifact(path.join(workspace, 'target/sunset-certificate/ruff.log'))},
      pinnedAttempt: {exit: 1, version: '0.16.3',
        command: 'nix develop -c /nix/store/sddv779f9zwcm72i0nsa8sgy1zjcl2n4-ruff-0.16.3/bin/ruff check tests/resources/python/tests/unit/test_process_generation.py',
        log: artifact(path.join(workspace, 'target/sunset-certificate/ruff-pinned.log'))},
      inheritedFindings: '4 in unchanged historical Python test: I001 import order, F401 Path, F401 parse_python_expression, F841 exception binding',
    },
  },
  liveInputs: {modelSha256: input.modelSha256, parametersSha256: input.parametersSha256,
    normalizedSettingsSha256: input.normalizedSettingsSha256,
    processDefinitions: Object.fromEntries(['all-raw', 'all-grouped', 'nonfactorizable-raw',
      'nonfactorizable-grouped'].map((row) => [row, input.stages.current[row].definitionSha256]))},
  summaries: {source43: joined.native.source43, current: joined.native.current},
  proofArtifacts: Object.fromEntries(fs.readdirSync(proofRoot).filter((name) => name.endsWith('.json'))
    .map((name) => [name, artifact(path.join(proofRoot, name))])),
  proofScripts: Object.fromEntries(fs.readdirSync(scripts).filter((name) => name.endsWith('.mjs'))
    .map((name) => [name, artifact(path.join(scripts, name))])),
  captures: Object.fromEntries(Object.entries(roots).map(([route, root]) => [route, capture(root)])),
  internalFrames: {
    source43: allFiles(path.join(historicalWorkspace,
      'target/sunset-certificate/frames/source43')).map(artifact),
    current: allFiles(path.join(workspace, 'target/sunset-certificate/frames/current')).map(artifact),
  },
  limits: [
    'Only two grouped two-loop vacuum rows and paired raw controls; other source44 rows excluded.',
    'Original whole tests fail at the first old grouped label and never reach the later row.',
    'Source44 historical labels are provenance, not accepted expectations.',
    'Finite production 4D/SU(3) local-factor and tensor certificate; no whole-graph expansion.',
    'Preparation parent only; ordered-head validation and new Astra proof review remain.',
    'No production, assertion, lock, NixCI, upload, push, fetch or PR change.',
  ],
};
fs.writeFileSync(output, JSON.stringify(result, null, 2) + '\n');
console.log(`receipt binds ${result.captures.source43.length} source43 and ` +
  `${result.captures.current.length} current capture files`);
