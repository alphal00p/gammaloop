import fs from 'node:fs';
import {liveCards} from './live-rule.mjs';
import path from 'node:path';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';

const [historical, current, output] = process.argv.slice(2);
if (!historical || !current || !output) {
  throw new Error('usage: node inputs.mjs <historical-capture> <current-capture> <output.json>');
}
const hash = (bytes) => crypto.createHash('sha256').update(bytes).digest('hex');
const rows = ['all-raw', 'all-grouped', 'nonfactorizable-raw', 'nonfactorizable-grouped'];
const bytes = (root, name) => fs.readFileSync(path.join(root, name));
const model = bytes(historical, 'post-import-model.json');
const parameters = bytes(historical, 'post-import-parameters.json');
if (!model.equals(bytes(current, 'post-import-model.json'))) {
  throw new Error('live post-import model or input parameter card differs between producers');
}
const normalizeSettings = (buffer) => {
  const value = JSON.parse(buffer);
  delete value.cli_settings.state.folder;
  return value;
};
const settings = normalizeSettings(bytes(historical, 'post-import-settings.json'));
const expectedSettings = {cli_settings: {global: {n_cores: {feyngen: 16}}, state: {}},
  default_runtime_settings: {},
  model_fingerprint: '6e16ebaf9300927398e0e73855f43d72268a46d09c8337974efacf854ea11d9b',
  restriction: 'default'};
assert.deepStrictEqual(settings, expectedSettings, 'actual post-import runtime policy');
assert.deepStrictEqual(normalizeSettings(bytes(current, 'post-import-settings.json')),
  expectedSettings, 'current post-import runtime policy');
function expectedDefinition(row) {
  const filtered = row.endsWith('-raw');
  const nonfactorizable = row.startsWith('nonfactorizable');
  return {
    request: {
      process: {incoming: [], outgoing_alternatives: [[]]},
      kind: 'amplitude', loops: {total: {start: 2, end: 2}},
      restrictions: {
        particle_allow: [['d', 'g', 'ghG', 't', 'u'].map((value) => ({kind: 'name', value}))],
        particle_exclude: [], vertex_allow: [], vertex_exclude: [],
      },
      orders: {plain: {QED: {inclusive: {min: 0, max: 0}}}, squared: {}, perturbative: {}},
      topology: {
        bridges: nonfactorizable ? {AtMost: 0} : 'Any',
        self_loops: true, zero_flow_edges: false, max_vertices: null,
        self_energies: {massive: false, massless: false},
        tadpoles: {massive: false, massless: false},
        snails: {massive: false, massless: false},
        fermion_loops: 'unrestricted',
        factorized_loops: nonfactorizable ? {inclusive: {min: 1, max: 1}} : 'unrestricted',
      },
      cuts: {amplitude: {orders: {plain: {}, squared: {}, perturbative: {}},
        loops: 'unrestricted'}, blobs: 'unrestricted', spectators: 'unrestricted',
      sewn_tadpoles: false, symmetrize_left_right: false},
      symmetrization: {initial: false, final_state: false, external_fermions: false},
      classes: 'Auto',
    },
    numerators: {prefactor: '1', projector: null, naming: {named: {prefix: 'GL'}},
      filter_zero_color: filtered, loop_momentum_bases: {}, named_loop_momentum_bases: {}},
    selection: {selected_diagram_ids: null, selected_diagram_names: null,
      vetoed_diagram_ids: [], vetoed_diagram_names: [], forced_cuts: []},
    grouping: filtered ? 'none' : {up_to_scalar: {numerical_sample_seed: 3,
      number_of_numerical_samples: 5, differentiate_particle_masses_only: true,
      fully_numerical_substitution: false, check_canonical_numerator: false,
      symmetric_polarizations: false}},
    folder_name: 'empty_empty', process_id: 0,
  };
}
const canonicalCard = (buffer) => Object.fromEntries(Object.entries(JSON.parse(buffer)).sort());
assert.deepStrictEqual(canonicalCard(parameters), canonicalCard(bytes(current, 'post-import-parameters.json')),
  'complete cross-producer input card');
const cardOwners = Object.fromEntries([['source43', historical], ['current', current]].map(([route, root]) =>
  [route, liveCards(root, ['post-import', ...rows.map((row) => `${row}-after`)]
    .map((stage) => path.join(root, `${stage}-parameters.json`)))]));
const result = {cardOwners, modelSha256: hash(model), parametersSha256: hash(parameters),
  normalizedSettingsSha256: hash(JSON.stringify(settings)), stages: {}};
for (const [route, root] of [['source43', historical], ['current', current]]) {
  result.stages[route] = {};
  for (const row of rows) {
    const stage = `${row}-after`;
    if (!model.equals(bytes(root, `${stage}-model.json`)) ||
        JSON.stringify(canonicalCard(parameters)) !== JSON.stringify(canonicalCard(bytes(root, `${stage}-parameters.json`))) ||
        JSON.stringify(settings) !== JSON.stringify(normalizeSettings(bytes(root, `${stage}-settings.json`)))) {
      throw new Error(`${route}/${row}: CLIState changed after model import`);
    }
    const definition = bytes(root, `${row}-definition.json`);
    assert.deepStrictEqual(JSON.parse(definition), expectedDefinition(row),
      `${route}/${row}: complete resolved generation policy`);
    if (!definition.equals(bytes(historical, `${row}-definition.json`))) {
      throw new Error(`${route}/${row}: resolved process definition differs between producers`);
    }
    result.stages[route][row] = {
      definitionSha256: hash(definition),
      graphWitnessSha256: hash(bytes(root, `${row}-graphs.tsv`)),
      nativeSourceSha256: row.endsWith('grouped') ? hash(bytes(root, `${row}-native-source.json`)) : null,
      nativeFinalSha256: row.endsWith('grouped') ? hash(bytes(root, `${row}-native-final.json`)) : null,
      nativeGroupsSha256: row.endsWith('grouped') ? hash(bytes(root, `${row}-native-groups.json`)) : null,
    };
  }
}
fs.writeFileSync(output, JSON.stringify(result, null, 2) + '\n');
console.log('actual live post-import state is common and unchanged through all four rows');
