import fs from 'node:fs';
import path from 'node:path';
import assert from 'node:assert/strict';
import {cases, parseDot, canonical, cutChecks, componentData, parseFactorExpression,
  verifyFactor, readCase, matched, mutationControls} from '../analyze-graphs.mjs';

const current = 'target/ghost-weight-certificate/current/states/ghost_weight_certificate';
const output = 'target/ghost-weight-certificate/round2/controls.json';
const results = [];
function control(name, expected, action) {
  let outcome = 'accepted';
  let error = null;
  try { action(); } catch (failure) { outcome = 'rejected'; error = String(failure.message).slice(0, 320); }
  results.push({name, expected, outcome, error});
  assert.equal(outcome, expected, `control failed: ${name}`);
}
function witness(family, label, folder, name, kind) {
  const root = path.join(current, family);
  const row = fs.readFileSync(path.join(root, `${label}-graph-witness.tsv`), 'utf8')
    .trimEnd().split('\n').slice(1).map((line) => line.split('\t')).find((parts) => parts[1] === name);
  assert(row);
  const [index, graphName, factor, weight, group, initial, cuts, thresholds] = row;
  const record = {index: Number(index), name: graphName, factor, weight, group,
    initial: JSON.parse(initial), cuts: JSON.parse(cuts), thresholds: JSON.parse(thresholds)};
  const dotPath = path.join(root, `${label}-dot/processes/${kind}/${folder}/default/${name}.dot`);
  const dot = fs.readFileSync(dotPath, 'utf8');
  return {record, graph: parseDot(dot, kind), dot, dotPath, kind};
}
function validate(item) {
  cutChecks(item.graph, item.record, item.kind);
  const graph = canonical(item.graph, item.record);
  verifyFactor(item.record, componentData(item.graph), graph.automorphisms, item.kind);
  return graph;
}
const electron = witness('epem_bbh', 'bbh_qcd2', 'epem_bbxh', 'GL0000', 'cross_sections');
const vacuum = witness('vacuum_nonfactorizable', 'l2', 'empty_empty', 'GL1', 'amplitudes');
const gluon = witness('vacuum_nonfactorizable', 'l2', 'empty_empty', 'GL0', 'amplitudes');
const baselineElectron = validate(electron);
const baselineVacuum = validate(vacuum);
const factor = vacuum.record.factor;
const invalidFactors = [
  ['trailing coefficient', `${factor}*2`],
  ['trailing bare symbol', `${factor}*Bogus`],
  ['unsupported functional head and noninteger argument', `${factor}*Bogus(x)`],
  ['trailing additive term', `${factor}+7`],
  ['squared complete factor', `(${factor})^2`],
  ['multiply instead of divide by AutG', factor.replace('/AutG(2)', '*AutG(2)')],
  ['zero denominator', `${factor}/0`],
  ['duplicated AutG head with compensating division', `${factor}*AutG(2)/AutG(2)`],
];
for (const [name, changed] of invalidFactors) {
  control(name, 'rejected', () => validate({...vacuum, record: {...vacuum.record, factor: changed}}));
}
for (const [name, equivalent] of [
  ['parenthesized complete factor', `((${factor}))`],
  ['neutral scalar and whitespace', ` ( ${factor} ) * 1 / 1 `],
  ['qualified DOT spelling', '(feynkit_generator_factor::AutG(2))^(-1)*feynkit_generator_factor::ExternalFermionOrderingSign(1)*feynkit_generator_factor::InternalFermionLoopSign(-1)'],
]) {
  control(name, 'accepted', () => validate({...vacuum, record: {...vacuum.record, factor: equivalent}}));
}
control('reversed incoming positron arrow fails locally', 'rejected', () => {
  const changed = structuredClone(electron);
  changed.graph.edges.find((edge) => edge.particle === 'e+' && edge.pu === null).direction = 'forward';
  validate(changed);
});
control('valid whole ghost circuit edge-occurrence reversal', 'accepted', () => {
  const changed = structuredClone(vacuum);
  for (const edge of changed.graph.edges) if (edge.particle === 'ghG') edge.direction = 'back';
  const normalized = validate(changed);
  assert.equal(normalized.key, baselineVacuum.key);
  assert.equal(normalized.automorphisms, baselineVacuum.automorphisms);
});
control('neutral node and half-edge relabel, edge order and cut-list reorder', 'accepted', () => {
  const changed = structuredClone(electron);
  const nodeNames = [...changed.graph.nodes.keys()];
  const nodes = new Map(nodeNames.map((name, index) => [name, `neutral${nodeNames.length - index}`]));
  const halfEdges = [...changed.graph.ports.keys()];
  const hedges = new Map(halfEdges.map((hedge) => [hedge, halfEdges.length - 1 - hedge]));
  changed.graph.nodes = new Map([...changed.graph.nodes].reverse().map(([name, label]) => [nodes.get(name), label]));
  changed.graph.edges = changed.graph.edges.map((edge) => ({...edge,
    u: nodes.get(edge.u), v: nodes.get(edge.v),
    pu: edge.pu === null ? null : hedges.get(edge.pu),
    pv: edge.pv === null ? null : hedges.get(edge.pv)})).reverse();
  changed.graph.ports = new Map();
  changed.graph.seams = new Map();
  for (const [index, edge] of changed.graph.edges.entries()) {
    for (const [hedge, node] of [[edge.pu, edge.u], [edge.pv, edge.v]]) {
      if (hedge !== null) changed.graph.ports.set(hedge, {edge: index, node});
    }
    if (edge.pair !== null) {
      if (!changed.graph.seams.has(edge.pair)) changed.graph.seams.set(edge.pair, []);
      changed.graph.seams.get(edge.pair).push(edge);
    }
  }
  const remap = (selection) => Object.fromEntries(Object.entries(selection)
    .map(([key, values]) => [key, values.map((hedge) => hedges.get(hedge)).reverse()]));
  changed.record.initial = remap(changed.record.initial);
  changed.record.cuts = changed.record.cuts.map(remap).reverse();
  changed.record.thresholds = changed.record.thresholds.map(remap).reverse();
  const neutral = validate(changed);
  assert.deepEqual([neutral.key, neutral.cutKey, neutral.automorphisms],
    [baselineElectron.key, baselineElectron.cutKey, baselineElectron.automorphisms]);
});
control('proper-subgraph threshold is valid', 'accepted', () => {
  const selection = electron.record.thresholds[0];
  assert(selection.left.length + selection.right.length < electron.graph.ports.size);
  assert.equal(validate(electron).cutKey, baselineElectron.cutKey);
});
control('proper-subgraph threshold missing crossing mate fails locally', 'rejected', () => {
  const changed = structuredClone(electron);
  changed.record.thresholds[0].right_cut.pop();
  cutChecks(changed.graph, changed.record, changed.kind);
});
control('duplicate selected cut incidence fails locally', 'rejected', () => {
  const changed = structuredClone(electron);
  changed.record.cuts[0].left_cut.push(changed.record.cuts[0].left_cut[0]);
  cutChecks(changed.graph, changed.record, changed.kind);
});
control('malformed closed ghost arrow fails locally', 'rejected', () => {
  const changed = structuredClone(vacuum);
  changed.graph.edges.find((edge) => edge.particle === 'ghG').direction = 'back';
  validate(changed);
});
control('removed ghost directedness fails locally', 'rejected', () => {
  const changed = structuredClone(vacuum);
  changed.graph.edges.find((edge) => edge.particle === 'ghG').direction = 'none';
  validate(changed);
});
control('unsupported physical self-loop substitution fails parser', 'rejected', () => {
  const changed = vacuum.dot.replace('v1:0\t-> v0:1', 'v1:0\t-> v1:1');
  assert.notEqual(changed, vacuum.dot);
  parseDot(changed, 'amplitudes');
});
control('parallel gluon multiplicity changes directed AutG', 'rejected', () => {
  const changed = structuredClone(gluon);
  changed.graph.edges.pop();
  const graph = canonical(changed.graph, changed.record);
  verifyFactor(changed.record, componentData(changed.graph), graph.automorphisms, changed.kind);
});
control('wrong AutG head argument', 'rejected', () => validate({...vacuum,
  record: {...vacuum.record, factor: factor.replace('AutG(2)', 'AutG(4)')}}));
control('external role mutation fails locally', 'rejected', () => {
  const changed = structuredClone(electron);
  const entry = [...changed.graph.nodes.values()].find((label) => label[0] === 'X' && label[1] === 'in');
  entry[1] = 'out';
  validate(changed);
});
control('external-order sign head mutation', 'rejected', () => validate({...electron,
  record: {...electron.record, factor: electron.record.factor.replace('ExternalFermionOrderingSign(-1)',
    'ExternalFermionOrderingSign(1)')}}));
control('incoming-antifermion sign head mutation', 'rejected', () => validate({...electron,
  record: {...electron.record, factor: electron.record.factor.replace('AntiFermionSpinSumSign(-1)',
    'AntiFermionSpinSumSign(1)')}}));
control('initial-state side mutation', 'rejected', () => {
  const changed = structuredClone(electron);
  changed.record.initial.left[0] = changed.record.initial.right[0];
  validate(changed);
});
control('selected physical cut missing incidence', 'rejected', () => {
  const changed = structuredClone(electron);
  changed.record.cuts[0].left_cut.pop();
  validate(changed);
});
control('recorded per-graph rational weight flip', 'rejected', () => validate({...vacuum,
  record: {...vacuum.record, weight: '1/2'}}));
control('raw expression parser consumes all tokens', 'rejected', () => parseFactorExpression(`${factor}*Bogus`));
const configuration = cases.find((entry) => entry.id === 'vacuum_l2');
const witnessPath = path.join(current, 'vacuum_nonfactorizable/l2-graph-witness.tsv');
const originalWitness = fs.readFileSync(witnessPath, 'utf8');
const virtualFiles = (overrides) => ({...fs, readFileSync(location, ...args) {
  const content = overrides.get(path.resolve(String(location)));
  return content === undefined ? fs.readFileSync(location, ...args)
    : typeof args[0] === 'string' ? content : Buffer.from(content);
}});
control('mirrored TSV and DOT factor corruption fails full readCase without hash guards', 'rejected', () => {
  const lines = originalWitness.trimEnd().split('\n');
  const fields = lines[2].split('\t');
  fields[2] += '*2';
  lines[2] = fields.join('\t');
  const modifiedDot = vacuum.dot.replace(/(overall_factor = ")([^"]+)(";)/,
    (_, prefix, text, suffix) => `${prefix}${text}*2${suffix}`);
  const overrides = new Map([[path.resolve(witnessPath), `${lines.join('\n')}\n`],
    [path.resolve(vacuum.dotPath), modifiedDot]]);
  readCase(current, configuration, 'l2', virtualFiles(overrides));
});
control('DOT-only raw factor corruption fails full readCase', 'rejected', () => {
  const changedDot = vacuum.dot.replace(/(overall_factor = ")([^"]+)(";)/,
    (_, prefix, text, suffix) => `${prefix}${text}*2${suffix}`);
  readCase(current, configuration, 'l2', virtualFiles(new Map([[path.resolve(vacuum.dotPath), changedDot]])));
});
control('DOT-only evaluated scalar corruption fails full readCase', 'rejected', () => {
  const changedDot = vacuum.dot.replace('overall_factor_evaluated = "-1/2"',
    'overall_factor_evaluated = "-2"');
  assert.notEqual(changedDot, vacuum.dot);
  readCase(current, configuration, 'l2', virtualFiles(new Map([[path.resolve(vacuum.dotPath), changedDot]])));
});
control('same-count and same-sum physical replacement fails full join', 'rejected', () => {
  const lines = originalWitness.trimEnd().split('\n');
  const replacement = lines[2].split('\t');
  const removed = lines[3].split('\t');
  assert.equal(replacement[3], removed[3]);
  replacement[0] = removed[0];
  replacement[1] = removed[1];
  lines[3] = replacement.join('\t');
  const removedDot = witness('vacuum_nonfactorizable', 'l2', 'empty_empty', 'GL2', 'amplitudes');
  const overrides = new Map([[path.resolve(witnessPath), `${lines.join('\n')}\n`],
    [path.resolve(removedDot.dotPath), vacuum.dot.replace(/^digraph GL1/, 'digraph GL2')]]);
  const original = readCase(current, configuration, 'l2');
  const changed = readCase(current, configuration, 'l2', virtualFiles(overrides));
  matched(original, changed, 'vacuum_l2', 'same-count-sum-physical-replacement');
});
const original = mutationControls(current);
assert.equal(original.filter((entry) => entry.outcome === 'rejected').length, 6);
assert.equal(original.filter((entry) => entry.outcome === 'accepted').length, 1);
fs.writeFileSync(output, `${JSON.stringify({schema: 1, round1ReviewerFactorReproducer:
  'target/ghost-weight-certificate/review-astra-scope-b-8wU1pd/factor-repro.json',
  hashGuards: false, results, originalControls: original}, null, 2)}\n`);
console.log(`controls: ${results.length} passed (${results.filter((entry) => entry.outcome === 'rejected').length} rejected, ${results.filter((entry) => entry.outcome === 'accepted').length} accepted); original seven passed`);
