import fs from 'node:fs';
import path from 'node:path';
import assert from 'node:assert/strict';
import {cases, parseDot, canonical, cutChecks, componentData, parseFactorExpression,
  verifyFactor, readCase, matched} from '../analyze-graphs.mjs';

const full = JSON.parse(fs.readFileSync('target/ghost-weight-certificate/round3/full-physical-joins.json'));
const resultPath = 'target/ghost-weight-certificate/round3/controls.json';
const results = [];
function control(name, expected, action) {
  let outcome = 'accepted';
  let error = null;
  let detail = null;
  try { detail = action(); } catch (failure) { outcome = 'rejected'; error = String(failure.message).slice(0, 400); }
  results.push({name, expected, outcome, error, detail});
  assert.equal(outcome, expected, `control failed: ${name}`);
}
const configuration = cases.find((entry) => entry.id === 'vacuum_l2');
function data(route, name = 'GL1') {
  const root = full[`${route}Root`];
  const family = path.join(root, 'vacuum_nonfactorizable');
  const witnessPath = path.join(family, 'l2-graph-witness.tsv');
  const source = fs.readFileSync(witnessPath, 'utf8');
  const fields = source.trimEnd().split('\n').slice(1).map((line) => line.split('\t'))
    .find((row) => row[1] === name);
  assert(fields);
  const [index, graphName, factor, weight, group, initial, cuts, thresholds] = fields;
  const record = {index: Number(index), name: graphName, factor, weight, group,
    initial: JSON.parse(initial), cuts: JSON.parse(cuts), thresholds: JSON.parse(thresholds)};
  const dotPath = full.cases['vacuum_l2/l2'][route].find((row) => row.name === name).dotPath;
  const dot = fs.readFileSync(dotPath, 'utf8');
  const graph = parseDot(dot, configuration.kind);
  return {route, root, witnessPath, source, record, dotPath, dot, graph};
}
const current = data('current');
const historical = data('historical');
const suffixes = ['*-2^0', '/-2^0'];
const virtualFiles = (overrides) => ({...fs, readFileSync(location, ...args) {
  const text = overrides.get(path.resolve(String(location)));
  return text === undefined ? fs.readFileSync(location, ...args)
    : typeof args[0] === 'string' ? text : Buffer.from(text);
}});
function overlay(witness, suffix, lane) {
  const changed = new Map();
  if (lane !== 'DOT') {
    const lines = witness.source.trimEnd().split('\n');
    const index = lines.findIndex((line) => line.split('\t')[1] === witness.record.name);
    const fields = lines[index].split('\t');
    fields[2] += suffix;
    lines[index] = fields.join('\t');
    changed.set(path.resolve(witness.witnessPath), `${lines.join('\n')}\n`);
  }
  if (lane !== 'TSV') {
    const dot = witness.dot.replace(/(overall_factor = ")([^"\n]+)(";)/,
      (_, prefix, factor, end) => `${prefix}${factor}${suffix}${end}`);
    assert.notEqual(dot, witness.dot);
    changed.set(path.resolve(witness.dotPath), dot);
  }
  return virtualFiles(changed);
}
const original = current.record.factor;
const physical = componentData(current.graph);
const aut = canonical(current.graph, current.record).automorphisms;
for (const suffix of suffixes) {
  control(`raw parser gives Symbolica value for ${suffix}`, 'accepted', () => {
    const value = parseFactorExpression(`${original}${suffix}`).value;
    assert.deepEqual(value, {n: 1n, d: 2n});
    return '1/2 rather than recorded -1/2';
  });
  control(`local factor checker rejects ${suffix}`, 'rejected', () =>
    verifyFactor({...current.record, factor: `${original}${suffix}`}, physical, aut, 'amplitudes'));
  for (const route of [historical, current]) for (const lane of ['TSV', 'DOT', 'both']) {
    control(`${route.route} full readCase ${lane} ${suffix}`, 'rejected', () =>
      readCase(route.root, configuration, 'l2', overlay(route, suffix, lane)));
  }
  control(`mirrored historical/current full join ${suffix}`, 'rejected', () => {
    const old = readCase(historical.root, configuration, 'l2', overlay(historical, suffix, 'both'));
    const next = readCase(current.root, configuration, 'l2', overlay(current, suffix, 'both'));
    matched(old, next, 'vacuum_l2', 'mirrored-invalid-power');
  });
}
for (const [name, expression, expected] of [
  ['negative base parenthesized zero power', `${original}*(-2)^0`, '-1/2'],
  ['negative base nested parentheses', `${original}*((-2)^0)`, '-1/2'],
  ['explicit negative power value', `${original}*(-2)^(-1)*(-2)`, '-1/2'],
]) control(name, 'accepted', () => {
  const value = parseFactorExpression(expression).value;
  const text = `${value.n}/${value.d}`;
  assert.equal(text, expected);
  verifyFactor({...current.record, factor: expression}, physical, aut, 'amplitudes');
  return text;
});
const electronConfig = cases.find((entry) => entry.id === 'epem_bbh_qcd2');
const electronRoot = full.currentRoot;
const electronWitness = path.join(electronRoot, electronConfig.family, `${electronConfig.label}-graph-witness.tsv`);
const fields = fs.readFileSync(electronWitness, 'utf8').split('\n')[1].split('\t');
const electronRecord = {index: Number(fields[0]), name: fields[1], factor: fields[2], weight: fields[3],
  group: fields[4], initial: JSON.parse(fields[5]), cuts: JSON.parse(fields[6]), thresholds: JSON.parse(fields[7])};
const electronDot = full.cases[`${electronConfig.id}/${electronConfig.label}`].current[0].dotPath;
const electronGraph = parseDot(fs.readFileSync(electronDot, 'utf8'), electronConfig.kind);
const originalKey = canonical(electronGraph, electronRecord).cutKey;
control('initial left duplicate occurrence', 'rejected', () => {
  const changed = structuredClone(electronRecord);
  changed.initial.left.push(changed.initial.left[0]);
  cutChecks(electronGraph, changed, electronConfig.kind);
});
control('initial right duplicate occurrence', 'rejected', () => {
  const changed = structuredClone(electronRecord);
  changed.initial.right.push(changed.initial.right[0]);
  cutChecks(electronGraph, changed, electronConfig.kind);
});
for (const side of ['left', 'right']) {
  control(`initial ${side} unknown half-edge`, 'rejected', () => {
    const changed = structuredClone(electronRecord);
    changed.initial[side][0] = 123456;
    cutChecks(electronGraph, changed, electronConfig.kind);
  });
  control(`initial ${side} missing seam mate`, 'rejected', () => {
    const changed = structuredClone(electronRecord);
    changed.initial[side].pop();
    cutChecks(electronGraph, changed, electronConfig.kind);
  });
  control(`initial ${side} wrong in-domain mate`, 'rejected', () => {
    const changed = structuredClone(electronRecord);
    const taken = new Set([...changed.initial.left, ...changed.initial.right]);
    changed.initial[side][0] = [...electronGraph.ports.keys()].find((hedge) => !taken.has(hedge));
    cutChecks(electronGraph, changed, electronConfig.kind);
  });
}
control('initial seam-side reordering remains neutral', 'accepted', () => {
  const changed = structuredClone(electronRecord);
  changed.initial.left.reverse();
  changed.initial.right.reverse();
  cutChecks(electronGraph, changed, electronConfig.kind);
  assert.equal(canonical(electronGraph, changed).cutKey, originalKey);
});
const prior = JSON.parse(fs.readFileSync('target/ghost-weight-certificate/round3/prior-controls-replay.json'));
assert.equal(prior.originalSeven, 7);
assert.equal(prior.round2Controls, 33);
assert.equal(prior.reviewerPriorControls, 31);
assert.equal(prior.reviewerExtraControls, 111);
assert(results.every((entry) => entry.expected === entry.outcome));
fs.writeFileSync(resultPath, `${JSON.stringify({schema: 1, hashGuardsForMutations: false,
  priorSuite: {originalSeven: 7, round2: 33, reviewerPrior: 31, reviewerExtra: 111},
  results}, (_, value) => typeof value === 'bigint' ? String(value) : value, 2)}\n`);
console.log(`round3 controls: ${results.length} passed (${results.filter((entry) => entry.outcome === 'rejected').length} rejected, ${results.filter((entry) => entry.outcome === 'accepted').length} accepted)`);
