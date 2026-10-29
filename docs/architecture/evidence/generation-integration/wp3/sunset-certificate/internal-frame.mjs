import fs from 'node:fs';
import {liveRule, liveModel} from './live-rule.mjs';
import {complete} from './factor.mjs';
import {checkAggregate, nativeCliExpression} from './tensor-product.mjs';
import path from 'node:path';
import crypto from 'node:crypto';

const [historical, current, historicalCapture, currentCapture, output] = process.argv.slice(2);
if (!historical || !current || !historicalCapture || !currentCapture || !output) {
  throw new Error('usage: node internal-frame.mjs <historical-frames> <current-frames> <historical-capture> <current-capture> <output.json>');
}
const hash = (file) => crypto.createHash('sha256').update(fs.readFileSync(file)).digest('hex');
const equal = (left, right, where) => {
  if (JSON.stringify(left) !== JSON.stringify(right)) throw new Error(`${where}: exact mismatch`);
};
const swapHedge = (source) => source.replaceAll('hedge(0,1)', 'hedge(X,1)')
  .replaceAll('hedge(1,1)', 'hedge(0,1)').replaceAll('hedge(X,1)', 'hedge(1,1)');
const metricHead = 'spenso::{symmetric,linear,real,spenso::tensor}::g(';

function normalizedSymmetricMetrics(source) {
  let result = source;
  let position = 0;
  while ((position = result.indexOf(metricHead, position)) >= 0) {
    const start = position + metricHead.length;
    let parentheses = 0;
    let braces = 0;
    let comma = -1;
    let end = -1;
    for (let index = start; index < result.length; index++) {
      const character = result[index];
      if (character === '(') parentheses++;
      else if (character === ')') {
        if (parentheses === 0 && braces === 0) {end = index; break;}
        parentheses--;
      } else if (character === '{') braces++;
      else if (character === '}') braces--;
      else if (character === ',' && parentheses === 0 && braces === 0) comma = index;
    }
    if (comma < 0 || end < 0) throw new Error('incomplete symmetric metric expression');
    const argumentsInOrder = [result.slice(start, comma), result.slice(comma + 1, end)].sort();
    const argumentsText = argumentsInOrder.join(',');
    result = result.slice(0, start) + argumentsText + result.slice(end);
    position = start + argumentsText.length + 1;
  }
  return result;
}

function splitTopLevel(source, separator) {
  let parentheses = 0;
  let braces = 0;
  let start = 0;
  const parts = [];
  for (let index = 0; index < source.length; index++) {
    const character = source[index];
    if (character === '(') parentheses++;
    else if (character === ')') parentheses--;
    else if (character === '{') braces++;
    else if (character === '}') braces--;
    else if (character === separator && parentheses === 0 && braces === 0) {
      parts.push(source.slice(start, index));
      start = index + 1;
    }
    if (parentheses < 0 || braces < 0) throw new Error('unbalanced factorized expression');
  }
  if (parentheses !== 0 || braces !== 0) throw new Error('unbalanced factorized expression');
  parts.push(source.slice(start));
  return parts;
}

function signedTopMassSum(source) {
  if (!source.startsWith('(')) throw new Error('top edge lacks its local mass/slash sum');
  let depth = 0;
  let end = -1;
  for (let index = 0; index < source.length; index++) {
    if (source[index] === '(') depth++;
    else if (source[index] === ')' && --depth === 0) {end = index; break;}
  }
  if (end < 0) throw new Error('top local sum did not close');
  const terms = splitTopLevel(source.slice(1, end), '+');
  if (terms.length !== 2 || terms.filter((term) => term.includes('::Q(0,')).length !== 1 ||
      terms.filter((term) => term.includes('::MT*')).length !== 1) {
    throw new Error('top local sum is not the exact mass-plus-slash pair');
  }
  return '(' + terms.map((term) => term.includes('::Q(0,') ? `-1*${term}` : term)
    .sort().join('+') + ')' + source.slice(end + 1);
}

function transported(source, graph, owner, index) {
  let result = normalizedSymmetricMetrics(
    ['GL1', 'GL2', 'GL3', 'GL4'].includes(graph) ? swapHedge(source) : source);
  if (graph === 'GL1' && owner === 'vertex' && index === 0) result = `-1*${result}`;
  if (graph === 'GL1' && owner === 'aggregate') {
    if (!result.startsWith('(-𝑖)*')) throw new Error('ghost aggregate scalar differs');
    result = '(𝑖)*' + result.slice('(-𝑖)*'.length);
  }
  if (['GL2', 'GL4'].includes(graph) &&
      (owner === 'aggregate' || owner === 'edge' && index === 0)) {
    if (!result.startsWith('(𝑖)*')) throw new Error('light fermion scalar differs');
    result = '(-𝑖)*' + result.slice('(𝑖)*'.length);
  }
  if (graph === 'GL3' && (owner === 'aggregate' || owner === 'edge' && index === 0)) {
    result = signedTopMassSum(result);
  }
  return result;
}

function fields(file) {
  const result = {vertex_map: {comparison: [], representative: []},
    edge_map: {comparison: [], representative: []}, final_edge_map: [],
    graph_vertex: {source: [], comparison: [], representative: []},
    graph_edge: {source: [], comparison: [], representative: []}};
  for (const line of fs.readFileSync(file, 'utf8').trimEnd().split('\n')) {
    const [kind, ...parts] = line.split('\t');
    if (['vertex_map', 'edge_map', 'graph_vertex', 'graph_edge'].includes(kind)) {
      const [stage, ...data] = parts;
      result[kind][stage].push(data);
    } else if (kind === 'final_edge_map') result[kind].push(parts.map(Number));
    else if (Object.hasOwn(result, kind)) throw new Error(`duplicate frame field ${kind}`);
    else result[kind] = parts.join('\t');
  }
  return result;
}

function diagram(file) {return JSON.parse(fs.readFileSync(file, 'utf8'));}
function vector(basis, edge) {
  const signature = basis.edge_signatures[String(edge)];
  if (!signature || signature.external.length) throw new Error('missing exact vacuum momentum signature');
  return signature.loops;
}

function checkBasis(source, target, edgeMap, where) {
  const a = source.loop_momentum_basis;
  const b = target.loop_momentum_basis;
  if (a.external_edges.length || b.external_edges.length ||
      a.dependent_externals.length || b.dependent_externals.length ||
      a.loop_edges.length !== 2 || b.loop_edges.length !== 2) {
    throw new Error(`${where}: wrong signed affine dimension or shift`);
  }
  const entries = edgeMap.map(([old, next, flipped]) =>
    [Number(old), Number(next), Number(flipped) ? -1 : 1]);
  equal(entries.map(([old]) => old).sort(), [...source.edges.keys()], `${where}: source edge map`);
  equal(entries.map(([, next]) => next).sort(), [...target.edges.keys()], `${where}: target edge map`);
  const byOld = new Map(entries.map((entry) => [entry[0], entry]));
  const matrix = a.loop_edges.map((edge) => {
    const [, next, sign] = byOld.get(edge);
    return vector(b, next).map((coefficient) => sign * coefficient);
  });
  for (const [old, next, sign] of entries) {
    const transformed = [0, 1].map((column) => vector(a, old)
      .reduce((sum, coefficient, row) => sum + coefficient * matrix[row][column], 0));
    equal(transformed, vector(b, next).map((coefficient) => sign * coefficient),
      `${where}: signed affine edge ${old}`);
  }
  const determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
  if (Math.abs(determinant) !== 1) throw new Error(`${where}: frame map is not unimodular`);
  return matrix;
}

function checkGraphRecords(frame, source, comparison, representative, model, fingerprint, where) {
  equal(frame.graph_vertex.source, frame.graph_vertex.representative,
    `${where}: provisional representative vertices`);
  equal(frame.graph_edge.source, frame.graph_edge.representative,
    `${where}: provisional representative edges`);
  for (const [stage, diagram] of [['source', source], ['comparison', comparison],
    ['representative', representative]]) {
    complete(diagram.overall_factor);
    checkAggregate(diagram, nativeCliExpression, `${where}/${stage}`);
    equal(diagram.model, fingerprint, `${where}/${stage}: actual post-import model owner`);
    const vertices = frame.graph_vertex[stage];
    const edges = frame.graph_edge[stage];
    if (vertices.length !== diagram.vertices.length || edges.length !== diagram.edges.length) {
      throw new Error(`${where}: provisional ${stage} graph size differs`);
    }
    for (const [index, [id, rule]] of vertices.entries()) {
      const match = rule.match(/^VertexRuleId\((\d+)\)$/);
      if (Number(id) !== index || !match || Number(match[1]) !== diagram.vertices[index].interaction ||
          !model.vertex_rules[Number(match[1])]) {
        throw new Error(`${where}: provisional ${stage} vertex rule differs`);
      }
    }
    for (const [index, [id, from, to, particle, directed]] of edges.entries()) {
      const match = particle.match(/^ParticleId\((\d+)\)$/);
      const [endpoints, edge] = diagram.edges[index];
      if (Number(id) !== index || !match || !model.particles[Number(match[1])] ||
          Number(match[1]) !== edge.particle || Number(directed) !== Number(edge.directed) ||
          from !== `v:${endpoints.source}` || to !== `v:${endpoints.target}`) {
        throw new Error(`${where}: provisional ${stage} oriented edge differs`);
      }
      const carrier = model.particles[edge.particle];
      const sourceRule = model.vertex_rules[diagram.vertices[endpoints.source].interaction];
      const targetRule = model.vertex_rules[diagram.vertices[endpoints.target].interaction];
      if (!Number.isSafeInteger(edge.source_slot) || !Number.isSafeInteger(edge.target_slot) ||
          sourceRule.particles[edge.source_slot] !== carrier.antiname ||
          targetRule.particles[edge.target_slot] !== carrier.name) {
        throw new Error(`${where}/${stage}/edge/${index}: live rule half-edge slot or orientation differs`);
      }
    }
    equal(diagram.half_edge_order,
      diagram.edges.flatMap((_, edge) => [{edge, endpoint: 'source'},
        {edge, endpoint: 'target'}]), `${where}: finalized half-edge owners`);
  }
  for (const [old, next, flipped] of frame.edge_map.comparison.map((row) => row.map(Number))) {
    const a = frame.graph_edge.source[old];
    const b = frame.graph_edge.comparison[next];
    const oldParticle = Number(a[3].match(/^ParticleId\((\d+)\)$/)?.[1]);
    const newParticle = Number(b[3].match(/^ParticleId\((\d+)\)$/)?.[1]);
    const expectedName = flipped ? model.particles[oldParticle].antiname :
      model.particles[oldParticle].name;
    if (b[1] !== (flipped ? a[2] : a[1]) ||
        b[2] !== (flipped ? a[1] : a[2]) ||
        model.particles[newParticle].name !== expectedName || b[4] !== a[4]) {
      throw new Error(`${where}: comparison reversed edge or antiparticle owner differs`);
    }
    const sourceEdge = source.edges[old][1];
    const targetEdge = comparison.edges[next][1];
    if (targetEdge.source_slot !== (flipped ? sourceEdge.target_slot : sourceEdge.source_slot) ||
        targetEdge.target_slot !== (flipped ? sourceEdge.source_slot : sourceEdge.target_slot) ||
        targetEdge.orientation !== (flipped ?
          sourceEdge.orientation === 'Default' ? 'Reversed' : 'Default' :
          sourceEdge.orientation)) {
      throw new Error(`${where}: comparison half-edge slot bijection differs`);
    }
  }
}

function checkOwner(source, comparison, graph, owner, index, where) {
  const lhs = source.numerator;
  const rhs = comparison.numerator;
  const transformed = transported(lhs, graph, owner, index);
  const expected = normalizedSymmetricMetrics(rhs);
  const left = owner === 'aggregate' ? splitTopLevel(transformed, '*').sort() : [transformed];
  const right = owner === 'aggregate' ? splitTopLevel(expected, '*').sort() : [expected];
  equal(left, right, `${where}: factor-preserving ${owner} ${index} tensor`);
  const leftAxes = source.num_interface?.map(swapHedge) ?? [];
  const rightAxes = comparison.num_interface ?? [];
  const permutation = leftAxes.map((axis) => rightAxes.indexOf(axis));
  if (new Set(leftAxes).size !== leftAxes.length ||
      new Set(rightAxes).size !== rightAxes.length ||
      permutation.includes(-1) || new Set(permutation).size !== rightAxes.length) {
    throw new Error(`${where}: complete logical ports ${owner} ${index}: bijection`);
  }
  const expectedPermutation = leftAxes.map((_, position) =>
    owner === 'edge' && index === 0 && ['GL1', 'GL2', 'GL3', 'GL4'].includes(graph) &&
      position < 2 ? 1 - position : position);
  equal(permutation, expectedPermutation,
    `${where}: logical port layout is inconsistent with edge reversal and tensor transport`);
  return permutation;
}

const frameRoots = {source43: historical, current};
const captureRoots = {source43: historicalCapture, current: currentCapture};
const result = {routes: {}};
for (const [route, root] of Object.entries(frameRoots)) {
  result.routes[route] = {};
  const model = JSON.parse(fs.readFileSync(path.join(captureRoots[route], 'post-import-model.json')));
  const fingerprintHex = liveModel(captureRoots[route]);
  if (!/^[0-9a-f]{64}$/.test(fingerprintHex)) throw new Error('invalid captured model fingerprint');
  const fingerprint = [...Buffer.from(fingerprintHex, 'hex')];
  for (const [filter, first, replay, count] of [['all', 0, 1, 6],
    ['nonfactorizable', 2, 3, 5]]) {
    const native = JSON.parse(fs.readFileSync(path.join(captureRoots[route],
      `${filter}-grouped-native-source.json`)));
    const records = [];
    for (let id = 0; id < count; id++) {
      const prefix = (run) => path.join(root, `run-${run}`, `source-${id}`);
      const primary = prefix(first);
      const duplicate = prefix(replay);
      for (const stage of ['source', 'comparison', 'representative']) {
        const name = `${stage}.json`;
        if (!fs.readFileSync(`${primary}-${name}`).equals(fs.readFileSync(`${duplicate}-${name}`))) {
          throw new Error(`${route}/${filter}/${id}: CLI and native replay ${stage} differ`);
        }
      }
      const mainText = fs.readFileSync(`${primary}.txt`, 'utf8')
        .replace(/^run_index\t\d+/m, 'run_index\tX');
      const replayText = fs.readFileSync(`${duplicate}.txt`, 'utf8')
        .replace(/^run_index\t\d+/m, 'run_index\tX');
      if (mainText !== replayText) throw new Error(`${route}/${filter}/${id}: CLI/native frame maps differ`);
      const frame = fields(`${primary}.txt`);
      const source = diagram(`${primary}-source.json`);
      const comparison = diagram(`${primary}-comparison.json`);
      const representative = diagram(`${primary}-representative.json`);
      const where = `${route}/${filter}/${source.name}`;
      if (frame.source_index !== String(id) || frame.source_name !== source.name ||
          source.name !== native[id].name ||
          JSON.stringify(source) !== JSON.stringify(native[id]) ||
          !fs.readFileSync(`${primary}-source.json`).equals(
            fs.readFileSync(`${primary}-representative.json`))) {
        throw new Error(`${where}: original source and selected representative are not bound`);
      }
      equal(frame.vertex_map.comparison.map((row) => row.map(Number)),
        source.vertices.map((_, vertex) => [vertex, vertex]), `${where}: comparison vertex map`);
      equal(frame.vertex_map.representative.map((row) => row.map(Number)),
        source.vertices.map((_, vertex) => [vertex, vertex]), `${where}: representative vertex map`);
      equal(frame.final_edge_map, source.edges.map((_, edge) => [edge, edge]),
        `${where}: final edge owner map`);
      const comparisonMap = frame.edge_map.comparison;
      const representativeMap = frame.edge_map.representative;
      const distinct = ['GL1', 'GL2', 'GL3', 'GL4'].includes(source.name);
      if (Number(frame.distinct) !== Number(distinct) ||
          comparisonMap.length !== source.edges.length ||
          representativeMap.length !== source.edges.length) {
        throw new Error(`${where}: pair kind or edge map differs`);
      }
      equal(comparisonMap.map((row) => row.map(Number)),
        source.edges.map((_, edge) => [edge, edge, Number(distinct && edge === 0)]),
        `${where}: complete comparison edge map`);
      equal(representativeMap.map((row) => row.map(Number)),
        source.edges.map((_, edge) => [edge, edge, 0]), `${where}: representative edge map`);
      checkGraphRecords(frame, source, comparison, representative, model, fingerprint, where);
      if (frame.requested_prefactor !== '1' || frame.comparison_prefactor !== '1' ||
          frame.representative_prefactor !== '1' || source.numerator_prefactor !== '1' ||
          comparison.numerator_prefactor !== '1' || source.projector !== '1' ||
          comparison.projector !== '1' || source.overall_factor !== comparison.overall_factor ||
          source.symmetry_factor !== comparison.symmetry_factor ||
          source.cuts.length || comparison.cuts.length ||
          source.topology_threshold_candidates.length || comparison.topology_threshold_candidates.length) {
        throw new Error(`${where}: requested prefactor, cuts or source weight differ`);
      }
      const matrix = checkBasis(source, comparison, comparisonMap, `${where}/comparison`);
      const representativeMatrix = checkBasis(source, representative, representativeMap,
        `${where}/representative`);
      equal(representativeMatrix, [[1, 0], [0, 1]], `${where}: representative frame`);
      equal(matrix, distinct ? [[-1, 0], [0, 1]] : [[1, 0], [0, 1]],
        `${where}: observed signed comparison frame`);
      const logicalPortPermutations = {vertices: [], edges: []};
      for (const [index, vertex] of source.vertices.entries()) {
        logicalPortPermutations.vertices.push(checkOwner(vertex, comparison.vertices[index],
          source.name, 'vertex', index, where));
      }
      for (const [index, [, edge]] of source.edges.entries()) {
        logicalPortPermutations.edges.push(checkOwner(edge, comparison.edges[index][1],
          source.name, 'edge', index, where));
      }
      checkOwner(source, comparison, source.name, 'aggregate', 0, where);
      records.push({sourceIndex: id, sourceName: source.name, distinct, matrix,
        logicalPortPermutations,
        comparisonEdgeMap: comparisonMap.map((row) => row.map(Number)),
        sourceSha256: hash(`${primary}-source.json`),
        comparisonSha256: hash(`${primary}-comparison.json`),
        representativeSha256: hash(`${primary}-representative.json`),
        frameSha256: hash(`${primary}.txt`)});
    }
    if (count > 4) {
      const file = (id) => path.join(root, `run-${first}`, `source-${id}-comparison.json`);
      const u = diagram(file(2));
      const d = diagram(file(4));
      equal(u.vertices.map((vertex) => [vertex.numerator, vertex.num_interface]),
        d.vertices.map((vertex) => [vertex.numerator, vertex.num_interface]),
        `${route}/${filter}: comparison d/u vertex tensors`);
      equal(u.edges.map(([, edge]) => [edge.numerator, edge.num_interface]),
        d.edges.map(([, edge]) => [edge.numerator, edge.num_interface]),
        `${route}/${filter}: comparison d/u edge tensors`);
      if (u.numerator !== d.numerator) {
        throw new Error(`${route}/${filter}: comparison d/u full factorized tensor differs`);
      }
    }
    result.routes[route][filter] = records;
  }
}
for (const filter of ['all', 'nonfactorizable']) {
  const old = result.routes.source43[filter];
  const fresh = result.routes.current[filter];
  equal(old.map((row) => [row.sourceName, row.matrix, row.comparisonEdgeMap]),
    fresh.map((row) => [row.sourceName, row.matrix, row.comparisonEdgeMap]),
    `${filter}: source43/current physical frame maps`);
}
result.liveRules = Object.fromEntries(Object.entries(frameRoots).map(([route, root]) => {
  const files = fs.readdirSync(root, {recursive: true}).filter((file) => file.endsWith('.json'))
    .map((file) => path.join(root, file)).sort();
  return [route, liveRule(captureRoots[route], files)];
}));
fs.writeFileSync(output, JSON.stringify(result, null, 2) + '\n');
console.log('all original comparison and representative frames, local factors and d/u comparison tensors verified');
