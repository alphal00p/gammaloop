import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';

const [currentRoot, historicalRoot, fullOutput, summaryOutput] = process.argv.slice(2);
const model = JSON.parse(fs.readFileSync('docs/architecture/evidence/generation-integration/wp3/pocket-weight-diagnosis/live-post-import-model.json'));
const particle = new Map(model.particles.map((entry) => [entry.name, entry]));
const sha = (value) => crypto.createHash('sha256').update(value).digest('hex');
const encode = JSON.stringify;
const sorted = (items) => [...items].sort((a, b) => encode(a).localeCompare(encode(b), 'en'));
const equal = (a, b, description) => assert.deepEqual(a, b, description);
const cases = [
  {id: 'epem_bbh_qcd2', family: 'epem_bbh', label: 'bbh_qcd2', kind: 'cross_sections', folder: 'epem_bbxh'},
  {id: 'slow_aa_ttx_qcd2', family: 'slow_aa_ttx', label: 'qcd2', kind: 'cross_sections', folder: 'aa_ttx'},
  {id: 'vacuum_l2', family: 'vacuum_nonfactorizable', label: 'l2', kind: 'amplitudes', folder: 'empty_empty'},
  {id: 'vacuum_l3', family: 'vacuum_nonfactorizable', label: 'l3', kind: 'amplitudes', folder: 'empty_empty'},
  {id: 'vacuum_l4', family: 'vacuum_nonfactorizable', label: 'l4', kind: 'amplitudes', folder: 'empty_empty'},
];
function attributes(source) {
  const result = {};
  for (const match of source.matchAll(/(\w+)=((?:"(?:[^"\\]|\\.)*")|[^\s\]]+)/g)) {
    assert(!(match[1] in result), 'duplicate DOT attribute');
    result[match[1]] = match[2].startsWith('"') ? JSON.parse(match[2]) : match[2];
  }
  return result;
}
function parseDot(source, expectedKind) {
  const nodes = new Map();
  const edges = [];
  const ports = new Map();
  for (const line of source.split('\n')) {
    let match = line.match(/^\s*(\w+)(?::(\d+))?\s*->\s*(\w+)(?::(\d+))?\s*\[(.*)\];$/);
    if (match) {
      const attr = attributes(match[5]);
      assert(attr.particle && attr.id !== undefined, 'unparsed DOT edge');
      const edge = {
        u: match[1], v: match[3],
        pu: match[2] === undefined ? null : Number(match[2]),
        pv: match[4] === undefined ? null : Number(match[4]),
        particle: attr.particle,
        direction: attr.dir ?? 'forward',
        pair: attr.is_cut ?? null,
        momentum: attr.lmb_rep ?? null,
        id: Number(attr.id),
      };
      assert(edge.u !== edge.v, 'self-loop needs separate endpoint treatment');
      assert(particle.has(edge.particle), `particle absent from live model: ${edge.particle}`);
      edges.push(edge);
      continue;
    }
    match = line.match(/^\s*(\w+)\s*\[(.*)\];$/);
    if (match) {
      const attr = attributes(match[2]);
      assert(!nodes.has(match[1]), 'duplicate DOT node');
      nodes.set(match[1], attr.int_id ? ['V', attr.int_id] : null);
      continue;
    }
    assert(!/^\s*\w+(?::\d+)?\s*->/.test(line), 'unparsed DOT edge syntax');
  }
  assert.equal(new Set(edges.map((edge) => edge.id)).size, edges.length);
  for (const [index, edge] of edges.entries()) {
    for (const [hedge, node] of [[edge.pu, edge.u], [edge.pv, edge.v]]) {
      if (hedge === null) continue;
      assert(!ports.has(hedge), 'duplicate half-edge');
      ports.set(hedge, {edge: index, node});
    }
  }
  const seams = new Map();
  for (const edge of edges) {
    if (edge.pair === null) continue;
    if (!seams.has(edge.pair)) seams.set(edge.pair, []);
    seams.get(edge.pair).push(edge);
  }
  for (const group of seams.values()) {
    assert.equal(group.length, 2, 'sewn external pair must have two edges');
    const incoming = group.find((edge) => edge.pu === null);
    const outgoing = group.find((edge) => edge.pv === null);
    assert(incoming && outgoing, 'sewn external endpoints missing');
    const externalIndex = incoming.momentum?.match(/^P\((\d+),a___\)$/);
    assert(externalIndex, 'fixed external momentum label missing');
    equal(incoming.particle, outgoing.particle, 'sewn external particle mismatch');
    equal(Number(incoming.pair), outgoing.pu, 'sewn half-edge partner mismatch');
    nodes.set(incoming.u, ['X', 'in', Number(externalIndex[1]), incoming.particle]);
    nodes.set(outgoing.v, ['X', 'out', Number(externalIndex[1]), outgoing.particle]);
  }
  equal(nodes.size, new Set(edges.flatMap((edge) => [edge.u, edge.v])).size, 'node inventory');
  assert([...nodes.values()].every(Boolean), 'unresolved physical node');
  if (expectedKind === 'amplitudes') {
    assert.equal(seams.size, 0, 'vacuum amplitude has a seam');
    assert.equal([...nodes.values()].filter((label) => label[0] === 'X').length, 0);
  } else {
    assert.equal(seams.size, 2, 'expected two labelled initial-state seams');
    assert.equal([...nodes.values()].filter((label) => label[0] === 'X').length, 4);
  }
  return {nodes, edges, ports, seams};
}
function orientation(edge) {
  if (edge.direction === 'none') return {directed: false, source: edge.u, target: edge.v};
  if (edge.direction === 'back') return {directed: true, source: edge.v, target: edge.u};
  assert(['forward', 'both'].includes(edge.direction), `unknown direction ${edge.direction}`);
  return {directed: true, source: edge.u, target: edge.v};
}
function permutations(items) {
  if (items.length < 2) return [items];
  return items.flatMap((item, index) => permutations(items.filter((_, other) => other !== index))
    .map((tail) => [item, ...tail]));
}
function rank(values) {
  const unique = [...new Set(values)].sort();
  return values.map((value) => unique.indexOf(value));
}
function canonical(graph, record) {
  const nodes = [...graph.nodes.keys()];
  let colors = rank(nodes.map((node) => encode(graph.nodes.get(node))));
  for (let round = 0; round < nodes.length; round++) {
    const lookup = new Map(nodes.map((node, index) => [node, colors[index]]));
    const next = rank(nodes.map((node, index) => encode([colors[index], sorted(graph.edges
      .filter((edge) => edge.u === node || edge.v === node)
      .map((edge) => {
        const direction = orientation(edge);
        const flow = direction.directed ? direction.source === node ? 'out' : 'in' : 'both';
        return [edge.particle, flow, edge.pair !== null, lookup.get(edge.u === node ? edge.v : edge.u)];
      }))])));
    if (encode(next) === encode(colors)) break;
    colors = next;
  }
  const classes = [...new Set(colors)].sort((a, b) => a - b)
    .map((color) => nodes.filter((_, index) => colors[index] === color));
  let orders = [[]];
  for (const group of classes) orders = orders.flatMap((order) => permutations(group).map((tail) => [...order, ...tail]));
  assert(orders.length < 100000, 'physical node permutation budget exceeded');
  let key = null;
  let cutKey = null;
  let representative = null;
  let automorphisms = 0;
  for (const order of orders) {
    const map = new Map(order.map((node, index) => [node, index]));
    const rows = graph.edges.map((edge, index) => {
      const direction = orientation(edge);
      let source = map.get(direction.source);
      let target = map.get(direction.target);
      if (!direction.directed && source > target) [source, target] = [target, source];
      return {row: [source, target, edge.particle, direction.directed, edge.pair !== null], index};
    }).sort((a, b) => encode(a.row).localeCompare(encode(b.row), 'en'));
    const candidate = encode([order.map((node) => graph.nodes.get(node)), rows.map((item) => item.row)]);
    if (key !== null && candidate > key) continue;
    if (key === null || candidate < key) {
      key = candidate;
      cutKey = null;
      automorphisms = 0;
    }
    const classesOfEdges = [];
    for (const item of rows) {
      if (classesOfEdges.length && encode(classesOfEdges.at(-1)[0].row) === encode(item.row)) {
        classesOfEdges.at(-1).push(item);
      } else classesOfEdges.push([item]);
    }
    let edgeOrders = [[]];
    for (const group of classesOfEdges) {
      edgeOrders = edgeOrders.flatMap((orderOfEdges) => permutations(group).map((tail) => [...orderOfEdges, ...tail]));
    }
    automorphisms += edgeOrders.length;
    for (const edgeOrder of edgeOrders) {
      const edgeMap = new Map(edgeOrder.map((item, index) => [item.index, index]));
      const selected = (hedges) => sorted(hedges.map((hedge) => {
        const owner = graph.ports.get(hedge);
        assert(owner, 'half-edge has no physical owner');
        return [edgeMap.get(owner.edge), map.get(owner.node)];
      }));
      const cut = (selection) => Object.fromEntries(['left', 'right', 'left_cut', 'right_cut']
        .map((field) => [field, selected(selection[field])]));
      const seams = sorted([...graph.seams.values()].map((group) =>
        sorted(group.map((edge) => edgeMap.get(graph.edges.indexOf(edge))))));
      const candidateCuts = encode({
        initial: {left: selected(record.initial.left), right: selected(record.initial.right)},
        physical: sorted(record.cuts.map(cut)),
        thresholds: sorted(record.thresholds.map(cut)),
        seams,
      });
      if (cutKey === null || candidateCuts < cutKey) {
        cutKey = candidateCuts;
        representative = {
          nodes: Object.fromEntries(map),
          edges: Object.fromEntries(edgeMap),
        };
      }
    }
  }
  return {key, cutKey, representative, automorphisms, nodeSearchOrders: orders.length};
}
function cutChecks(graph, record, kind) {
  const hedges = [...graph.ports.keys()].sort((a, b) => a - b);
  let properSubgraphThresholds = 0;
  equal(hedges, Array.from({length: hedges.length}, (_, index) => index), 'noncontiguous physical half-edges');
  assert(Array.isArray(record.initial.left) && Array.isArray(record.initial.right),
    'initial cut sides must be arrays');
  const leftInitial = new Set(record.initial.left);
  const rightInitial = new Set(record.initial.right);
  equal(leftInitial.size, record.initial.left.length, 'duplicate initial left half-edge');
  equal(rightInitial.size, record.initial.right.length, 'duplicate initial right half-edge');
  assert([...leftInitial, ...rightInitial].every((hedge) => graph.ports.has(hedge)),
    'initial cut uses unknown physical half-edge');
  const expectedInitial = {left: [], right: []};
  for (const group of graph.seams.values()) {
    const incoming = group.find((edge) => edge.pu === null);
    const outgoing = group.find((edge) => edge.pv === null);
    expectedInitial.left.push(incoming.pv);
    expectedInitial.right.push(outgoing.pu);
  }
  equal(sorted([...leftInitial]), sorted(expectedInitial.left), 'initial left seam incidence');
  equal(sorted([...rightInitial]), sorted(expectedInitial.right), 'initial right seam incidence');
  for (const [type, selections] of [['physical', record.cuts], ['threshold', record.thresholds]]) for (const selection of selections) {
    const left = new Set(selection.left);
    const right = new Set(selection.right);
    const leftCut = new Set(selection.left_cut);
    const rightCut = new Set(selection.right_cut);
    equal(left.size, selection.left.length);
    equal(right.size, selection.right.length);
    equal(leftCut.size, selection.left_cut.length, 'duplicate left cut half-edge');
    equal(rightCut.size, selection.right_cut.length, 'duplicate right cut half-edge');
    assert([...left, ...right].every((hedge) => graph.ports.has(hedge)), 'selection uses unknown half-edge');
    equal(new Set([...left, ...right]).size, left.size + right.size, 'selection sides overlap');
    assert([...leftCut].every((hedge) => left.has(hedge)), 'left cut outside selected left subgraph');
    assert([...rightCut].every((hedge) => right.has(hedge)), 'right cut outside selected right subgraph');
    if (type === 'physical') {
      equal(sorted([...left, ...right]), sorted(hedges), 'physical cut is not a full half-edge partition');
      for (const hedge of leftInitial) assert(left.has(hedge), 'incoming seam moved to right');
      for (const hedge of rightInitial) assert(right.has(hedge), 'outgoing seam moved to left');
    } else if (left.size + right.size < hedges.length) properSubgraphThresholds++;
    const vertexSide = new Map();
    for (const hedge of [...left, ...right]) {
      const owner = graph.ports.get(hedge);
      const onLeft = left.has(hedge);
      if (vertexSide.has(owner.node)) equal(vertexSide.get(owner.node), onLeft, 'cut splits physical vertex');
      else vertexSide.set(owner.node, onLeft);
    }
    const crossingLeft = [];
    const crossingRight = [];
    for (const edge of graph.edges) {
      if (edge.pu === null || edge.pv === null) continue;
      if (!left.has(edge.pu) && !right.has(edge.pu)) continue;
      if (!left.has(edge.pv) && !right.has(edge.pv)) continue;
      if (left.has(edge.pu) === left.has(edge.pv)) continue;
      crossingLeft.push(left.has(edge.pu) ? edge.pu : edge.pv);
      crossingRight.push(right.has(edge.pu) ? edge.pu : edge.pv);
    }
    equal(sorted([...leftCut]), sorted(crossingLeft), `${type} left crossing incidence`);
    equal(sorted([...rightCut]), sorted(crossingRight), `${type} right crossing incidence`);
  }
  if (kind === 'amplitudes') {
    equal(record.cuts, []);
    equal(record.thresholds, []);
  } else {
    assert(record.cuts.length > 0, 'cross section has no selected physical cut');
  }
  return properSubgraphThresholds;
}
function componentData(graph) {
  const anticommuting = (name) => {
    const spin = particle.get(name)?.spin;
    return spin === 2 || spin === -1;
  };
  const adjacency = new Map();
  const incident = new Map();
  for (const edge of graph.edges) {
    if (!anticommuting(edge.particle)) continue;
    assert(orientation(edge).directed, 'anticommuting edge has no arrow');
    for (const [node, next] of [[edge.u, edge.v], [edge.v, edge.u]]) {
      if (!adjacency.has(node)) adjacency.set(node, []);
      adjacency.get(node).push(next);
      if (!incident.has(node)) incident.set(node, []);
      incident.get(node).push(edge);
    }
  }
  const visited = new Set();
  const components = [];
  for (const start of adjacency.keys()) {
    if (visited.has(start)) continue;
    const pending = [start];
    const nodes = new Set();
    while (pending.length) {
      const node = pending.pop();
      if (visited.has(node)) continue;
      visited.add(node);
      nodes.add(node);
      pending.push(...adjacency.get(node));
    }
    const external = [...nodes].filter((node) => graph.nodes.get(node)[0] === 'X');
    const edges = [...new Set([...nodes].flatMap((node) => incident.get(node)))];
    const ghosts = edges.map((edge) => particle.get(edge.particle).spin === -1);
    if (external.length === 0) {
      assert(ghosts.every((value) => value === ghosts[0]), 'mixed ghost and fermion closed component');
    } else {
      equal(external.length, 2, 'open external Grassmann chain endpoints');
      assert(ghosts.every((value) => !value), 'external ghost chain');
    }
    let sources = 0;
    let sinks = 0;
    for (const node of nodes) {
      const flows = incident.get(node).map((edge) => orientation(edge).source === node ? 1 : -1);
      const label = graph.nodes.get(node);
      if (label[0] === 'V') {
        equal(sorted(flows), [-1, 1], 'internal Grassmann arrow flow');
      } else {
        equal(flows.length, 1, 'external Grassmann endpoint degree');
        const externalParticle = particle.get(label[3]);
        assert.equal(externalParticle.spin, 2, 'open external Grassmann endpoint is not a fermion');
        const expectedSource = (label[1] === 'in') !== (externalParticle.pdg_code < 0);
        equal(flows[0], expectedSource ? 1 : -1, 'external particle/role arrow flow');
        if (flows[0] === 1) sources++; else sinks++;
      }
    }
    if (external.length) equal([sources, sinks], [1, 1], 'open Grassmann source/sink count');
    components.push({external, ghost: ghosts[0] ?? false});
  }
  const closed = components.filter((component) => component.external.length === 0);
  const quarkLoops = closed.filter((component) => !component.ghost).length;
  const ghostLoops = closed.filter((component) => component.ghost).length;
  const open = components.filter((component) => component.external.length !== 0);
  const externalAdjacency = new Map();
  const connect = (left, right) => {
    if (!externalAdjacency.has(left)) externalAdjacency.set(left, []);
    if (!externalAdjacency.has(right)) externalAdjacency.set(right, []);
    externalAdjacency.get(left).push(right);
    externalAdjacency.get(right).push(left);
  };
  for (const component of open) connect(...component.external);
  for (const pair of graph.seams.values()) {
    const first = pair.find((edge) => edge.pu === null);
    const second = pair.find((edge) => edge.pv === null);
    if (anticommuting(first.particle)) connect(first.u, second.v);
  }
  for (const neighbors of externalAdjacency.values()) equal(neighbors.length, 2, 'external closure degree');
  const externalVisited = new Set();
  let externalCycles = 0;
  for (const start of externalAdjacency.keys()) {
    if (externalVisited.has(start)) continue;
    externalCycles++;
    const pending = [start];
    while (pending.length) {
      const node = pending.pop();
      if (externalVisited.has(node)) continue;
      externalVisited.add(node);
      pending.push(...externalAdjacency.get(node));
    }
  }
  const externalOrderingSign = (open.length + externalCycles) % 2 ? -1 : 1;
  const incomingAntifermions = [...graph.nodes.values()].filter((label) => label[0] === 'X' && label[1] === 'in'
    && particle.get(label[3]).spin === 2 && particle.get(label[3]).pdg_code < 0).length;
  return {quarkLoops, ghostLoops, openChains: open.length, externalCycles,
    externalOrderingSign, antiSpinSumSign: incomingAntifermions % 2 ? -1 : 1,
    incomingAntifermions};
}
function rational(value) {
  assert(/^-?\d+(\/[1-9]\d*)?$/.test(value), `nonrational or zero-denominator weight ${value}`);
  const [numerator, denominator = '1'] = value.split('/');
  return reduce({n: BigInt(numerator), d: BigInt(denominator)});
}
function gcd(first, second) {
  let a = first < 0n ? -first : first;
  let b = second;
  while (b) [a, b] = [b, a % b];
  return a;
}
function reduce(value) {
  assert(value.d !== 0n, 'zero denominator');
  const sign = value.d < 0n ? -1n : 1n;
  const common = gcd(value.n, value.d * sign);
  return {n: sign * value.n / common, d: sign * value.d / common};
}
function multiply(first, second) {
  return reduce({n: first.n * second.n, d: first.d * second.d});
}
function divide(first, second) {
  assert(second.n !== 0n, 'division by zero');
  return reduce({n: first.n * second.d, d: first.d * second.n});
}
function add(first, second) {
  return reduce({n: first.n * second.d + second.n * first.d, d: first.d * second.d});
}
function printable(value) {
  assert(Number.isSafeInteger(Number(value.n)) && Number.isSafeInteger(Number(value.d)),
    'JSON integer exceeds exact safe range');
  return {numerator: Number(value.n), denominator: Number(value.d)};
}
function parseFactorExpression(source) {
  // Grammar: product -> unary (*|/ unary)*; unary -> (- unary) | power;
  // power -> primary [^ signed-integer]; primary -> integer | head | (product).
  // Power binds before unary signs, as in pinned Symbolica. Additive expressions,
  // chained powers and non-rational powers are outside this bounded scalar grammar.
  assert(source.length <= 4096, 'scalar factor exceeds bounded grammar size');
  let cursor = 0;
  const heads = [];
  const skip = () => { while (/\s/.test(source[cursor] ?? '')) cursor++; };
  const peek = () => { skip(); return source[cursor]; };
  const integer = (signed = false, allowPlus = false) => {
    skip();
    const pattern = signed ? allowPlus ? '^[+-]?\\d+' : '^-?\\d+' : '^\\d+';
    const match = new RegExp(pattern).exec(source.slice(cursor));
    assert(match, `integer expected at ${cursor}`);
    cursor += match[0].length;
    return BigInt(match[0]);
  };
  const power = (base, exponent) => {
    assert(exponent >= -1n && exponent <= 1n, 'unsupported scalar power');
    if (exponent === -1n) return divide({n: 1n, d: 1n}, base);
    if (exponent === 0n) {
      assert(base.n !== 0n, 'indeterminate zero to zero power');
      return {n: 1n, d: 1n};
    }
    return base;
  };
  function primary() {
    let value;
    if (peek() === '(') {
      cursor++;
      value = product();
      assert(peek() === ')', 'unclosed factor parenthesis');
      cursor++;
    } else {
      skip();
      const head = /^(?:feynkit_generator_factor::)?([A-Za-z]+)\(/.exec(source.slice(cursor));
      if (head) {
        cursor += head[0].length;
        const argument = integer(true);
        assert(peek() === ')', 'factor head arity or argument');
        cursor++;
        assert(['AutG', 'InternalFermionLoopSign', 'ExternalFermionOrderingSign',
          'AntiFermionSpinSumSign'].includes(head[1]), `unsupported factor head ${head[1]}`);
        heads.push([head[1], argument]);
        value = {n: argument, d: 1n};
      } else {
        assert(/\d/.test(peek() ?? ''), `unsupported factor atom at ${cursor}`);
        value = {n: integer(), d: 1n};
      }
    }
    return value;
  }
  function powered() {
    let value = primary();
    if (peek() === '^') {
      cursor++;
      const parenthesized = peek() === '(';
      if (parenthesized) cursor++;
      const exponent = integer(true, true);
      if (parenthesized) {
        assert(peek() === ')', 'unclosed scalar exponent');
        cursor++;
      }
      value = power(value, exponent);
      assert(peek() !== '^', 'chained scalar powers are unsupported');
    }
    return value;
  }
  function unary() {
    assert(peek() !== '+', 'unary plus is outside the bounded factor grammar');
    if (peek() === '-') {
      cursor++;
      const value = unary();
      return {n: -value.n, d: value.d};
    }
    return powered();
  }
  function product() {
    let value = unary();
    while (peek() === '*' || peek() === '/') {
      const operator = source[cursor++];
      const next = unary();
      value = operator === '*' ? multiply(value, next) : divide(value, next);
    }
    return value;
  }
  const value = product();
  assert(cursor === source.length, `unsupported or unconsumed factor syntax at ${cursor}`);
  return {value, heads};
}
function verifyFactor(record, physical, automorphisms, kind) {
  const parsed = parseFactorExpression(record.factor);
  const internal = (physical.quarkLoops + physical.ghostLoops) % 2 ? -1 : 1;
  const expectedHeads = [['AutG', BigInt(automorphisms)],
    ['ExternalFermionOrderingSign', BigInt(physical.externalOrderingSign)]];
  if (kind !== 'amplitudes') expectedHeads.push(['AntiFermionSpinSumSign', BigInt(physical.antiSpinSumSign)]);
  if (internal < 0) expectedHeads.push(['InternalFermionLoopSign', -1n]);
  const ordered = (heads) => heads.map(([name, value]) => [name, String(value)])
    .sort((left, right) => left[0].localeCompare(right[0], 'en'));
  equal(ordered(parsed.heads), ordered(expectedHeads), 'complete factor head occurrences and arguments');
  const evaluated = rational(record.weight);
  const signed = BigInt(internal * physical.externalOrderingSign *
    (kind === 'amplitudes' ? 1 : physical.antiSpinSumSign));
  const expectedWeight = reduce({n: signed, d: BigInt(automorphisms)});
  equal(parsed.value, expectedWeight, 'complete raw scalar factor versus graph-derived weight');
  equal(evaluated, expectedWeight, 'recorded evaluated scalar versus graph-derived weight');
  const ghostless = reduce({n: BigInt((physical.quarkLoops % 2 ? -1 : 1) * physical.externalOrderingSign
    * (kind === 'amplitudes' ? 1 : physical.antiSpinSumSign)), d: BigInt(automorphisms)});
  return {evaluated, ghostless};
}
function dotScalarField(dot, field) {
  const matches = [...dot.matchAll(new RegExp(`^\\s*${field} = ("(?:[^"\\\\]|\\\\.)*");$`, 'gm'))];
  equal(matches.length, 1, `DOT ${field} occurrence count`);
  return JSON.parse(matches[0][1]);
}
function readCase(root, configuration, label, fileSystem = fs) {
  const family = path.join(root, configuration.family);
  const witness = path.join(family, `${label}-graph-witness.tsv`);
  const source = fileSystem.readFileSync(witness, 'utf8');
  const lines = source.trimEnd().split('\n');
  const header = lines.shift();
  const definitionPath = path.join(family, `${label}-process-definition.json`);
  const definition = JSON.parse(fileSystem.readFileSync(definitionPath));
  equal(definition.grouping, 'none');
  const graphDirectory = path.join(family, `${label}-dot/processes`, configuration.kind, configuration.folder, 'default');
  const rows = [];
  let sum = {n: 0n, d: 1n};
  let ghostless = {n: 0n, d: 1n};
  for (const line of lines) {
    const fields = line.split('\t');
    equal(fields.length, 8, 'graph witness column count');
    const [index, name, factor, weight, group, initial, cuts, thresholds] = fields;
    const record = {index: Number(index), name, factor, weight, group,
      initial: JSON.parse(initial), cuts: JSON.parse(cuts), thresholds: JSON.parse(thresholds)};
    const dotPath = path.join(graphDirectory, `${name}.dot`);
    const dot = fileSystem.readFileSync(dotPath, 'utf8');
    const graph = parseDot(dot, configuration.kind);
    const properSubgraphThresholds = cutChecks(graph, record, configuration.kind);
    const neutral = canonical(graph, record);
    const components = componentData(graph);
    const factors = verifyFactor(record, components, neutral.automorphisms, configuration.kind);
    const dotFactor = dotScalarField(dot, 'overall_factor');
    const dotWeight = dotScalarField(dot, 'overall_factor_evaluated');
    const dotFactors = verifyFactor({...record, factor: dotFactor, weight: dotWeight}, components,
      neutral.automorphisms, configuration.kind);
    equal(dotFactors.evaluated, factors.evaluated, 'DOT and TSV parsed raw scalar values differ');
    equal(dotWeight, weight, 'DOT evaluated factor differs from graph witness');
    sum = add(sum, factors.evaluated);
    ghostless = add(ghostless, factors.ghostless);
    rows.push({name, index: record.index, key: neutral.key, keySha256: sha(neutral.key),
      cutKey: neutral.cutKey, cutSha256: sha(neutral.cutKey), map: neutral.representative,
      factor, weight, dotFactor, dotWeight, aut: neutral.automorphisms, ...components,
      cuts: record.cuts.length, thresholds: record.thresholds.length, properSubgraphThresholds,
      dotPath, dotSha256: sha(dot)});
  }
  equal(rows.length, Number(header.split(' | ')[0]), 'header graph count');
  const headerWeight = rational(header.split(' | ')[1].split(' = ')[0]);
  equal(sum.n * headerWeight.d, headerWeight.n * sum.d, 'header rational sum');
  const names = fileSystem.readdirSync(graphDirectory).filter((name) => name.endsWith('.dot'));
  equal(names.length, rows.length, 'DOT count does not match graph rows');
  equal(new Set(rows.map((row) => row.name)).size, rows.length, 'duplicate graph name');
  return {header, definition, definitionPath, witness, witnessSha256: sha(source),
    rows, count: rows.length, sum: printable(sum), ghostless: printable(ghostless)};
}
function histogram(rows, projection) {
  const result = {};
  for (const row of rows) {
    const key = String(projection(row));
    result[key] = (result[key] ?? 0) + 1;
  }
  return result;
}
function matched(left, right, id, label) {
  equal(left.count, right.count, `${id}/${label} count mismatch`);
  equal(left.definition, right.definition, `${id}/${label} resolved inputs differ`);
  const leftKeys = new Map(left.rows.map((row) => [row.key, row]));
  const rightKeys = new Map(right.rows.map((row) => [row.key, row]));
  equal(leftKeys.size, left.rows.length, `${id}/${label} duplicate historical physical key`);
  equal(rightKeys.size, right.rows.length, `${id}/${label} duplicate current physical key`);
  equal([...leftKeys.keys()].sort(), [...rightKeys.keys()].sort(), `${id}/${label} physical membership differs`);
  const joined = [];
  for (const [key, old] of leftKeys) {
    const current = rightKeys.get(key);
    equal(old.cutKey, current.cutKey, `${id}/${label} cut incidence differs`);
    equal(old.factor, current.factor, `${id}/${label} raw factor differs`);
    equal(old.weight, current.weight, `${id}/${label} evaluated weight differs`);
    equal(old.aut, current.aut, `${id}/${label} automorphism differs`);
    equal(old.properSubgraphThresholds, current.properSubgraphThresholds,
      `${id}/${label} proper-subgraph threshold count differs`);
    equal([old.quarkLoops, old.ghostLoops, old.externalOrderingSign, old.antiSpinSumSign],
      [current.quarkLoops, current.ghostLoops, current.externalOrderingSign, current.antiSpinSumSign],
      `${id}/${label} physical sign data differs`);
    joined.push({keySha256: sha(key), cutSha256: sha(old.cutKey),
      historicalName: old.name, currentName: current.name,
      historicalMap: old.map, currentMap: current.map,
      factor: old.factor, weight: old.weight, aut: old.aut,
      quarkLoops: old.quarkLoops, ghostLoops: old.ghostLoops,
      externalOrderingSign: old.externalOrderingSign, antiSpinSumSign: old.antiSpinSumSign,
      selectedCuts: old.cuts, thresholds: old.thresholds});
  }
  joined.sort((a, b) => a.keySha256.localeCompare(b.keySha256));
  return joined;
}
function mutationControls(root) {
  const directory = path.join(root, 'epem_bbh');
  const line = fs.readFileSync(path.join(directory, 'bbh_qcd2-graph-witness.tsv'), 'utf8').split('\n')[1];
  const [index, name, factor, weight, group, initial, cuts, thresholds] = line.split('\t');
  const record = {index: Number(index), name, factor, weight, group,
    initial: JSON.parse(initial), cuts: JSON.parse(cuts), thresholds: JSON.parse(thresholds)};
  const dot = fs.readFileSync(path.join(directory, 'bbh_qcd2-dot/processes/cross_sections/epem_bbxh/default', `${name}.dot`), 'utf8');
  const graph = parseDot(dot, 'cross_sections');
  const original = canonical(graph, record);
  const physical = componentData(graph);
  verifyFactor(record, physical, original.automorphisms, 'cross_sections');
  const rejected = (description, operation) => {
    assert.throws(operation, undefined, `semantic mutation escaped: ${description}`);
    return {description, outcome: 'rejected'};
  };
  const results = [];
  results.push(rejected('particle label with graph count fixed', () => {
    const changed = structuredClone(graph);
    changed.edges.find((edge) => edge.pu !== null && edge.pv !== null).particle = 'g';
    equal(canonical(changed, record).key, original.key);
  }));
  results.push(rejected('selected physical cut incidence', () => {
    const changed = structuredClone(record);
    changed.cuts[0].left_cut.pop();
    cutChecks(graph, changed, 'cross_sections');
  }));
  results.push(rejected('threshold incidence while graph and total are fixed', () => {
    const changed = structuredClone(record);
    changed.thresholds[0].right_cut.pop();
    equal(canonical(graph, changed).cutKey, original.cutKey);
  }));
  results.push(rejected('fermion factor with graph and evaluated weight fixed', () => {
    const changed = {...record, factor: record.factor.replace('InternalFermionLoopSign(-1)*', '')};
    verifyFactor(changed, physical, original.automorphisms, 'cross_sections');
  }));
  results.push(rejected('per-graph rational weight under unchanged aggregate count', () => {
    const changed = {...record, weight: '1'};
    verifyFactor(changed, physical, original.automorphisms, 'cross_sections');
  }));
  const second = fs.readFileSync(path.join(directory, 'bbh_qcd2-graph-witness.tsv'), 'utf8').split('\n')[2].split('\t');
  equal(second[3], record.weight, 'same-weight membership control needs equal scalar weights');
  const secondRecord = {initial: JSON.parse(second[5]), cuts: JSON.parse(second[6]), thresholds: JSON.parse(second[7])};
  const secondGraph = parseDot(fs.readFileSync(path.join(directory,
    'bbh_qcd2-dot/processes/cross_sections/epem_bbxh/default', `${second[1]}.dot`), 'utf8'), 'cross_sections');
  const secondKey = canonical(secondGraph, secondRecord).key;
  assert.notEqual(original.key, secondKey, 'membership control requires distinct physical keys');
  results.push(rejected('replace a physical member with a same-weight duplicate, preserving graph count and sum', () => {
    const changedKeys = [original.key, original.key];
    equal(new Set(changedKeys).size, 2, 'duplicate physical key under same total');
  }));
  const renamed = structuredClone(graph);
  const nodeNames = new Map([...graph.nodes.keys()].map((node, i) => [node, `renamed${i}`]));
  const hedgeNames = new Map([...graph.ports.keys()].map((hedge, i) => [hedge, 10000 + i]));
  renamed.nodes = new Map([...graph.nodes].map(([node, label]) => [nodeNames.get(node), label]));
  renamed.edges = graph.edges.map((edge) => ({...edge, u: nodeNames.get(edge.u), v: nodeNames.get(edge.v),
    pu: edge.pu === null ? null : hedgeNames.get(edge.pu),
    pv: edge.pv === null ? null : hedgeNames.get(edge.pv)}));
  renamed.ports = new Map([...graph.ports].map(([hedge, owner]) =>
    [hedgeNames.get(hedge), {...owner, node: nodeNames.get(owner.node)}]));
  renamed.seams = new Map([...graph.seams].map(([key, edges]) =>
    [key, edges.map((edge) => renamed.edges[graph.edges.indexOf(edge)])]));
  const remap = (selection) => Object.fromEntries(Object.entries(selection)
    .map(([field, values]) => [field, values.map((hedge) => hedgeNames.get(hedge))]));
  const relabelled = {...record, initial: remap(record.initial),
    cuts: record.cuts.map(remap).reverse(), thresholds: record.thresholds.map(remap).reverse()};
  const neutral = canonical(renamed, relabelled);
  equal([neutral.key, neutral.cutKey, neutral.automorphisms],
    [original.key, original.cutKey, original.automorphisms], 'neutral relabel should preserve physical join');
  results.push({description: 'node and half-edge relabel plus cut-list permutation', outcome: 'accepted'});
  return results;
}

function main() {
if (!currentRoot || !historicalRoot || !fullOutput || !summaryOutput) {
  throw new Error('usage: node analyze-graphs.mjs <current-capture> <source43-capture> <full.json> <summary.json>');
}
const full = {schema: 1, currentRoot, historicalRoot, cases: {}};
const summary = {schema: 1, status: 'independent physical graph, cut and scalar-weight derivation; review pending',
  method: 'complete rule-coloured directed multigraph key with edge-occurrence permutations, labelled sewn external states and exact half-edge cut maps; closed quark/ghost components and external closure signs from the live model; BigInt rational factor evaluation',
  cases: {}};
for (const configuration of cases) {
  const labels = [configuration.label, `${configuration.label}-no-zero-filter`];
  const datasets = {};
  for (const label of labels) {
    const old = readCase(historicalRoot, configuration, label);
    const current = readCase(currentRoot, configuration, label);
    const joins = matched(old, current, configuration.id, label);
    datasets[label] = {old, current, joins};
    full.cases[`${configuration.id}/${label}`] = {
      historical: old.rows, current: current.rows, joins,
      historicalDefinition: old.definition, currentDefinition: current.definition,
      historicalWitnessSha256: old.witnessSha256, currentWitnessSha256: current.witnessSha256,
    };
    summary.cases[`${configuration.id}/${label}`] = {
      graphsPerRoute: old.count, uniqueNeutralKeys: joins.length,
      sum: old.sum, ghostlessSum: old.ghostless,
      historicalWitnessSha256: old.witnessSha256,
      currentWitnessSha256: current.witnessSha256,
      definitionsByteEqual: fs.readFileSync(old.definitionPath).equals(fs.readFileSync(current.definitionPath)),
      factors: histogram(old.rows, (row) => row.factor),
      ghostLoops: histogram(old.rows, (row) => row.ghostLoops),
      automorphisms: histogram(old.rows, (row) => row.aut),
      selectedCuts: old.rows.reduce((count, row) => count + row.cuts, 0),
      thresholdCandidates: old.rows.reduce((count, row) => count + row.thresholds, 0),
      properSubgraphThresholds: old.rows.reduce((count, row) => count + row.properSubgraphThresholds, 0),
    };
    console.log(configuration.id, label, old.count, `${old.sum.numerator}/${old.sum.denominator}`,
      `ghostless ${old.ghostless.numerator}/${old.ghostless.denominator}`);
  }
  const filtered = datasets[configuration.label].old;
  const unfiltered = datasets[`${configuration.label}-no-zero-filter`].old;
  const kept = new Set(filtered.rows.map((row) => row.key));
  for (const row of filtered.rows) {
    const retained = unfiltered.rows.find((candidate) => candidate.key === row.key);
    assert(retained, `${configuration.id} filtered graph absent without filter`);
    equal([row.factor, row.weight, row.cutKey], [retained.factor, retained.weight, retained.cutKey],
      `${configuration.id} retained graph changed under zero-filter toggle`);
  }
  const removed = unfiltered.rows.filter((row) => !kept.has(row.key));
  let removedWeight = {n: 0n, d: 1n};
  for (const row of removed) removedWeight = add(removedWeight, rational(row.weight));
  summary.cases[`${configuration.id}/${configuration.label}`].zeroFilterRemoved = {count: removed.length,
    weight: printable(removedWeight), keys: removed.map((row) => sha(row.key)).sort()};
}
summary.semanticMutationControls = mutationControls(currentRoot);
fs.writeFileSync(fullOutput, `${JSON.stringify(full)}\n`);
summary.fullWitnessPath = fullOutput;
summary.fullWitnessSha256 = sha(fs.readFileSync(fullOutput));
fs.writeFileSync(summaryOutput, `${JSON.stringify(summary)}\n`);
}
export {cases, parseDot, canonical, cutChecks, componentData, parseFactorExpression,
  verifyFactor, readCase, matched, mutationControls};
if (process.argv[1] && path.resolve(process.argv[1]) === new URL(import.meta.url).pathname) main();
