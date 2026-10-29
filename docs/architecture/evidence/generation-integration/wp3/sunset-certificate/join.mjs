import fs from 'node:fs';
import {liveRule, liveModel, liveGlobals, liveDots} from './live-rule.mjs';
import path from 'node:path';
import crypto from 'node:crypto';
import {complete, parse, shape, same as sameFactor} from './factor.mjs';
import {checkAggregate, tensorKey, nativeCliExpression} from './tensor-product.mjs';

const [historical, current, output] = process.argv.slice(2);
if (!historical || !current || !output) {
  throw new Error('usage: node join.mjs <historical-capture> <current-capture> <output.json>');
}
const digest = (value) => crypto.createHash('sha256').update(value).digest('hex');
const scalarContract = (actual, expected, where) => {
  const lhs = complete(actual);
  const rhs = complete(expected);
  const ordered = (items) => items.map((item) => JSON.stringify(item)).sort();
  if (!sameFactor(lhs.value, rhs.value) ||
      JSON.stringify(ordered(lhs.heads)) !== JSON.stringify(ordered(rhs.heads)) ||
      JSON.stringify(ordered(lhs.groups)) !== JSON.stringify(ordered(rhs.groups))) {
    throw new Error(`${where}: scalar head/argument/occurrence contract differs`);
  }
};
const labels = ['all-raw', 'all-grouped', 'nonfactorizable-raw', 'nonfactorizable-grouped'];
const gcd = (a, b) => b === 0n ? a : gcd(b, a % b);
const rational = (source) => {
  const [n, d = '1'] = source.split('/');
  const numerator = BigInt(n);
  const denominator = BigInt(d);
  if (denominator === 0n) throw new Error('zero rational denominator');
  const factor = gcd(numerator < 0n ? -numerator : numerator, denominator);
  return [numerator / factor, denominator / factor];
};
const add = ([a, b], [c, d]) => rational(`${a * d + c * b}/${b * d}`);
const multiply = ([a, b], [c, d]) => rational(`${a * c}/${b * d}`);
const fractionString = ([n, d]) => d === 1n ? String(n) : `${n}/${d}`;
const sum = (values) => values.reduce(add, [0n, 1n]);

function physical(vertices, edges) {
  if (!vertices.length || !edges.length) throw new Error('incomplete vacuum graph');
  const permutations = (items) => items.length < 2 ? [items] : items.flatMap((item, index) =>
    permutations(items.filter((_, offset) => offset !== index)).map((tail) => [item, ...tail]));
  if (vertices.length > 7) throw new Error('unbounded canonicalization');
  let best;
  for (const order of permutations(vertices)) {
    const index = Object.fromEntries(order.map((vertex, position) => [vertex.id, position]));
    const physicalEdges = edges.map((edge) => {
      let from = index[edge.from];
      let to = index[edge.to];
      if (edge.dir === 'back') [from, to] = [to, from];
      if (edge.dir === 'none' && from > to) [from, to] = [to, from];
      return [from, to, edge.particle, edge.dir === 'none' ? 'undirected' : 'directed'];
    }).sort((a, b) => JSON.stringify(a).localeCompare(JSON.stringify(b)));
    const candidate = JSON.stringify({vertices: order.map((vertex) => vertex.int_id),
      edges: physicalEdges});
    if (!best || candidate < best) best = candidate;
  }
  return {key: digest(best), canonical: JSON.parse(best)};
}

function graph(file, parsed) {
  const source = fs.readFileSync(file, 'utf8');
  const vertices = parsed.vertices.map(({id, owner, attributes, numInterface}) => ({
    id, owner, ...attributes, numInterface,
  }));
  const edges = parsed.edges.map(({id, edgeId, from, fromOwner, fromPort, fromPortLabel,
    to, toOwner, toPort, toPortLabel, dir, attributes, numInterface}) => ({
    ...attributes, id, edgeId, from, fromOwner, fromPort, fromPortLabel,
    to, toOwner, toPort, toPortLabel, numInterface,
    dir: dir === 'forward' ? undefined : dir,
  }));
  if (vertices.some((vertex, index) => !Number.isSafeInteger(vertex.owner) || vertex.owner !== index) ||
      edges.some((edge, index) => !Number.isSafeInteger(edge.edgeId) || edge.edgeId !== index ||
        !Number.isSafeInteger(edge.fromOwner) || !Number.isSafeInteger(edge.toOwner) ||
        !Number.isSafeInteger(edge.fromPort) || !Number.isSafeInteger(edge.toPort))) {
    throw new Error(`${file}: parsed DOT owner IDs are not a bounded exact inventory`);
  }
  const globalField = (name) => parsed.global[name];
  return {...physical(vertices, edges), vertices, edges,
    globalNumerator: globalField('num'),
    numeratorPrefactorInterface: parsed.globalInterfaces.numeratorPrefactor ?? undefined,
    projector: globalField('projector'),
    projectorInterface: parsed.globalInterfaces.projector ?? undefined,
    overallFactor: globalField('overall_factor'),
    dotSha256: digest(source)};
}

function nativePhysical(diagram, model) {
  const vertices = diagram.vertices.map((vertex, index) => ({
    id: `v${index}`, int_id: model.vertex_rules[vertex.interaction]?.name,
  }));
  const edges = diagram.edges.map(([endpoints, edge]) => ({
    from: `v${endpoints.source}`, to: `v${endpoints.target}`,
    particle: model.particles[edge.particle]?.name,
    dir: edge.directed ? edge.orientation === 'Reversed' ? 'back' : undefined : 'none',
  }));
  if (vertices.some((vertex) => !vertex.int_id) || edges.some((edge) => !edge.particle)) {
    throw new Error(`native model identity missing for ${diagram.name}`);
  }
  return physical(vertices, edges);
}

function bindCliPayload(diagram, record, model, where) {
  if (diagram.vertices.length !== record.vertices.length || diagram.edges.length !== record.edges.length ||
      new Set(record.vertices.map((vertex) => vertex.id)).size !== diagram.vertices.length ||
      diagram.name !== record.name || nativePhysical(diagram, model).key !== record.key) {
    throw new Error(`${where}: native and CLI graph inventory differ`);
  }
  for (const [index, vertex] of diagram.vertices.entries()) {
    const cli = record.vertices[index];
    if (cli.owner !== index || cli.int_id !== model.vertex_rules[vertex.interaction].name ||
        tensorKey(cli.num) !==
          tensorKey(nativeCliExpression(vertex.numerator)) ||
        JSON.stringify(cli.numInterface) !==
          JSON.stringify(vertex.num_interface ?? null)) {
      throw new Error(`${where}/vertex/${index}: CLI payload differs from native owner`);
    }
  }
  for (const [index, [endpoints, edge]] of diagram.edges.entries()) {
    const cli = record.edges[index];
    const sourcePort = diagram.half_edge_order.findIndex((halfEdge) =>
      halfEdge.edge === index && halfEdge.endpoint === 'source');
    const targetPort = diagram.half_edge_order.findIndex((halfEdge) =>
      halfEdge.edge === index && halfEdge.endpoint === 'target');
    if (cli.edgeId !== index || cli.fromOwner !== endpoints.source ||
        cli.toOwner !== endpoints.target ||
        cli.from !== record.vertices[endpoints.source].id ||
        cli.to !== record.vertices[endpoints.target].id || cli.fromPort !== sourcePort ||
        cli.toPort !== targetPort || cli.particle !== model.particles[edge.particle].name ||
        (cli.dir ?? undefined) !== (edge.directed ?
          edge.orientation === 'Reversed' ? 'back' : undefined : 'none') ||
        tensorKey(cli.num) !==
          tensorKey(nativeCliExpression(edge.numerator)) ||
        JSON.stringify(cli.numInterface) !==
          JSON.stringify(edge.num_interface ?? null)) {
      throw new Error(`${where}/edge/${index}: CLI payload differs from native owner`);
    }
  }
}

function bindLightLocalFactors(record, model) {
  const quarks = record.edges.filter((edge) => ['d', 'u'].includes(edge.particle));
  const gluons = record.edges.filter((edge) => edge.particle === 'g');
  if (record.vertices.length !== 2 || quarks.length !== 2 || gluons.length !== 1 ||
      quarks.some((edge) => edge.dir === 'back' || edge.dir === 'none')) {
    throw new Error(`${record.name}: light sunset incidence differs`);
  }
  const quark = quarks[0].particle;
  if (quarks.some((edge) => edge.particle !== quark)) throw new Error('mixed light quark loop');
  const ruleName = quark === 'd' ? 'V_74' : 'V_135';
  const rule = model.vertex_rules.find((vertex) => vertex.name === ruleName);
  const propagator = model.propagators.find((entry) => entry.particle === quark);
  const gluonPropagator = model.propagators.find((entry) => entry.particle === 'g');
  if (!rule || JSON.stringify(rule.couplings) !== '[["GC_11"]]' ||
      JSON.stringify(rule.lorentz_structures) !== '["FFV1"]' ||
      JSON.stringify(rule.color_structures) !== '["T(3,2,1)"]' ||
      !propagator || propagator.numerator !== '(𝑖)*PSlash(idx(1,2),idx(1,1))' ||
      propagator.denominator !== '(P(idx(1,1)))^2' ||
      !gluonPropagator ||
      gluonPropagator.numerator !== '(-𝑖)*Metric(idx(1,1),idx(1,2))' ||
      gluonPropagator.denominator !== '(P(idx(1,1)))^2') {
    throw new Error(`${record.name}: effective local model rule differs`);
  }
  const hedge = (index) => `hedge(${index},1)`;
  const bis = (index) => `spenso::bis(4,${hedge(index)})`;
  const cof = (index) => `spenso::cof(3,${hedge(index)})`;
  const coad = (index) => `spenso::coad(8,${hedge(index)})`;
  const mink = (index) => `spenso::mink(4,${hedge(index)})`;
  for (const vertex of record.vertices) {
    const incoming = quarks.filter((edge) => edge.to === vertex.id);
    const outgoing = quarks.filter((edge) => edge.from === vertex.id);
    if (incoming.length !== 1 || outgoing.length !== 1) {
      throw new Error(`${record.name}: directed quark incidence differs`);
    }
    const input = incoming[0].toPort;
    const output = outgoing[0].fromPort;
    const g = gluons[0].from === vertex.id ? gluons[0].fromPort : gluons[0].toPort;
    const expected = `UFO::GC_11*spenso::gamma(${bis(output)},${bis(input)},${mink(g)})` +
      `*spenso::t(${coad(g)},${cof(input)},spenso::dind(${cof(output)}))`;
    if (vertex.int_id !== ruleName || vertex.num !== expected) {
      throw new Error(`${record.name}: complete model-derived vertex coefficient differs`);
    }
  }
  for (const edge of quarks) {
    const source = edge.fromPort;
    const target = edge.toPort;
    const momentum = `spenso::mink(4,edge(${edge.id},1))`;
    const expected = `(𝑖)*Q(${edge.id},${momentum})*spenso::g(${cof(source)},` +
      `spenso::dind(${cof(target)}))*spenso::gamma(${bis(target)},${bis(source)},${momentum})`;
    if (edge.num !== expected) {
      throw new Error(`${record.name}: complete model-derived quark propagator differs`);
    }
  }
  const gluon = gluons[0];
  const expectedGluon = `(-𝑖)*spenso::g(${coad(gluon.fromPort)},${coad(gluon.toPort)})` +
    `*spenso::g(${mink(gluon.fromPort)},${mink(gluon.toPort)})`;
  if (gluon.num !== expectedGluon) {
    throw new Error(`${record.name}: complete model-derived gluon propagator differs`);
  }
}

function dotRouting(expression, loops) {
  const vector = Array(loops).fill(0);
  for (const term of expression.split('+')) {
    const match = term.match(/^(?:(-?\d+)\*)?K\((\d+),a___\)$/);
    if (!match || Number(match[2]) >= loops) throw new Error(`invalid exact DOT routing ${expression}`);
    vector[Number(match[2])] += match[1] === undefined ? 1 : Number(match[1]);
  }
  return vector;
}

function conserveVacuumMomentum(record) {
  const net = new Map(record.vertices.map((vertex) => [vertex.id, [0, 0]]));
  for (const edge of record.edges) {
    const vector = dotRouting(edge.lmb_rep, 2);
    for (let index = 0; index < 2; index++) {
      net.get(edge.from)[index] -= vector[index];
      net.get(edge.to)[index] += vector[index];
    }
  }
  if ([...net.values()].some((vector) => vector[0] !== 0 || vector[1] !== 0)) {
    throw new Error(`${record.name}: complete signed vacuum momentum conservation differs`);
  }
}

function bindDotRouting(diagram, record, context) {
  if (diagram.edges.length !== record.edges.length) throw new Error(`${context}: routing edge count`);
  const basis = diagram.loop_momentum_basis;
  if (basis.external_edges.length || basis.dependent_externals.length) {
    throw new Error(`${context}: unexpected external shift`);
  }
  for (const [index, edge] of record.edges.entries()) {
    const signature = basis.edge_signatures[String(index)];
    if (edge.edgeId !== index || !signature || signature.external.length ||
        JSON.stringify(dotRouting(edge.lmb_rep, basis.loop_edges.length)) !==
          JSON.stringify(signature.loops)) {
      throw new Error(`${context}: DOT and native signed edge frame differ at ${index}`);
    }
  }
}

function automorphisms(record) {
  const vertices = record.vertices;
  const edges = record.edges;
  const permutations = (items) => items.length < 2 ? [items] : items.flatMap((item, index) =>
    permutations(items.filter((_, offset) => offset !== index)).map((tail) => [item, ...tail]));
  const edgeKey = (edge, index) => {
    let from = index[edge.from];
    let to = index[edge.to];
    if (edge.dir === 'back') [from, to] = [to, from];
    if (edge.dir === 'none' && from > to) [from, to] = [to, from];
    return JSON.stringify([from, to, edge.particle, edge.dir === 'none']);
  };
  const factorial = (n) => n < 2 ? 1 : n * factorial(n - 1);
  const original = Object.fromEntries(vertices.map((vertex, index) => [vertex.id, index]));
  const expected = edges.map((edge) => edgeKey(edge, original)).sort().join('|');
  let count = 0;
  for (const order of permutations(vertices)) {
    if (order.some((vertex, index) => vertex.int_id !== vertices[index].int_id)) continue;
    const index = Object.fromEntries(order.map((vertex, position) => [vertex.id, position]));
    const keys = edges.map((edge) => edgeKey(edge, index)).sort();
    if (keys.join('|') !== expected) continue;
    const multiplicities = new Map();
    for (const key of keys) multiplicities.set(key, (multiplicities.get(key) ?? 0) + 1);
    const selfLoopFlips = edges.filter((edge) => edge.dir === 'none' &&
      edge.from === edge.to).length;
    count += 2 ** selfLoopFlips * [...multiplicities.values()].reduce((product, multiplicity) =>
      product * factorial(multiplicity), 1);
  }
  return count;
}

function grassmannLoops(record, model) {
  const particles = new Map(model.particles.map((particle) => [particle.name, particle]));
  const lines = record.edges.filter((edge) => {
    const particle = particles.get(edge.particle);
    if (!particle) throw new Error(`unresolved particle ${edge.particle}`);
    return particle.spin === 2 || particle.ghost_number !== 0;
  });
  const parent = new Map(record.vertices.map((vertex) => [vertex.id, vertex.id]));
  const find = (id) => parent.get(id) === id ? id : find(parent.get(id));
  const degree = new Map();
  const incoming = new Map();
  const outgoing = new Map();
  for (const line of lines) {
    if (line.dir === 'none') throw new Error(`${record.name}: undirected Grassmann carrier`);
    const source = line.dir === 'back' ? line.to : line.from;
    const target = line.dir === 'back' ? line.from : line.to;
    outgoing.set(source, (outgoing.get(source) ?? 0) + 1);
    incoming.set(target, (incoming.get(target) ?? 0) + 1);
    parent.set(find(line.from), find(line.to));
    degree.set(line.from, (degree.get(line.from) ?? 0) + 1);
    degree.set(line.to, (degree.get(line.to) ?? 0) + 1);
  }
  const components = new Map();
  for (const line of lines) {
    const root = find(line.from);
    if (!components.has(root)) components.set(root, new Set());
    components.get(root).add(line.from).add(line.to);
  }
  for (const vertices of components.values()) {
    if ([...vertices].some((vertex) => degree.get(vertex) !== 2)) {
      throw new Error(`branching anticommuting lines in ${record.name}`);
    }
    if ([...vertices].some((vertex) => incoming.get(vertex) !== 1 ||
        outgoing.get(vertex) !== 1)) {
      throw new Error(`${record.name}: directed Grassmann circuit is not closed`);
    }
  }
  return components.size;
}

function dotPath(root, label, name) {
  const directory = path.join(root, `${label}-dot/processes/amplitudes`);
  const matches = fs.readdirSync(directory, {recursive: true})
    .filter((entry) => entry === `${name}.dot` || entry.endsWith(`/${name}.dot`))
    .map((entry) => path.join(directory, entry));
  if (matches.length !== 1) throw new Error(`expected one ${name}.dot in ${directory}, found ${matches.length}`);
  return matches[0];
}

function read(root, label) {
  const lines = fs.readFileSync(path.join(root, `${label}-graphs.tsv`), 'utf8').trimEnd().split('\n');
  const summary = lines.shift();
  const files = lines.map((line) => dotPath(root, label, line.split('\t')[1]));
  const parsed = liveDots(files);
  const records = lines.map((line, index) => {
    const [sourceIndex, name, factor, weight, groupId, master, globalPrefactor] = line.split('\t');
    if (!/^(0|[1-9]\d*)$/.test(sourceIndex) ||
        !Number.isSafeInteger(Number(sourceIndex))) {
      throw new Error(`${label}/${name}: invalid source occurrence ID`);
    }
    const dot = files[index];
    if (!fs.existsSync(dot)) throw new Error(`missing graph DOT: ${dot}`);
    const captured = graph(dot, parsed[index]);
    if (parsed[index].name !== name) throw new Error(`${label}/${name}: parsed DOT name differs`);
    if (captured.edges.some((edge) => edge.is_cut !== undefined)) {
      throw new Error(`${label}/${name}: vacuum has unexpected cut carrier`);
    }
    return {sourceIndex: Number(sourceIndex), name, factor, weight, groupId, master,
      globalPrefactor, ...captured};
  });
  return {summary, records};
}

const routes = Object.fromEntries([['source43', historical], ['current', current]].map(([route, root]) =>
  [route, Object.fromEntries(labels.map((label) => [label, read(root, label)]))]));
const captureRoots = {source43: historical, current};
const native = {};
const globalOwners = [];
const owner = (expression, ports) => ({expression, ...(ports === undefined ? {} : {interface: ports})});
for (const route of Object.values(routes)) for (const row of Object.values(route)) {
  for (const record of row.records) globalOwners.push(
    owner(record.globalPrefactor),
    owner(record.globalNumerator, record.numeratorPrefactorInterface),
    owner(record.projector, record.projectorInterface));
}
const equal = (left, right, message) => {
  if (JSON.stringify(left) !== JSON.stringify(right)) throw new Error(message);
};
const capturedFingerprint = (root) => {
  const value = liveModel(root);
  if (!/^[0-9a-f]{64}$/.test(value)) throw new Error('invalid captured live model fingerprint');
  return [...Buffer.from(value, 'hex')];
};
for (const [route, root] of [['source43', historical], ['current', current]]) {
  const model = JSON.parse(fs.readFileSync(path.join(root, 'post-import-model.json')));
  for (const filter of ['all', 'nonfactorizable']) {
    for (const record of routes[route][`${filter}-raw`].records) {
      conserveVacuumMomentum(record);
      record.directedAutomorphisms = automorphisms(record);
      record.grassmannLoops = grassmannLoops(record, model);
      const aut = record.directedAutomorphisms;
      const independent = rational(`${record.grassmannLoops % 2 ? -1 : 1}/${aut}`);
      const expected = `${record.grassmannLoops % 2 ? 'InternalFermionLoopSign(-1)*' : ''}` +
        `ExternalFermionOrderingSign(1)/AutG(${aut})`;
      scalarContract(record.overallFactor, expected, `${route}/${filter}/${record.name}/DOT`);
      const parsed = complete(record.factor);
      const dot = complete(record.overallFactor);
      if (JSON.stringify(shape(parse(record.factor))) !== JSON.stringify(shape(parse(expected))) ||
          !sameFactor(parsed.value, independent) || !sameFactor(dot.value, independent)) {
        throw new Error(`${route}/${filter}/${record.name}: complete independent weight factor differs`);
      }
      if (record.weight !== fractionString(independent)) {
        throw new Error(`${route}/${filter}/${record.name}: independent rational weight differs`);
      }
    }
    for (const record of routes[route][`${filter}-grouped`].records) {
      conserveVacuumMomentum(record);
    }
  }
}
for (const [route, root] of Object.entries(captureRoots)) {
  native[route] = {};
  const model = JSON.parse(fs.readFileSync(path.join(root, 'post-import-model.json')));
  const liveModelFingerprint = capturedFingerprint(root);
  for (const filter of ['all', 'nonfactorizable']) {
    const label = `${filter}-grouped`;
    const source = JSON.parse(fs.readFileSync(path.join(root, `${label}-native-source.json`)));
    const final = JSON.parse(fs.readFileSync(path.join(root, `${label}-native-final.json`)));
    const grouping = JSON.parse(fs.readFileSync(path.join(root, `${label}-native-groups.json`)));
    const raw = routes[route][`${filter}-raw`].records;
    const grouped = routes[route][label].records;
    if (source.length !== raw.length || final.length !== grouped.length ||
        grouping.groups.length !== final.length) throw new Error(`${route}/${filter}: native inventory size`);
    for (const [stage, diagrams] of [['source', source], ['final', final]]) {
      for (const [index, diagram] of diagrams.entries()) {
        globalOwners.push(owner(diagram.numerator_prefactor, diagram.numerator_prefactor_interface),
          owner(diagram.projector, diagram.projector_interface));
        complete(diagram.overall_factor);
        equal(diagram.model, liveModelFingerprint,
          `${route}/${filter}/${stage}/${index}: diagram model is not the captured post-import model`);
      }
    }
    for (const [index, diagram] of source.entries()) {
      checkAggregate(diagram, nativeCliExpression, `${route}/${filter}/source/${index}`);
    }
    const sourceText = fs.readFileSync(path.join(root, `${label}-native-source.json`), 'utf8');
    const sourceIds = [...sourceText.matchAll(/"id":\s*(\d+)/g)].map((match) => match[1]);
    if (grouping.source.length !== source.length || sourceIds.length !== source.length ||
        new Set(sourceIds).size !== source.length ||
        new Set(source.map((diagram) => diagram.name)).size !== source.length ||
        new Set(raw.map((record) => record.key)).size !== raw.length) {
      throw new Error(`${route}/${filter}: raw occurrence IDs or physical provenance are not unique`);
    }
    const sourceOccurrences = grouping.groups.flatMap((group) =>
      group.members.map((member) => member.source_diagram));
    equal([...sourceOccurrences].sort((a, b) => a - b), source.map((_, index) => index),
      `${route}/${filter}: grouped member occurrences are not a bijection of raw sources`);
    equal(grouping.groups.map((group) => group.master).sort((a, b) => a - b),
      final.map((_, index) => index),
      `${route}/${filter}: final master indices are not a bijection`);
    for (const [index, diagram] of source.entries()) {
      if (raw[index].sourceIndex !== index ||
          raw[index].groupId !== `Some(GroupId(${index}))` ||
          raw[index].master !== 'true' || diagram.name !== raw[index].name ||
          nativePhysical(diagram, model).key !== raw[index].key ||
          diagram.cuts.length !== 0 || diagram.topology_threshold_candidates.length !== 0 ||
          grouping.source[index].index !== index || grouping.source[index].name !== diagram.name ||
          JSON.stringify(shape(parse(diagram.overall_factor))) !==
            JSON.stringify(shape(parse(raw[index].overallFactor))) ||
          JSON.stringify(shape(parse(grouping.source[index].factor))) !==
            JSON.stringify(shape(parse(raw[index].factor)))) {
        throw new Error(`${route}/${filter}: source ${index} does not bind to the raw physical graph`);
      }
      scalarContract(diagram.overall_factor, raw[index].factor, `${route}/${filter}/source/${index}`);
      complete(grouping.source[index].factor);
      bindDotRouting(diagram, raw[index], `${route}/${filter}/source/${index}`);
      bindCliPayload(diagram, raw[index], model, `${route}/${filter}/source/${index}`);
    }
    for (const [index, group] of grouping.groups.entries()) {
      const master = final[group.master];
      const output = grouped[group.master];
      const selected = group.members.filter((member) => member.source_name === master.name);
      if (selected.length !== 1 || output.sourceIndex !== group.master ||
          output.groupId !== `Some(GroupId(${index}))` || output.master !== 'true' ||
          master.name !== output.name ||
          master.cuts.length !== 0 || master.topology_threshold_candidates.length !== 0 ||
          nativePhysical(master, model).key !== output.key ||
          master.name !== source[selected[0].source_diagram].name ||
          JSON.stringify(shape(parse(master.overall_factor))) !==
            JSON.stringify(shape(parse(output.overallFactor)))) {
        throw new Error(`${route}/${filter}: final master ${index} is not bound`);
      }
      bindDotRouting(master, output, `${route}/${filter}/final/${index}`);
      bindCliPayload(master, output, model, `${route}/${filter}/final/${index}`);
      const {overall_factor: masterFactor, ...masterPayload} = master;
      const {overall_factor: sourceFactor, ...sourcePayload} = source[selected[0].source_diagram];
      equal(masterPayload, sourcePayload, `${route}/${filter}: final tensor, routing or frame differs from master`);
      const expectedFactor = group.members.length === 1 ?
        raw[group.members[0].source_diagram].factor :
        [...group.members].sort((a, b) => a.source_diagram - b.source_diagram)
          .map((member) => `NumeratorDependentGrouping(${member.source_diagram},${member.ratio},` +
            `${raw[member.source_diagram].factor})`).join('+');
      const groupedFactor = complete(output.factor);
      const dotFactor = complete(output.overallFactor);
      scalarContract(output.overallFactor, expectedFactor, `${route}/${filter}/final/${index}/DOT`);
      scalarContract(master.overall_factor, expectedFactor, `${route}/${filter}/final/${index}/native`);
      if (JSON.stringify(shape(parse(output.factor))) !== JSON.stringify(shape(parse(expectedFactor))) ||
          !sameFactor(groupedFactor.value, rational(output.weight)) ||
          !sameFactor(dotFactor.value, groupedFactor.value)) {
        throw new Error(`${route}/${filter}: complete final scalar factor differs at group ${index}`);
      }
      const groupIds = groupedFactor.groups.map((term) =>
        [term.sourceId, term.ratio, term.memberWeight])
        .sort((a, b) => a[0] - b[0]);
      if (group.members.length === 1 && groupIds.length ||
          group.members.length > 1 && groupIds.length !== group.members.length) {
        throw new Error(`${route}/${filter}: output group terms missing`);
      }
      for (const member of group.members) {
        complete(member.overall_factor);
        const input = source[member.source_diagram];
        if (!input || member.source_name !== input.name || member.ratio !== '1' ||
            member.overall_factor !== grouping.source[member.source_diagram].factor ||
            JSON.stringify(shape(parse(member.overall_factor))) !==
              JSON.stringify(shape(parse(raw[member.source_diagram].factor))) ||
            !raw.some((record) => record.name === input.name && record.key === nativePhysical(input, model).key)) {
          throw new Error(`${route}/${filter}: member identity, ratio or factor differs`);
        }
        if (group.members.length > 1 && !groupIds.some(([id, ratio, weight]) =>
            id === member.source_diagram && ratio === member.ratio &&
            sameFactor(rational(weight), rational(raw[id].weight)))) {
          throw new Error(`${route}/${filter}: source occurrence, ratio or weight absent from final factor`);
        }
      }
      const expectedWeight = sum(group.members.map((member) =>
        multiply(rational(raw[member.source_diagram].weight), rational(member.ratio))));
      if (output.weight !== fractionString(expectedWeight)) {
        throw new Error(`${route}/${filter}: final group ${index} rational weight differs`);
      }
      if (group.members.length > 1) {
        const first = source[group.members[0].source_diagram];
        for (const member of group.members.slice(1)) {
          const other = source[member.source_diagram];
          equal(other.vertices.map((vertex) => [vertex.numerator, vertex.num_interface]),
            first.vertices.map((vertex) => [vertex.numerator, vertex.num_interface]),
            `${route}/${filter}: member local vertex tensors differ`);
          equal(other.edges.map(([, edge]) => [edge.numerator, edge.num_interface]),
            first.edges.map(([, edge]) => [edge.numerator, edge.num_interface]),
            `${route}/${filter}: member local edge tensors differ`);
          if (other.numerator !== first.numerator ||
              JSON.stringify(other.num_interface) !== JSON.stringify(first.num_interface)) {
            throw new Error(`${route}/${filter}: complete factorized member tensor differs`);
          }
        }
      }
    }
    const light = raw.filter((record) => record.name === 'GL2' || record.name === 'GL4');
    if (light.length !== 2) throw new Error(`${route}/${filter}: missing d/u light pair`);
    for (const record of light) bindLightLocalFactors(record, model);
    const incidence = light.map((record) => {
      const expression = [...record.vertices.map((vertex) => vertex.num),
        ...record.edges.map((edge) => edge.num)].join('*');
      const axes = new Map();
      for (const match of expression.matchAll(
        /spenso::(bis|cof|coad|mink)\((\d+),(hedge|edge)\((\d+),1\)\)/g)) {
        const [, representation, dimension, owner, index] = match;
        if (!({bis: '4', cof: '3', coad: '8', mink: '4'})[representation] ||
            ({bis: '4', cof: '3', coad: '8', mink: '4'})[representation] !== dimension) {
          throw new Error(`${route}/${filter}/${record.name}: contracted dimension changed`);
        }
        const key = `${representation}:${dimension}:${owner}:${index}`;
        axes.set(key, (axes.get(key) ?? 0) + 1);
      }
      if (axes.size !== 14 || [...axes.values()].some((count) => count !== 2)) {
        throw new Error(`${route}/${filter}/${record.name}: complete 4D/SU(3) contraction differs`);
      }
      return {name: record.name, axes: Object.fromEntries([...axes].sort())};
    });
    equal(incidence[0].axes, incidence[1].axes,
      `${route}/${filter}: d/u contracted index domains or owners differ`);
    native[route][filter] = {sourceCount: source.length, finalCount: final.length,
      rawTotal: fractionString(sum(raw.map((record) => rational(record.weight)))),
      groupedTotal: fractionString(sum(grouped.map((record) => rational(record.weight)))),
      singletonConstant: fractionString(sum(grouping.groups.flatMap((group, index) =>
        group.members.length === 1 ? [rational(grouped[group.master].weight)] : []))),
      lightPairAxes: incidence,
      groups: grouping.groups.map((group) => ({master: final[group.master].name,
        members: group.members.map(({source_diagram, source_name, ratio}) =>
          ({source_diagram, source_name, ratio}))}))};
    for (const [kind, records, expected] of [
      ['raw', raw, native[route][filter].rawTotal],
      ['grouped', grouped, native[route][filter].groupedTotal],
    ]) {
      const match = routes[route][`${filter}-${kind}`].summary.match(/^(\d+) \| (.*) = (\S+)$/);
      const expectedExpression = kind === 'raw' ? expected :
        native[route][filter].singletonConstant + grouping.groups
          .filter((group) => group.members.length > 1)
          .flatMap((group) => [...group.members].sort((a, b) =>
            a.source_diagram - b.source_diagram)
            .map((member) => `Group(${member.source_diagram},${member.ratio},` +
              `${raw[member.source_diagram].weight})`)).map((term) => `+${term}`).join('');
      if (!match || Number(match[1]) !== records.length || match[3] !== expected ||
          JSON.stringify(shape(parse(match[2]))) !== JSON.stringify(shape(parse(expectedExpression))) ||
          !sameFactor(complete(match[2]).value, rational(expected))) {
        throw new Error(`${route}/${filter}/${kind}: independently summed total differs from CLI`);
      }
    }
  }
}
const graphKeys = (records) => records.map((record) => record.key).sort();
for (const filter of ['all', 'nonfactorizable']) {
  const a = graphKeys(routes.source43[`${filter}-raw`].records);
  const b = graphKeys(routes.current[`${filter}-raw`].records);
  if (JSON.stringify(a) !== JSON.stringify(b)) throw new Error(`${filter}: physical raw graph inventory changed`);
  for (const label of [`${filter}-raw`, `${filter}-grouped`]) {
    const old = routes.source43[label].records;
    const fresh = routes.current[label].records;
    const byKey = new Map(fresh.map((record) => [record.key, record]));
    for (const record of old) {
      const other = byKey.get(record.key);
      if (!other || record.weight !== other.weight || record.name !== other.name ||
          record.factor !== other.factor || record.globalPrefactor !== other.globalPrefactor ||
          record.globalNumerator !== other.globalNumerator || record.projector !== other.projector ||
          record.vertices.length !== other.vertices.length ||
          record.edges.length !== other.edges.length) {
        throw new Error(`${label}: producer output differs for ${record.name}`);
      }
      for (const [index, vertex] of record.vertices.entries()) {
        for (const field of ['int_id', 'num']) {
          if (vertex[field] !== other.vertices[index][field]) {
            throw new Error(`${label}/${record.name}: vertex ${index} ${field} differs`);
          }
        }
      }
      for (const [index, edge] of record.edges.entries()) {
        for (const field of ['from', 'to', 'fromPort', 'toPort', 'particle', 'dir', 'lmb_rep', 'num']) {
          if (edge[field] !== other.edges[index][field]) {
            throw new Error(`${label}/${record.name}: edge ${index} ${field} differs`);
          }
        }
      }
    }
  }
}
const globalBinding = liveGlobals(globalOwners);
const liveRules = Object.fromEntries(Object.entries(captureRoots).map(([route, root]) =>
  [route, liveRule(root, ['all', 'nonfactorizable'].flatMap((filter) =>
    ['source', 'final'].map((stage) => path.join(root, `${filter}-grouped-native-${stage}.json`))))]));
fs.writeFileSync(output, JSON.stringify({routes, native, liveRules, globalBinding}, null, 2) + '\n');
console.log('physical source IDs, local tensors, ratios, final masters and raw payloads verified');
