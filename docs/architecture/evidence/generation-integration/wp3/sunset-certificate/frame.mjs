import fs from 'node:fs';
import path from 'node:path';

const [historical, current, output] = process.argv.slice(2);
if (!historical || !current || !output) {
  throw new Error('usage: node frame.mjs <historical-capture> <current-capture> <output.json>');
}
const equal = (left, right, message) => {
  if (JSON.stringify(left) !== JSON.stringify(right)) throw new Error(message);
};
const normalizedParticle = (name) => ['d', 'u'].includes(name) ? 'massless-light-quark' : name;
const result = {};
for (const [route, root] of [['source43', historical], ['current', current]]) {
  result[route] = {};
  const model = JSON.parse(fs.readFileSync(path.join(root, 'post-import-model.json')));
  for (const filter of ['all', 'nonfactorizable']) {
    const label = `${filter}-grouped`;
    const source = JSON.parse(fs.readFileSync(path.join(root, `${label}-native-source.json`)));
    const final = JSON.parse(fs.readFileSync(path.join(root, `${label}-native-final.json`)));
    const groups = JSON.parse(fs.readFileSync(path.join(root, `${label}-native-groups.json`))).groups;
    const records = [];
    for (const group of groups) {
      const master = final[group.master];
      for (const member of group.members) {
        const diagram = source[member.source_diagram];
        const context = `${route}/${filter}/${diagram.name}->${master.name}`;
        const oldBasis = diagram.loop_momentum_basis;
        const newBasis = master.loop_momentum_basis;
        if (oldBasis.external_edges.length || newBasis.external_edges.length ||
            oldBasis.dependent_externals.length || newBasis.dependent_externals.length) {
          throw new Error(`${context}: unexpected external affine shift`);
        }
        if (diagram.edges.length !== master.edges.length ||
            oldBasis.loop_edges.length !== newBasis.loop_edges.length ||
            oldBasis.loop_edges.length !== 2) {
          throw new Error(`${context}: frame dimensions differ`);
        }
        const edgeMap = diagram.edges.map(([endpoints, edge], index) => {
          const [targetEndpoints, target] = master.edges[index];
          const a = model.particles[edge.particle]?.name;
          const b = model.particles[target.particle]?.name;
          if (normalizedParticle(a) !== normalizedParticle(b) || edge.directed !== target.directed ||
              edge.source_slot !== target.source_slot || edge.target_slot !== target.target_slot ||
              edge.orientation !== target.orientation ||
              edge.numerator !== target.numerator ||
              JSON.stringify(edge.num_interface) !== JSON.stringify(target.num_interface)) {
            throw new Error(`${context}: edge ${index} physical owner or local tensor differs`);
          }
          if (endpoints.source !== targetEndpoints.source || endpoints.target !== targetEndpoints.target) {
            throw new Error(`${context}: edge ${index} endpoint transport is not the selected map`);
          }
          return {source: index, target: index, sign: 1};
        });
        for (const [index, vertex] of diagram.vertices.entries()) {
          const target = master.vertices[index];
          if (vertex.numerator !== target.numerator ||
              JSON.stringify(vertex.num_interface) !== JSON.stringify(target.num_interface)) {
            throw new Error(`${context}: vertex ${index} tensor or axis differs`);
          }
        }
        const vector = (basis, edge) => {
          const signature = basis.edge_signatures[String(edge)];
          if (!signature || signature.external.length) throw new Error(`${context}: missing exact signature`);
          return signature.loops;
        };
        const matrix = oldBasis.loop_edges.map((edge) => vector(newBasis, edge));
        for (const {source: oldEdge, target: newEdge, sign} of edgeMap) {
          const old = vector(oldBasis, oldEdge);
          const transformed = [0, 1].map((column) => old.reduce((sum, coefficient, row) =>
            sum + coefficient * matrix[row][column], 0));
          const expected = vector(newBasis, newEdge).map((coefficient) => sign * coefficient);
          equal(transformed, expected, `${context}: signed affine routing differs on edge ${oldEdge}`);
        }
        records.push({source: diagram.name, master: master.name, edgeMap, matrix,
          sourceBasis: oldBasis, masterBasis: newBasis});
      }
    }
    result[route][filter] = records;
  }
}
fs.writeFileSync(output, JSON.stringify(result, null, 2) + '\n');
console.log('all within-producer signed affine maps and complete edge routing verified');
