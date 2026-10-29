import fs from 'node:fs';
import path from 'node:path';

const [capture, output] = process.argv.slice(2);
if (!capture || !output) throw new Error('usage: node model.mjs <capture> <output.json>');
const model = JSON.parse(fs.readFileSync(path.join(capture, 'post-import-model.json')));
const card = JSON.parse(fs.readFileSync(path.join(capture, 'post-import-parameters.json')));
const named = (items, name) => {
  const matches = items.filter((item) => item.name === name);
  if (matches.length !== 1) throw new Error(`expected one ${name}, found ${matches.length}`);
  return matches[0];
};
const particles = Object.fromEntries(['d', 'd~', 'u', 'u~'].map((name) =>
  [name, named(model.particles, name)]));
const zero = named(model.parameters, 'ZERO');
if (zero.nature !== 'internal' || JSON.stringify(zero.value) !== '[0,0]' ||
    Object.hasOwn(card, 'ZERO')) throw new Error('ZERO is not a fixed zero mass');
for (const [d, u] of [['d', 'u'], ['d~', 'u~']]) {
  const left = particles[d];
  const right = particles[u];
  for (const field of ['spin', 'color', 'mass', 'width', 'ghost_number']) {
    if (left[field] !== right[field]) throw new Error(`${d}/${u} ${field} differs`);
  }
  if (left.mass !== 'ZERO' || left.spin !== 2 || Math.abs(left.color) !== 3) {
    throw new Error(`${d}/${u} effective mass or spin/color domain differs`);
  }
  const prop = (name) => model.propagators.filter((candidate) => candidate.particle === name);
  const a = prop(d);
  const b = prop(u);
  if (a.length !== 1 || b.length !== 1 || a[0].numerator !== b[0].numerator ||
      a[0].denominator !== b[0].denominator) {
    throw new Error(`${d}/${u} propagator rules differ`);
  }
}
const dVertex = named(model.vertex_rules, 'V_74');
const uVertex = named(model.vertex_rules, 'V_135');
for (const field of ['color_structures', 'lorentz_structures', 'couplings']) {
  if (JSON.stringify(dVertex[field]) !== JSON.stringify(uVertex[field])) {
    throw new Error(`d/u gluon vertex ${field} differs`);
  }
}
if (JSON.stringify(dVertex.particles) !== '["d~","d","g"]' ||
    JSON.stringify(uVertex.particles) !== '["u~","u","g"]') {
  throw new Error('oriented quark and antiquark slots differ');
}
const gluon = named(model.particles, 'g');
const ghost = named(model.particles, 'ghG');
const top = named(model.particles, 't');
if (gluon.mass !== 'ZERO' || gluon.spin !== 3 || gluon.color !== 8 ||
    ghost.mass !== 'ZERO' || ghost.spin !== -1 || ghost.color !== 8 ||
    ghost.ghost_number !== 1 || top.mass !== 'MT' || top.spin !== 2 ||
    top.color !== 3 || JSON.stringify(card.MT) !== '[173,0]') {
  throw new Error('singleton gluon, ghost or top physical model input differs');
}
const ghostPropagator = named(model.propagators, 'ghG_propFeynman');
const ghostVertex = named(model.vertex_rules, 'V_35');
if (ghostPropagator.numerator !== '𝑖' ||
    ghostPropagator.denominator !== '(P(idx(1,1)))^2' ||
    JSON.stringify(ghostVertex.particles) !== '["ghG","ghG~","g"]' ||
    JSON.stringify(ghostVertex.couplings) !== '[["GC_10"]]' ||
    JSON.stringify(ghostVertex.lorentz_structures) !== '["UUV1"]' ||
    JSON.stringify(ghostVertex.color_structures) !== '["f(3,1,2)"]') {
  throw new Error('ghost propagator or vertex input differs');
}
const result = {restriction: model.restriction, zero, particles,
  propagators: ['d', 'd~', 'u', 'u~'].map((name) =>
    named(model.propagators.filter((propagator) => propagator.particle === name), `${name}_propFeynman`)),
  dVertex, uVertex,
  gluon, ghost, top, ghostPropagator, ghostVertex,
  gluonTripleVertex: named(model.vertex_rules, 'V_36'),
  gluonQuarticVertex: named(model.vertex_rules, 'V_37'),
  topGluonVertex: named(model.vertex_rules, 'V_137'),
  gluonCoupling: named(model.couplings, 'GC_11'),
  lorentz: named(model.lorentz_structures, 'FFV1')};
fs.writeFileSync(output, JSON.stringify(result, null, 2) + '\n');
console.log('live d/u masses, propagators, couplings and oriented vertex domains agree');
