import {rational, multiply, format} from './factor.mjs';

// Concrete indexed component factors commute. Keep every local sum intact: this
// normal form never distributes a product over a sum or performs a contraction.
function split(value, delimiter) {
  const parts = [];
  let start = 0;
  let depth = 0;
  for (let index = 0; index < value.length; index++) {
    if (value[index] === '(') depth++;
    else if (value[index] === ')') depth--;
    else if (value[index] === delimiter && depth === 0) {
      parts.push(value.slice(start, index));
      start = index + 1;
    }
    if (depth < 0) throw new Error('unbalanced tensor product');
  }
  if (depth) throw new Error('unbalanced tensor product');
  parts.push(value.slice(start));
  return parts;
}
const key = (value) => JSON.stringify(value);
const sorted = (values) => values.sort((a, b) => key(a).localeCompare(key(b)));
const scalar = (n, phase = 0) => ({coefficient: rational(n), phase, factors: []});
function product(values) {
  let coefficient = rational(1n);
  let phase = 0;
  const factors = new Map();
  for (const value of values) {
    coefficient = multiply(coefficient, value.coefficient);
    phase += value.phase;
    for (const [factor, power] of value.factors) {
      const id = key(factor);
      const previous = factors.get(id);
      factors.set(id, [factor, power + (previous?.[1] ?? 0)]);
    }
  }
  if (phase % 4 >= 2) coefficient = multiply(coefficient, rational(-1n));
  phase %= 2;
  return {coefficient, phase, factors: sorted([...factors.values()])};
}
const atomic = (value) => ({coefficient: rational(1n), phase: 0, factors: [[value, 1]]});
const packed = (value) => [format(value.coefficient), value.phase, value.factors];
function parse(source) {
  let value = source;
  while (value.startsWith('(') && value.endsWith(')')) {
    let depth = 0;
    let encloses = true;
    for (let index = 0; index < value.length - 1; index++) {
      if (value[index] === '(') depth++;
      if (value[index] === ')') depth--;
      if (depth === 0) {encloses = false; break;}
    }
    if (!encloses) break;
    value = value.slice(1, -1);
  }
  const sum = split(value, '+');
  if (sum.length > 1) return atomic(['sum', ...sorted(sum.map((term) => packed(parse(term))))]);
  const factors = split(value, '*');
  if (factors.length > 1) return product(factors.map(parse));
  if (value.startsWith('-')) return product([scalar(-1n), parse(value.slice(1))]);
  const power = split(value, '^');
  if (power.length === 2 && /^[1-4]$/.test(power[1])) {
    return product(Array.from({length: Number(power[1])}, () => parse(power[0])));
  }
  if (power.length > 1) throw new Error('unsupported tensor scalar power');
  if (/^\d+$/.test(value)) return scalar(BigInt(value));
  if (value === '𝑖') return scalar(1n, 1);
  const open = value.indexOf('(');
  if (open > 0 && value.endsWith(')')) {
    const head = value.slice(0, open);
    let args = split(value.slice(open + 1, -1), ',').map((arg) => packed(parse(arg)));
    if (head === 'spenso::g') args = sorted(args);
    if (head === 'spenso::f' && args.length === 3) {
      args = sorted([args, [args[1], args[2], args[0]], [args[2], args[0], args[1]]])[0];
    }
    return atomic(['call', head, args]);
  }
  if (!/^[A-Za-z_][A-Za-z0-9_:]*$/.test(value)) {
    throw new Error(`unsupported complete tensor atom: ${value}`);
  }
  return atomic(['atom', value]);
}
export const tensorKey = (source) => key(packed(parse(source)));

export function checkAggregate(diagram, normalize, where) {
  if (diagram.numerator_prefactor !== '1' || diagram.projector !== '1' ||
      (diagram.num_interface?.length ?? 0) !== 0) {
    throw new Error(`${where}: vacuum prefactor, projector or closed aggregate interface differs`);
  }
  const fragments = [...diagram.edges.map(([, edge]) => edge.numerator),
    ...diagram.vertices.map((vertex) => vertex.numerator)];
  const lhs = tensorKey(normalize(diagram.numerator));
  const rhs = tensorKey(fragments.map((source) => `(${normalize(source)})`).join('*'));
  if (lhs !== rhs) throw new Error(`${where}: aggregate differs from complete local product`);
}

const nativeCliHeads = [
  ['UFO::{}::', 'UFO::'],
  ['gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q', 'Q'],
  ['gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::Q', 'Q'],
  ['gammalooprs::{spenso::index,spenso::index-label:e}::', ''],
  ['gammalooprs::{spenso::index,spenso::index-label:v}::', ''],
  ['gammalooprs::{spenso::index}::', ''],
  ['spenso::{spenso::representation,spenso::self_dual}::', 'spenso::'],
  ['spenso::{spenso::representation,spenso::dualizable}::', 'spenso::'],
  ['spenso::{spenso::lower}::', 'spenso::'],
  ['spenso::{antisymmetric,real,spenso::tensor}::', 'spenso::'],
  ['spenso::{linear,spenso::tensor}::', 'spenso::'],
  ['spenso::{real,spenso::tensor}::', 'spenso::'],
  ['spenso::{symmetric,linear,real,spenso::tensor}::', 'spenso::'],
];
export function nativeCliExpression(source) {
  const result = nativeCliHeads.reduce((value, [native, cli]) => value.replaceAll(native, cli), source);
  if (result.includes('::{}::') || result.includes('::{')) {
    throw new Error(`unsupported native-to-CLI tensor head in ${source}`);
  }
  return result;
}
