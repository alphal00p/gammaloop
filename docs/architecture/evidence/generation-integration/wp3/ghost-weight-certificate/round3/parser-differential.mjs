import fs from 'node:fs';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {spawnSync} from 'node:child_process';
import {parseFactorExpression} from '../analyze-graphs.mjs';

const oracle = 'target/ghost-weight-certificate/round3/scalar-oracle';
const accepted = [
  '0', '1', '-1', '2', '-2', '--2',
  '2*3', '2*-3', '-2*3', '2/-3', '2/3*4', '2/(3*4)',
  '8/4/2', '(8/4)/2', '8/(4/2)', '2*(3/4)', '0/1',
  '2^0', '2^1', '2^-1', '2^(+1)', '2^(-1)', '2^+1',
  '-2^0', '(-2)^0', '(-2)^1', '(-2)^-1',
  '-2^-1', '-2^+1', '--2^0', '-(-2)^0',
  '-(2^0)', '-(-2)^(-1)', '(-2)^+1',
  '1/-2^0', '1/(-2)^0', '1/(-2^0)', '1/((-2)^0)',
  '0^1', '0^(+1)', '1^0', '1^(-1)',
  '9007199254740993/3', '-9007199254740993/3',
];
const unsupported = [
  '2^1^0', '2^(-1)^0', '2^(1^0)',
  '0^0', '0^-1', '1/0',
  '2^2', '2^-2',
  '1+2', '1-2', '2--3', '2+-3',
  '+1', '+2', '-+2', '+-2', '2/+3', '+2^0', '-+2^0', '+(-2)^0',
  '2.0', 'x', 'Bogus(2)', '2(3)',
];
const all = [...accepted, ...unsupported];
const result = spawnSync(oracle, all, {encoding: 'utf8', maxBuffer: 1024 * 1024,
  env: process.env});
assert.equal(result.status, 0, result.stderr);
const lines = result.stdout.trimEnd().split('\n').filter((line) => line.includes('\t'));
assert.equal(lines.length, all.length, 'pinned Symbolica scalar result count');
const ratio = (value) => `${value.n}/${value.d}`;
const rows = [];
for (const [index, expression] of all.entries()) {
  const [returned, authoritative] = lines[index].split('\t');
  assert.equal(returned, expression, 'oracle output order');
  if (index < accepted.length) {
    assert(/^-?\d+(\/\d+)?$/.test(authoritative), `oracle did not yield a rational: ${expression}`);
    const parsed = parseFactorExpression(expression).value;
    const expected = authoritative.includes('/') ? authoritative : `${authoritative}/1`;
    assert.equal(ratio(parsed), expected, `scalar precedence mismatch: ${expression}`);
    rows.push({expression, admitted: true, authoritative, checker: ratio(parsed)});
  } else {
    assert.throws(() => parseFactorExpression(expression), undefined,
      `unsupported scalar syntax escaped: ${expression}`);
    rows.push({expression, admitted: false, authoritative, checker: 'rejected explicitly'});
  }
}
const sha = (location) => crypto.createHash('sha256').update(fs.readFileSync(location)).digest('hex');
const report = {schema: 1, kernel: '1a1829b1f3d0a88b53b44085b3b9f34ef77945b3',
  oracleSource: {path: 'target/ghost-weight-certificate/round3/scalar-oracle.rs',
    sha256: sha('target/ghost-weight-certificate/round3/scalar-oracle.rs')},
  oracleBinary: {path: oracle, sha256: sha(oracle)},
  pinnedSymbolicaRlib: {path: '/common/build-cache/lcnbr/9c/5e5e2f7a370b0f/ci-optim/deps/libsymbolica-219a1a53d9363ef7.rlib',
    sha256: sha('/common/build-cache/lcnbr/9c/5e5e2f7a370b0f/ci-optim/deps/libsymbolica-219a1a53d9363ef7.rlib')},
  admittedCount: accepted.length, explicitlyRejectedCount: unsupported.length,
  parserScope: 'exact rational literals, unary signs, left-associated multiplication/division, parentheses, and one exponent -1/0/1; additive forms, chained powers, irrational/undefined powers, arbitrary symbols and implicit multiplication rejected',
  rows};
fs.writeFileSync('target/ghost-weight-certificate/round3/parser-differential.json',
  `${JSON.stringify(report, null, 2)}\n`);
console.log(`pinned scalar differential: ${accepted.length} admitted values equal Symbolica; ${unsupported.length} explicitly rejected`);
