import fs from 'node:fs';
import crypto from 'node:crypto';
import {spawnSync} from 'node:child_process';
import {complete, format} from './factor.mjs';

const [output] = process.argv.slice(2);
if (!output) throw new Error('usage: node scalar-differential.mjs <output.json>');
const oracleRoot = '/shared/lcnbr/dev/gammaloop/generation-pocket-weight-diagnosis/target/ghost-weight-certificate/round3';
const oracle = `${oracleRoot}/scalar-oracle`;
const oracleSource = `${oracleRoot}/scalar-oracle.rs`;
const hash = (file) => crypto.createHash('sha256').update(fs.readFileSync(file)).digest('hex');
if (hash(oracle) !== 'ebd36d06ccdd1f21877f715feed76bc9d329e3b4643a493eba77e1d687eff343' ||
    hash(oracleSource) !== 'fac6f5522314ce3057af3a3b87de6911283d46264ace241371a60cf2fb4b5ae4') {
  throw new Error('pinned Symbolica scalar oracle identity changed');
}
const accepted = [
  '0', '1', '-1', '--2', '2*3', '2*-3', '-2*3', '2/-3', '2/3*4', '8/4/2',
  '(8/4)/2', '8/(4/2)', '2*(3/4)', '0/1',
  '2^0', '2^1', '2^-1', '2^(+1)', '2^(-1)', '2^+1',
  '-2^0', '(-2)^0', '(-2)^1', '(-2)^-1', '-2^-1', '--2^0',
  '-(2^0)', '-(-2)^(-1)', '1/-2^0', '1/(-2)^0', '1/(-2^0)',
  '0^1', '1^0', '1^(-1)', '1+2', '1-2', '1+-2', '1+7',
  '-19/24+1/8', '9007199254740993/3',
];
const unsupported = [
  '2^1^0', '2^0^0', '2^(-1)^0', '2^(1^0)', '0^0', '0^-1', '1/0',
  '2^2', '2^-2', '+1', '2/+3', '2.0', 'x', 'Bogus(2)',
  'Bogus::AutG(2)', '1 trailing', '1(2)',
];
const all = [...accepted, ...unsupported];
const run = spawnSync(oracle, all, {encoding: 'utf8', env: process.env,
  maxBuffer: 1024 * 1024});
if (run.status !== 0) throw new Error(`pinned oracle failed: ${run.stderr.slice(0, 600)}`);
const lines = run.stdout.trimEnd().split('\n').filter((line) => line.includes('\t'));
if (lines.length !== all.length) throw new Error('pinned oracle row count differs');
const rows = all.map((expression, index) => {
  const [returned, authoritative] = lines[index].split('\t');
  if (returned !== expression || authoritative === undefined) {
    throw new Error(`pinned oracle row ${index} differs`);
  }
  if (index < accepted.length) {
    const actual = format(complete(expression).value);
    if (actual !== authoritative) {
      throw new Error(`scalar mismatch for ${expression}: ${actual} versus ${authoritative}`);
    }
    return {expression, admitted: true, authoritative, checker: actual};
  }
  let rejection;
  try {complete(expression);} catch (error) {rejection = error.message;}
  if (!rejection) throw new Error(`unsupported scalar syntax escaped: ${expression}`);
  return {expression, admitted: false, authoritative, rejection};
});
const result = {schema: 1, kernel: '1a1829b1f3d0a88b53b44085b3b9f34ef77945b3',
  oracle: {path: oracle, sha256: hash(oracle)},
  oracleSource: {path: oracleSource, sha256: hash(oracleSource)},
  accepted: accepted.length, explicitlyRejected: unsupported.length, rows};
fs.writeFileSync(output, JSON.stringify(result, null, 2) + '\n');
console.log(`pinned scalar differential: ${accepted.length} exact matches, ` +
  `${unsupported.length} explicit rejections`);
