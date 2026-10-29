import assert from 'node:assert/strict';
import crypto from 'node:crypto';
import fs from 'node:fs';
import path from 'node:path';

const [currentRoot, approvedRoot, outputPath] = process.argv.slice(2);
assert(currentRoot && approvedRoot && outputPath,
  'usage: compare-current-sunset.mjs <current-output-root> <approved-output-root> <output.json>');
const comparisons = [];
for (const name of ['inputs', 'model', 'join']) {
  const currentPath = path.join(currentRoot, `${name}.json`);
  const approvedPath = path.join(approvedRoot, `sunset-${name}.json`);
  const currentBytes = fs.readFileSync(currentPath);
  const approvedBytes = fs.readFileSync(approvedPath);
  assert.deepEqual(JSON.parse(currentBytes), JSON.parse(approvedBytes), `${name} semantic mismatch`);
  assert.equal(Buffer.compare(currentBytes, approvedBytes), 0, `${name} byte mismatch`);
  comparisons.push({
    name,
    currentPath,
    approvedPath,
    bytes: currentBytes.length,
    sha256: crypto.createHash('sha256').update(currentBytes).digest('hex'),
  });
}

const approvedJoin = JSON.parse(fs.readFileSync(path.join(approvedRoot, 'sunset-join.json')));
const mutation = structuredClone(JSON.parse(fs.readFileSync(path.join(currentRoot, 'join.json'))));
mutation.routes.current['all-grouped'].summary = '__DISPOSABLE_CURRENT_MUTATION__';
let mutationRejected = false;
try {
  assert.deepEqual(mutation, approvedJoin);
} catch {
  mutationRejected = true;
}
assert(mutationRejected, 'disposable current sunset mutation was not detected');
const report = {
  schema: 1,
  exact: true,
  comparisons,
  mutation: {
    target: 'current all-grouped summary',
    frozenPreparationInputUnchanged: true,
    rejected: mutationRejected,
  },
};
fs.writeFileSync(outputPath, `${JSON.stringify(report, null, 2)}\n`);
console.log('PASS current inputs/model/join equal approved records; mutation rejected');
