import assert from 'node:assert/strict';
import crypto from 'node:crypto';
import fs from 'node:fs';

const [currentPath, approvedPath, outputPath] = process.argv.slice(2);
assert(currentPath && approvedPath && outputPath,
  'usage: compare-current-ghost.mjs <current-full.json> <approved-full.json> <output.json>');
const current = JSON.parse(fs.readFileSync(currentPath));
const approved = JSON.parse(fs.readFileSync(approvedPath));
const observedRoots = {current: current.currentRoot, historical: current.historicalRoot};
let normalizedPaths = 0;
for (const rows of Object.values(current.cases)) {
  for (const row of rows.current) {
    assert(row.dotPath.startsWith(observedRoots.current), 'current DOT path outside capture root');
    row.dotPath = approved.currentRoot + row.dotPath.slice(observedRoots.current.length);
    normalizedPaths += 1;
  }
}
current.currentRoot = approved.currentRoot;
current.historicalRoot = approved.historicalRoot;
assert.deepEqual(current, approved, 'current combined capture differs from approved physical join');

const mutation = structuredClone(current);
const firstCase = Object.values(mutation.cases)[0];
firstCase.current[0].evaluatedWeight = '__DISPOSABLE_CURRENT_MUTATION__';
let mutationRejected = false;
try {
  assert.deepEqual(mutation, approved);
} catch {
  mutationRejected = true;
}
assert(mutationRejected, 'disposable current-side mutation was not detected');
const canonical = `${JSON.stringify(current)}\n`;
const report = {
  schema: 1,
  exact: true,
  currentPath,
  approvedPath,
  observedRoots,
  normalizedPaths,
  normalizedSha256: crypto.createHash('sha256').update(canonical).digest('hex'),
  mutation: {
    target: 'first current evaluatedWeight',
    frozenPreparationInputUnchanged: true,
    rejected: mutationRejected,
  },
};
fs.writeFileSync(outputPath, `${JSON.stringify(report, null, 2)}\n`);
console.log(`PASS ${normalizedPaths} current physical rows equal approved capture; mutation rejected`);
