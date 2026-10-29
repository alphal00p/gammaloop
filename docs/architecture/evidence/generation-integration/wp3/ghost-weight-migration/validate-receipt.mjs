import fs from 'node:fs';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';

const evidence = 'docs/architecture/evidence/generation-integration/wp3/ghost-weight-migration';
const guard = execFileSync('node', [`${evidence}/validate-proof.mjs`],
  {encoding: 'utf8', maxBuffer: 8 * 1024 * 1024}).trim();
const receipt = JSON.parse(fs.readFileSync(`${evidence}/receipt.json`));
assert.equal(receipt.status, 'prepared five-literal migration; NEW independent Astra migration/test review required');
assert.equal(receipt.literalMigration.rows.length, 5);
assert.equal(receipt.testExecution.oldTargetFailures, 5);
assert.equal(receipt.testExecution.newTargetPasses, 5);
assert.equal(receipt.testExecution.completeNamedCrossSectionTestPasses, 0);
assert.equal(receipt.testExecution.originalVacuumWholeTestPasses, 0);
const sha = (location) => crypto.createHash('sha256').update(fs.readFileSync(location)).digest('hex');
let bindings = 0;
function visit(value) {
  if (Array.isArray(value)) return value.forEach(visit);
  if (!value || typeof value !== 'object') return;
  if ('path' in value && 'sha256' in value) {
    assert.equal(sha(value.path), value.sha256, `changed migration binding: ${value.path}`);
    bindings++;
  }
  Object.values(value).forEach(visit);
}
visit(receipt);
const docs = JSON.parse(fs.readFileSync('target/ghost-weight-certificate/migration/migration-doc-audit.json'));
assert.deepEqual(docs, {errors: [], table_count: 1, th_count: 3});
const diff = JSON.parse(fs.readFileSync('target/ghost-weight-certificate/migration/five-row-diff.json'));
assert.deepEqual(diff.changedLines, [2334, 2366, 2402, 2404, 2407]);
const note = fs.readFileSync('docs/architecture/generation-ghost-weight-migration.typ', 'utf8');
const registry = fs.readFileSync('docs/developers.toml', 'utf8');
const title = note.split('\n')[0].slice(2);
assert.equal(title, 'Five certified ungrouped source44 snapshot migrations');
assert(registry.includes(`title = "${title}"`));
assert(registry.includes('source = "docs/architecture/generation-ghost-weight-migration.typ"'));
console.log(`${guard}; ${bindings} migration bindings; five literals and scoped docs valid`);
