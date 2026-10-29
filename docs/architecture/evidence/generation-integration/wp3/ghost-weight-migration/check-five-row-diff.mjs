import fs from 'node:fs';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';

const path = 'tests/tests/test_feyngen.rs';
const old = execFileSync('jj', ['file', 'show', '-r', 'generation-ghost-weight-certificate-approved', path],
  {encoding: 'utf8', maxBuffer: 8 * 1024 * 1024});
const current = fs.readFileSync(path, 'utf8');
const before = old.split('\n');
const after = current.split('\n');
assert.equal(before.length, after.length, 'test source line count changed');
const source44 = JSON.parse(fs.readFileSync(
  'docs/architecture/evidence/generation-integration/wp3/ghost-weight-certificate/source44-rows.json'));
assert.equal(source44.cases.length, 5);
const changed = before.flatMap((line, index) => line === after[index] ? [] : [index + 1]);
const expectedLines = source44.cases.map((row) => row.current.line);
assert.deepEqual(changed, expectedLines, 'only five target assertion lines may change');
const sha = (text) => crypto.createHash('sha256').update(text).digest('hex');
const rows = source44.cases.map((row) => {
  const index = row.current.line - 1;
  const original = before[index];
  const replacement = after[index];
  const oldLiteral = `@"${row.source43.expected}"`;
  const newLiteral = `@"${row.source44.expected}"`;
  assert(original.includes(row.exact_process_fragment), `${row.id} original command`);
  assert(replacement.includes(row.exact_process_fragment), `${row.id} new command`);
  assert(original.includes(oldLiteral) && replacement.includes(newLiteral), `${row.id} literals`);
  assert.equal(original.replace(oldLiteral, newLiteral), replacement,
    `${row.id} contains nonliteral edit`);
  const owner = after.slice(0, index + 1).reverse().find((line) => /^\s*fn test_/.test(line));
  assert(owner?.includes(`${row.qualified_test.split('::')[1]}(`), `${row.id} module-qualified owner`);
  assert(after.slice(0, index + 1).some((line) => line === 'mod slow {'), `${row.id} slow module`);
  return {id: row.id, qualifiedTest: row.qualified_test, line: row.current.line,
    originalLiteral: row.source43.expected, migratedLiteral: row.source44.expected,
    command: row.exact_process_fragment,
    originalLineSha256: sha(original), migratedLineSha256: sha(replacement)};
});
const failingStart = after.findIndex((line) => line.includes('fn test_generate_aa_ttx_cross_section_slow_filter()'));
assert(failingStart > 0);
const failingEnd = after.findIndex((line, index) => index > failingStart && /^\s*fn test_/.test(line));
const failing = after.slice(failingStart, failingEnd > 0 ? failingEnd : undefined).join('\n');
assert(failing.includes('@"874 | -186 = -186"'), 'approved distinct failing-row literal moved');
const summary = {
  schema: 1,
  approvedCertificateCommit: '121ac233d5ab014f4fadefccce702041bbbdf487',
  source44Commit: source44.sources.source44,
  path,
  originalSourceSha256: sha(old), migratedSourceSha256: sha(current),
  lineCount: before.length, changedLineCount: changed.length, changedLines: changed,
  rows,
  negativeInventory: {
    allOtherLinesByteEqual: true,
    distinctFailingAaRowStillApproved: true,
    groupedSunsetAndOtherGroupedSnapshotLinesChanged: 0,
    productionSourceFilesChanged: 0,
    setupModelCommandFilterWorkerAndStateLinesChanged: 0,
    earlierAndOtherAssertionsChanged: 0,
  },
};
fs.writeFileSync('target/ghost-weight-certificate/migration/five-row-diff.json',
  `${JSON.stringify(summary, null, 2)}\n`);
console.log(`exact five-row diff: ${changed.join(', ')}; all other ${before.length - changed.length} source lines byte equal`);
