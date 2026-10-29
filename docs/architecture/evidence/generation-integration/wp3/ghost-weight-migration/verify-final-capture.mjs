import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';

const approved = 'target/ghost-weight-certificate/current/states/ghost_weight_certificate';
const migrated = 'target/ghost-weight-certificate/migration/final-diagnostic-state/ghost_weight_certificate';
const sha = (bytes) => crypto.createHash('sha256').update(bytes).digest('hex');
function selected(root) {
  const files = [];
  function walk(folder) {
    for (const entry of fs.readdirSync(folder, {withFileTypes: true})) {
      const location = path.join(folder, entry.name);
      if (entry.isDirectory()) walk(location);
      else if (location.endsWith('.dot') || location.endsWith('-graph-witness.tsv')
        || location.endsWith('-process-definition.json') || location.endsWith('-live-model.json')
        || location.endsWith('-live-parameters.json') || location.endsWith('-live-settings.json')) {
        files.push(path.relative(root, location));
      }
    }
  }
  walk(root);
  return files.sort();
}
const oldFiles = selected(approved);
const newFiles = selected(migrated);
assert.deepEqual(newFiles, oldFiles, 'final binary capture file membership differs');
const counts = {dots: 0, graphWitnesses: 0, processDefinitions: 0, live: 0};
const oldHashLines = [];
const newHashLines = [];
let normalizedSettings = 0;
for (const relative of oldFiles) {
  const oldBytes = fs.readFileSync(path.join(approved, relative));
  const newBytes = fs.readFileSync(path.join(migrated, relative));
  if (relative.endsWith('-live-settings.json')) {
    const normalize = (bytes) => {
      const settings = JSON.parse(bytes);
      delete settings.cli_settings.state.folder;
      return settings;
    };
    const oldSettings = normalize(oldBytes);
    const newSettings = normalize(newBytes);
    assert.deepEqual(newSettings, oldSettings, `live settings differ ${relative}`);
    oldHashLines.push(`${relative}\0${sha(JSON.stringify(oldSettings))}\n`);
    newHashLines.push(`${relative}\0${sha(JSON.stringify(newSettings))}\n`);
    normalizedSettings++;
  } else {
    assert.equal(sha(newBytes), sha(oldBytes), `captured bytes differ ${relative}`);
    oldHashLines.push(`${relative}\0${sha(oldBytes)}\n`);
    newHashLines.push(`${relative}\0${sha(newBytes)}\n`);
  }
  if (relative.endsWith('.dot')) counts.dots++;
  else if (relative.endsWith('-graph-witness.tsv')) counts.graphWitnesses++;
  else if (relative.endsWith('-process-definition.json')) counts.processDefinitions++;
  else counts.live++;
}
assert.deepEqual(counts, {dots: 4239, graphWitnesses: 15, processDefinitions: 15, live: 63});
assert.equal(normalizedSettings, 21);
const targets = [
  ['epem_bbh/bbh_qcd2-graph-witness.tsv', '594 | -201 = -201'],
  ['slow_aa_ttx/qcd2-graph-witness.tsv', '874 | -186 = -186'],
  ['vacuum_nonfactorizable/l2-graph-witness.tsv', '5 | -23/12 = -23/12'],
  ['vacuum_nonfactorizable/l3-graph-witness.tsv', '30 | -13/12 = -13/12'],
  ['vacuum_nonfactorizable/l4-graph-witness.tsv', '200 | 457/48 = 457/48'],
];
for (const [file, expected] of targets) {
  const header = fs.readFileSync(path.join(migrated, file), 'utf8').split('\n')[0];
  assert.equal(header, expected, `${file} scalar sum`);
}
const result = {schema: 1, approvedRoot: approved, finalRoot: migrated,
  counts, normalizedSettings, totalSelectedFiles: oldFiles.length,
  oldByteTreeSha256IgnoringSettingsFolder: sha(oldHashLines.join('')),
  finalByteTreeSha256IgnoringSettingsFolder: sha(newHashLines.join('')),
  perTarget: Object.fromEntries(targets.map(([file, expected]) => [file, expected])),
  comparison: 'all DOT, TSV, process definitions, live model/card bytes exact; only state.folder normalized in 21 live settings captures'};
fs.writeFileSync('target/ghost-weight-certificate/migration/final-capture-comparison.json',
  `${JSON.stringify(result, null, 2)}\n`);
console.log(`final capture: ${oldFiles.length} selected files equal; 4239 DOTs; 15 witnesses; 15 definitions; 63 live inputs; five exact sums`);
