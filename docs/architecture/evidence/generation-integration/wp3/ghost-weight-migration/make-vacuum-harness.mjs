import fs from 'node:fs';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';

const evidence = 'target/ghost-weight-certificate/migration';
const sourcePath = 'tests/tests/test_feyngen.rs';
const approvedRef = 'generation-ghost-weight-certificate-approved';
const oldSource = execFileSync('jj', ['file', 'show', '-r', approvedRef, sourcePath],
  {encoding: 'utf8', maxBuffer: 8 * 1024 * 1024});
const newSource = fs.readFileSync(sourcePath, 'utf8');
const rows = JSON.parse(fs.readFileSync(
  'docs/architecture/evidence/generation-integration/wp3/ghost-weight-certificate/source44-rows.json')).cases
  .filter((row) => row.id.startsWith('vacuum_nonfactorizable_'));
assert.equal(rows.length, 3);
const sha = (text) => crypto.createHash('sha256').update(text).digest('hex');
const between = (source, first, after) => {
  const start = source.indexOf(first);
  const end = source.indexOf(after, start);
  assert(start >= 0 && end > start, `missing helper source range ${first}`);
  return source.slice(start, end).trimEnd();
};
const helpers = (source) => [
  between(source, 'fn evaluate_sign_origin(', '\nfn setup_tth_select_cli('),
  between(source, 'fn split_before_flags(', '\n// Diagnostic-only graph witness.'),
];
assert.deepEqual(helpers(oldSource), helpers(newSource), 'generation helper source changed');
const imports = `use color_eyre::Result;
use feynkit_graph::expressions::evaluate_overall_factor;
use gammaloop_integration_tests::{CLIState, get_test_cli, get_tests_workspace_path};
use gammalooprs::processes::ProcessCollection;
use insta::assert_snapshot;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::printer::CanonicalOrderingSettings;
use symbolica::{function, symbol};
use tracing::debug;
`;
function extract(source, route) {
  const anchor = source.indexOf('fn test_vacuum_amplitude_generation()');
  assert(anchor > source.indexOf('mod slow {'), 'vacuum test not in slow module');
  const block = source.slice(anchor, source.indexOf('\n    #[test]', anchor));
  const lines = block.split('\n');
  const setup = lines.filter((line) => line.includes('let mut cli = get_test_cli(')
    || line.includes('cli.run_command("import model sm.json")?;')
    || line.includes('cli.run_command("set global kv global.n_cores.feyngen=16")?;'));
  assert.equal(setup.length, 3, 'vacuum setup count');
  const assertions = rows.map((row) => {
    const found = lines.filter((line) => line.includes(row.exact_process_fragment)
      && line.includes('assert_snapshot!(feyngen_str(&mut cli, "amp"'));
    assert.equal(found.length, 1, `${route}/${row.id} assertion owner`);
    const expected = route === 'old' ? row.source43.expected : row.source44.expected;
    assert(found[0].includes(`@"${expected}"`), `${route}/${row.id} expected literal`);
    return {id: row.id, sourceLine: found[0], expected,
      command: row.exact_process_fragment, sourceLineSha256: sha(found[0])};
  });
  return {setup, assertions};
}
const old = extract(oldSource, 'old');
const fresh = extract(newSource, 'new');
assert.deepEqual(old.setup, fresh.setup, 'vacuum model/import/worker setup changed');
for (const [index, a] of old.assertions.entries()) {
  const b = fresh.assertions[index];
  assert.equal(a.command, b.command);
  assert.equal(a.sourceLine.replace(`@"${a.expected}"`, '@"VALUE"'),
    b.sourceLine.replace(`@"${b.expected}"`, '@"VALUE"'),
    'vacuum assertion has nonliteral change');
}
function harness(sourceRoute, content) {
  const tests = content.assertions.map((row, index) => `#[test]
#[rustfmt::skip]
fn scoped_vacuum_l${index + 2}() -> Result<()> {
${content.setup.join('\n')}
${row.sourceLine}
    Ok(())
}`).join('\n\n');
  const text = `${imports}\n${helpers(sourceRoute).join('\n\n')}\n\n${tests}\n`;
  return text;
}
const oldHarness = harness(oldSource, old);
const newHarness = harness(newSource, fresh);
const oldPath = `${evidence}/vacuum-old-assertions.rs`;
const newPath = `${evidence}/vacuum-new-assertions.rs`;
fs.writeFileSync(oldPath, oldHarness);
fs.writeFileSync(newPath, newHarness);
const report = {schema: 1, owner: sourcePath, approvedRef,
  oldOwnerSha256: sha(oldSource), newOwnerSha256: sha(newSource),
  helperRanges: helpers(oldSource).map((text) => ({sha256: sha(text), bytes: Buffer.byteLength(text)})),
  old: {...old, source: {path: oldPath, sha256: sha(oldHarness)}},
  new: {...fresh, source: {path: newPath, sha256: sha(newHarness)}},
  adaptation: 'Standalone rustc --test harness copies unchanged owner scalar helpers, exact original model/import/16-worker setup and each original assertion line into separate tests; earlier unrelated grouped assertions are not executed.'};
fs.writeFileSync(`${evidence}/vacuum-harness-source.json`, `${JSON.stringify(report, null, 2)}\n`);
console.log(`vacuum harness source: three independently runnable old/new assertion cells; helpers and setup exact`);
