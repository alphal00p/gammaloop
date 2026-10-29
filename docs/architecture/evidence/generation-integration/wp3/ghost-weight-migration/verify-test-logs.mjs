import fs from 'node:fs';
import crypto from 'node:crypto';
import assert from 'node:assert/strict';

const root = 'target/ghost-weight-certificate/migration';
const sha = (file) => crypto.createHash('sha256').update(fs.readFileSync(file)).digest('hex');
const log = (name) => {
  const file = `${root}/${name}.log`;
  return {file, sha256: sha(file), text: fs.readFileSync(file, 'utf8')};
};
const targets = [
  {id: 'epem_bbh_qcd2', old: '594 | -261 = -261', current: '594 | -201 = -201',
    oldTest: 'generate_epem_a_qqh-4', later: 'generate_epem_a_qqh-7', names: ['old-epem', 'new-epem']},
  {id: 'slow_aa_ttx_qcd2', old: '874 | -266 = -266', current: '874 | -186 = -186',
    oldTest: 'generate_aa_ttx_cross_section-3', later: 'generate_aa_ttx_cross_section-4', names: ['old-aa', 'new-aa']},
];
const cases = {};
for (const item of targets) {
  const old = log(item.names[0]);
  const fresh = log(item.names[1]);
  assert(old.text.includes(`-${item.old}`) && old.text.includes(`+${item.current}`), `${item.id} old-fail values`);
  assert(old.text.includes(`snapshot assertion for '${item.oldTest}' failed`), `${item.id} old row not reached`);
  assert(/test result: FAILED\. 0 passed; 1 failed/.test(old.text), `${item.id} old test count`);
  assert(fresh.text.includes(`snapshot assertion for '${item.later}' failed`), `${item.id} later boundary`);
  assert(!fresh.text.includes(`snapshot assertion for '${item.oldTest}' failed`), `${item.id} corrected row still failing`);
  assert(/test result: FAILED\. 0 passed; 1 failed/.test(fresh.text), `${item.id} complete-test outcome`);
  cases[item.id] = {old: {path: old.file, sha256: old.sha256, exitCode: 101, target: 'failed',
    testCount: {passed: 0, failed: 1}},
    new: {path: fresh.file, sha256: fresh.sha256, exitCode: 101, target: 'passed; subsequent unchanged grouped row failed',
      laterFailure: item.later, testCount: {passed: 0, failed: 1}}};
}
for (const loop of [2, 3, 4]) {
  const old = log(`vacuum-old-l${loop}`);
  const fresh = log(`vacuum-new-l${loop}`);
  const expected = {2: ['5 | -11/12 = -11/12', '5 | -23/12 = -23/12'],
    3: ['30 | -17/12 = -17/12', '30 | -13/12 = -13/12'],
    4: ['200 | -199/48 = -199/48', '200 | 457/48 = 457/48']}[loop];
  assert(old.text.includes(`-${expected[0]}`) && old.text.includes(`+${expected[1]}`), `L${loop} old failure values`);
  assert(old.text.includes(`snapshot assertion for 'scoped_vacuum_l${loop}' failed`), `L${loop} old assertion owner`);
  assert(/test result: FAILED\. 0 passed; 1 failed/.test(old.text), `L${loop} old test count`);
  assert(/test result: ok\. 1 passed; 0 failed/.test(fresh.text), `L${loop} new test count`);
  for (const text of [old.text, fresh.text]) {
    assert(text.includes('Running command import model sm.json'), `L${loop} model import`);
    assert(text.includes('Running command set global kv global.n_cores.feyngen=16'), `L${loop} original worker choice`);
    assert(text.includes(`Running command generate amp {} > {} | g ghg t u d QED==0 [{${loop}}] --filter-zero-color --numerator-grouping no_grouping --max-n-bridges 0 --number-of-factorized-loop-subtopologies 1 1 --clear-existing-processes --only-diagrams`), `L${loop} exact command`);
  }
  cases[`vacuum_l${loop}`] = {old: {path: old.file, sha256: old.sha256, exitCode: 101,
    target: 'failed at old literal', testCount: {passed: 0, failed: 1}},
  new: {path: fresh.file, sha256: fresh.sha256, exitCode: 0,
    target: 'passed', testCount: {passed: 1, failed: 0}}};
}
const fullVacuum = log('new-original-vacuum');
assert(fullVacuum.text.includes("snapshot assertion for 'vacuum_amplitude_generation-2' failed in line 2391"));
assert(/test result: FAILED\. 0 passed; 1 failed/.test(fullVacuum.text));
const diagnostic = log('final-diagnostic');
assert(/test result: ok\. 1 passed; 0 failed/.test(diagnostic.text));
const report = {schema: 1, cases, originalVacuumFunction: {
  path: fullVacuum.file, sha256: fullVacuum.sha256, exitCode: 101,
  firstFailure: 'unchanged grouped two-loop assertion at line 2391, before nonfactorizable rows',
  testCount: {passed: 0, failed: 1}},
  finalDiagnostic: {path: diagnostic.file, sha256: diagnostic.sha256, exitCode: 0,
    testCount: {passed: 1, failed: 0}},
  scope: 'Two exact named cross-section tests; six independent scoped vacuum assertion tests; one original vacuum-prefix control; one ignored physical witness capture.'};
fs.writeFileSync(`${root}/test-results.json`, `${JSON.stringify(report, null, 2)}\n`);
console.log('test logs: two old cross target failures; two new cross targets passed before later grouped failures; three old vacuum failures; three new vacuum passes; original grouped prefix failure; final diagnostic pass');
