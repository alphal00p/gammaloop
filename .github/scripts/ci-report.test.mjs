import assert from 'node:assert/strict';
import { execFile } from 'node:child_process';
import { mkdtemp, readFile, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { promisify } from 'node:util';
import test from 'node:test';
import { comparePairs, createReport, parseLog, summarizeActions, summarizeSuite } from './ci-report.mjs';

const overlap = await readFile(new URL('./fixtures/ci-overlap.ndjson', import.meta.url), 'utf8');
const wrapper = await readFile(new URL('./fixtures/ci-cached-wrapper.ndjson', import.meta.url), 'utf8');
const sha = 'a'.repeat(40);
const suiteUrl = 'https://nix-ci.com/gh:example:repo/ci-test/' + sha;
const spec = { layout: 'main', variant: 'baseline', scenario: 'warm', sha, suiteUrl };
const attr = 'packages.x86_64-linux.nix-ci-check-gammaloop-nextest-core';

test('completed transfers union within workers; Cargo and test evidence stay distinct', () => {
  const result = parseLog(overlap, { type: 'test', status: 'success', attribute: attr });
  assert.equal(result.workerSeconds, 14);
  assert.equal(result.downloadCount, 2);
  assert.equal(result.downloadReportedBytes, 3 * 1024 ** 2);
  assert.equal(result.downloadActiveSeconds, 6);
  assert.equal(result.uploadActiveSeconds, 3);
  assert.equal(result.transferActiveSeconds, 7);
  assert.equal(result.uploadReportedBytes, null);
  assert.equal(result.cargoFinishedSeconds, 4);
  assert.equal(result.compilationMessages, 1);
  assert.equal(parseLog(overlap.replace('Compiling example', 'Checking example')).checkingMessages, 1);
  assert.equal(result.testExecution, 'executed');
  assert.deepEqual(result.testSummaries.map(summary => [summary.executed, summary.skipped]), [[2, 1], [1, 0]]);
  assert.deepEqual(result.tests.map(item => item.status), ['PASS', 'SKIP', 'PASS']);
  assert.equal(result.builds.length, 2);
});

test('wrapper substitution is unknown; actual result substitution establishes reuse', () => {
  const job = { type: 'test', status: 'success', attribute: attr };
  assert.equal(parseLog(wrapper, job).testExecution, 'unknown');
  const result = wrapper.replaceAll('nix-ci-check-gammaloop-nextest-core', 'gammaloop-nextest-core');
  assert.equal(parseLog(result, job).testExecution, 'reused');
  assert.equal(parseLog(wrapper, { ...job, status: 'cached' }).testExecution, 'reused');
  assert.throws(() => parseLog(''), /empty log/);
  assert.throws(() => parseLog('{"log_message":"partial"}'), /invalid log record/);
  assert.throws(() => parseLog(overlap + '{'), /invalid NDJSON/);
});

test('suite clocks exclude queued deployment and unrelated attempts, while duplicate work counts once per job', () => {
  const job = {
    ...parseLog(overlap, { type: 'test', attribute: attr }), type: 'test', attribute: attr,
    url: suiteUrl + '/test', status: 'success', logStatus: 'ok', checkConclusion: 'success',
    checkStartedAt: '2026-09-06T00:00:02Z', checkCompletedAt: '2026-09-06T00:00:20Z',
  };
  const jobs = [
    { type: 'config', status: 'success', checkStartedAt: '2026-09-06T00:00:00Z', checkCompletedAt: '2026-09-06T00:00:01Z' },
    job,
    { ...job, type: 'build', attribute: 'packages.x86_64-linux.crate-test-binaries-core', url: suiteUrl + '/producer',
      checkCompletedAt: '2026-09-06T00:00:30Z', testExecution: 'not-applicable' },
    { type: 'deploy', status: 'queued', checkStartedAt: '2026-09-06T00:00:00Z', checkCompletedAt: null },
    { ...job, type: 'build', attribute: 'packages.other', status: 'cancelled', url: suiteUrl + '/cancelled', checkCompletedAt: '2026-09-06T02:00:00Z', builds: [] },
  ];
  const summary = summarizeSuite(spec, { status: 'failure' }, [{}], jobs);
  assert.equal(summary.suiteSeconds, 30);
  assert.equal(summary.requiredSeconds, 20);
  const withDocs = summarizeSuite(spec, { status: 'success' }, [{}], [...jobs,
    { ...job, type: 'build', attribute: 'checks.x86_64-linux.gammaloop-doc', checkCompletedAt: '2026-09-06T00:01:00Z' }]);
  assert.equal(withDocs.requiredSeconds, 20);
  assert.equal(summary.finalAggregateTailSeconds, null);
  assert.equal(summary.successfulTimingSample, false);
  assert.equal(summary.repeatedDerivationCount, 1);
  assert.equal(summary.repeatedBuildWorkers, 1);
  assert.equal(summary.repeatedDerivations[0].workerCount, 2);
  assert.equal(summarizeSuite(spec, { status: 'cancelled' }, [{}], jobs).suiteSeconds, null);
  assert.equal(summarizeSuite(spec, { status: 'running' }, [{}], jobs).evidenceComplete, false);
});

test('successful retries preserve abandoned work as lower bounds and cannot produce savings', () => {
  const config = { type: 'config', status: 'success', checkStartedAt: '2026-09-06T00:00:00Z',
    checkCompletedAt: '2026-09-06T00:00:01Z' };
  const core = { ...parseLog(overlap, { type: 'test', status: 'success', attribute: attr }),
    type: 'test', attribute: attr, status: 'success', logStatus: 'ok', url: suiteUrl + '/retry',
    checkStartedAt: '2026-09-06T00:00:10Z', checkCompletedAt: '2026-09-06T00:00:30Z' };
  const abandoned = { ...core, status: 'abandoned', url: suiteUrl + '/abandoned',
    checkStartedAt: '2026-09-06T00:00:01Z', checkCompletedAt: '2026-09-06T00:01:00Z' };
  const baseline = summarizeSuite(spec, { status: 'success' }, [{}], [config, core]);
  const candidate = summarizeSuite({ ...spec, variant: 'candidate' }, { status: 'success' }, [{}], [config, abandoned, core]);
  assert.equal(candidate.requiredPassed, true);
  assert.equal(candidate.requiredSeconds, 30);
  assert.equal(candidate.suiteSeconds, 30);
  assert.equal(candidate.observedWorkers, 2);
  assert.equal(candidate.observedWorkerMinutes, 28 / 60);
  assert.equal(candidate.downloadReportedBytes, 6 * 1024 ** 2);
  assert.equal(candidate.interruptedJobs, 1);
  assert.equal(candidate.missingLogs, 0);
  assert.equal(candidate.observedResourceLowerBound, true);
  assert.equal(candidate.evidenceComplete, false);
  assert.equal(baseline.observedResourceLowerBound, false);
  assert.equal(summarizeSuite(spec, { status: 'failure' }, [{}],
    [config, { ...core, status: 'failure' }]).evidenceComplete, true);
  assert.equal(comparePairs([baseline, candidate])[0].workerMinutesChangePercent, null);
});

test('Actions attempts retain their own clock and cancelled work without completed timing', () => {
  const run = {
    id: 1, run_attempt: 2, name: 'Nix', head_sha: sha, html_url: 'https://github.com/example/repo/actions/runs/1',
    status: 'completed', conclusion: 'success', created_at: '2026-09-06T00:00:00Z',
    run_started_at: '2026-09-06T01:00:00Z', updated_at: '2026-09-06T01:00:31Z',
    jobs: [{ status: 'completed', conclusion: 'success', started_at: '2026-09-06T01:00:10Z', completed_at: '2026-09-06T01:00:30Z' }],
  };
  const summary = summarizeActions(run);
  assert.equal(summary.createdToUpdatedSeconds, 3631);
  assert.equal(summary.attemptToLastJobSeconds, 30);
  assert.equal(summary.observedWorkerMinutes, 20 / 60);
  const cancelled = summarizeActions({ ...run, conclusion: 'cancelled' });
  assert.equal(cancelled.createdToUpdatedSeconds, null);
  assert.equal(cancelled.attemptToLastJobSeconds, null);
  assert.equal(cancelled.observedWorkerMinutes, 20 / 60);
});

test('paired deltas require complete successful evidence and never infer coverage equivalence', () => {
  const baseline = { ...spec, pair: '1', evidenceComplete: true, successfulTimingSample: true, requiredPassed: true,
    observedWorkerMinutes: 10, intermediateDownloadReportedBytes: 100, suiteStartedAt: '2026-09-06T00:00:00Z',
    jobs: [{ type: 'test', attribute: attr, checkStartedAt: '2026-09-06T00:00:05Z', checkCompletedAt: '2026-09-06T00:00:20Z', checkSeconds: 15 }] };
  const candidate = { ...baseline, variant: 'candidate', observedWorkerMinutes: 5, intermediateDownloadReportedBytes: 20,
    jobs: [{ ...baseline.jobs[0], checkStartedAt: '2026-09-06T00:00:15Z', checkCompletedAt: '2026-09-06T00:00:30Z' }] };
  const [pair] = comparePairs([baseline, candidate]);
  assert.equal(pair.comparable, true);
  assert.equal(pair.workerMinutesChangePercent, -50);
  assert.equal(pair.intermediateRestoreChangePercent, -80);
  assert.equal(pair.coverageEquivalent, null);
  assert.equal(pair.groups[0].baselineCompletedSeconds, 20);
  assert.equal(pair.groups[0].candidateCheckSeconds, 15);
  assert.equal(pair.groups[0].completionChangeSeconds, 10);
  assert.equal(comparePairs([baseline, { ...candidate, evidenceComplete: false }])[0].workerMinutesChangePercent, null);
});

test('offline collector scopes exact URLs, emits all formats and returns a failing CLI status for missing logs', async t => {
  const dir = await mkdtemp(join(tmpdir(), 'ci-report-test-'));
  t.after(() => rm(dir, { recursive: true, force: true }));
  const runs = [
    { type: 'config', status: 'success', url: suiteUrl + '/config' },
    { type: 'test', status: 'success', attribute: attr, url: suiteUrl + '/test' },
    { type: 'build', status: 'success', attribute: 'packages.x86_64-linux.ci-test-inputs', url: suiteUrl + '/producer' },
    { type: 'deploy', status: 'queued', url: suiteUrl + '/deploy' },
  ];
  const checks = [{ check_runs: [
    { id: 1, details_url: runs[0].url, started_at: '2026-09-06T00:00:00Z', completed_at: '2026-09-06T00:00:01Z', conclusion: 'success' },
    { id: 2, details_url: runs[1].url, started_at: '2026-09-06T00:00:01Z', completed_at: '2026-09-06T00:00:20Z', conclusion: 'success' },
    { id: 4, details_url: runs[2].url, started_at: '2026-09-06T00:00:01Z', completed_at: '2026-09-06T00:00:15Z', conclusion: 'success' },
    { id: 3, details_url: suiteUrl + '/other-attempt', started_at: '2026-09-06T12:00:00Z', completed_at: '2026-09-06T13:00:00Z', conclusion: 'success' },
  ] }];
  const manifest = { repository: 'example/repo', runs: [{ ...spec,
    offline: { suite: 'suite.json', checks: 'checks.json', jobs: 'jobs.json', actions: 'actions.json' } }] };
  for (const [name, value] of Object.entries({
    'manifest.json': manifest, 'suite.json': { commit: sha, status: 'success', runs },
    'checks.json': checks, 'jobs.json': [{ url: runs[0].url, file: 'wrapper.ndjson' }, { url: runs[1].url, file: 'overlap.ndjson' }, { url: runs[2].url, file: 'producer.ndjson' }],
    'actions.json': [],
  })) await writeFile(join(dir, name), JSON.stringify(value));
  await writeFile(join(dir, 'wrapper.ndjson'), wrapper);
  await writeFile(join(dir, 'overlap.ndjson'), overlap);
  await writeFile(join(dir, 'producer.ndjson'), overlap);
  const report = await createReport(join(dir, 'manifest.json'), join(dir, 'out'));
  assert.equal(report.complete, true);
  assert.equal(report.suites[0].suiteSeconds, 20);
  assert.equal(report.suites[0].executedTestJobs, 1);
  assert.equal(report.suites[0].intermediateDownloadReportedBytes, 3 * 1024 ** 2);
  assert.equal(report.suites[0].jobs[1].checkId, 2);
  assert.equal(report.suites[0].jobs[1].completedFromSuiteSeconds, 20);
  for (const name of ['manifest.json', 'report.json', 'report.md', 'suites.csv', 'jobs.csv', 'actions.csv', 'transfers.csv', 'compilations.csv', 'pairs.csv'])
    assert.ok((await readFile(join(dir, 'out', name), 'utf8')).length);
  runs[2].status = 'abandoned';
  await writeFile(join(dir, 'suite.json'), JSON.stringify({ commit: sha, status: 'success', runs }));
  const recovered = await createReport(join(dir, 'manifest.json'), join(dir, 'recovered'));
  assert.equal(recovered.complete, false);
  assert.equal(recovered.suites[0].requiredPassed, true);
  assert.equal(recovered.suites[0].interruptedJobs, 1);
  assert.equal(recovered.suites[0].jobs[2].logStatus, 'ok');
  assert.equal(recovered.suites[0].jobs[2].downloadReportedBytes, 3 * 1024 ** 2);
  assert.equal(recovered.suites[0].jobs[2].completedFromSuiteSeconds, null);
  assert.equal(recovered.suites[0].jobs[2].postWorkerSeconds, null);
  assert.ok(recovered.errors.some(error => error.stage === 'interrupted-worker'));
  runs[2].status = 'success';
  await writeFile(join(dir, 'suite.json'), JSON.stringify({ commit: sha, status: 'success', runs }));
  const incident = { jobUrl: runs[2].url, observedAt: '2026-09-06T00:00:05Z', evidenceFile: 'incident.json' };
  manifest.runs[0].observedInterruptions = [incident];
  await writeFile(join(dir, 'old-log.ndjson'), overlap);
  await writeFile(join(dir, 'incident.json'), JSON.stringify({ jobUrl: incident.jobUrl, status: 'abandoned',
    file: 'old-log.ndjson', workerSeconds: 14, downloadReportedBytes: 3 * 1024 ** 2 }));
  await writeFile(join(dir, 'manifest.json'), JSON.stringify(manifest));
  const replaced = await createReport(join(dir, 'manifest.json'), join(dir, 'replaced'));
  assert.equal(replaced.complete, false);
  assert.equal(replaced.suites[0].interruptedJobs, 1);
  assert.equal(replaced.suites[0].evidenceComplete, false);
  assert.equal(comparePairs([report.suites[0], { ...replaced.suites[0], variant: 'candidate' }])[0].workerMinutesChangePercent, null);
  assert.equal(replaced.suites[0].observedResourceLowerBound, true);
  assert.equal(replaced.suites[0].observedWorkerMinutes, report.suites[0].observedWorkerMinutes);
  assert.equal(replaced.suites[0].downloadReportedBytes, report.suites[0].downloadReportedBytes);
  assert.equal(replaced.suites[0].observedInterruptions[0].evidence.workerSeconds, 14);
  assert.ok((await readFile(join(dir, 'replaced/report.md'), 'utf8')).includes('Prior interruption:'));
  await rm(join(dir, 'old-log.ndjson'));
  await rm(join(dir, 'incident.json'));
  incident.evidenceFile = 'replaced/raw/1/interruptions/1.json';
  await writeFile(join(dir, 'manifest.json'), JSON.stringify(manifest));
  const replayed = await createReport(join(dir, 'manifest.json'), join(dir, 'replayed'));
  assert.equal(replayed.complete, false);
  assert.equal(replayed.errors.some(error => error.stage === 'interruption-evidence'), false);
  assert.equal(await readFile(join(dir, 'replayed/raw/1/interruptions/1.ndjson'), 'utf8'), overlap);
  incident.jobUrl = suiteUrl + '/other-suite/job';
  await writeFile(join(dir, 'manifest.json'), JSON.stringify(manifest));
  await assert.rejects(createReport(join(dir, 'manifest.json'), join(dir, 'invalid')), /exact-suite jobUrl/);
  delete manifest.runs[0].observedInterruptions;
  await writeFile(join(dir, 'manifest.json'), JSON.stringify(manifest));
  await rm(join(dir, 'overlap.ndjson'));
  const partial = await createReport(join(dir, 'manifest.json'), join(dir, 'partial'));
  assert.equal(partial.complete, false);
  assert.equal(partial.suites[0].missingLogs, 1);
  assert.equal(partial.suites[0].unknownTestJobs, 1);
  await assert.rejects(promisify(execFile)(process.execPath,
    [fileURLToPath(new URL('./ci-report.mjs', import.meta.url)), join(dir, 'manifest.json'), join(dir, 'cli')]),
    error => error.code === 2);
});

test('omitting an explicitly required group cannot look like a worker-time saving', () => {
  const doctest = 'packages.x86_64-linux.nix-ci-check-gammaloop-doctest';
  const expectation = { ...spec, requiredAttributes: [attr, doctest] };
  const config = { type: 'config', status: 'success', checkStartedAt: '2026-09-06T00:00:00Z',
    checkCompletedAt: '2026-09-06T00:00:01Z' };
  const core = { type: 'test', attribute: attr, status: 'success', checkConclusion: 'success', logStatus: 'ok',
    checkStartedAt: '2026-09-06T00:00:01Z', checkCompletedAt: '2026-09-06T00:00:21Z', workerSeconds: 20 };
  const docs = { ...core, attribute: doctest, workerSeconds: 21 };
  const baseline = summarizeSuite(expectation, { status: 'success' }, [{}], [config, core, docs]);
  const candidate = summarizeSuite({ ...expectation, variant: 'candidate' }, { status: 'success' }, [{}], [config, core]);
  assert.equal(baseline.evidenceComplete, true);
  assert.equal(candidate.requiredPassed, false);
  assert.equal(candidate.evidenceComplete, false);
  assert.deepEqual(candidate.requiredAttributes, [attr, doctest]);
  assert.deepEqual(candidate.observedRequiredAttributes, [attr]);
  assert.deepEqual(candidate.missingRequiredAttributes, [doctest]);
  const [pair] = comparePairs([baseline, candidate]);
  assert.equal(pair.comparable, false);
  assert.equal(pair.workerMinutesChangePercent, null);
});
