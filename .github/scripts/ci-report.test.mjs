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

test('transfer timeouts retain numeric retry evidence without inventing completed transfers', () => {
  const result = parseLog(JSON.stringify({ utc_time: '2026-09-08T15:00:00Z', relative_nanoseconds: 0,
    log_message: "warning: unable to upload 'https://cache.example/nar?token=private-sentinel': Timeout was reached (28) Operation too slow. Less than 1 bytes/sec transferred the last 300 seconds; retrying in 341 ms (attempt 1/5)\n"
      + "error: unable to download 'https://cache.example/nar': Operation too slow. Less than 2 bytes/sec transferred the last 60 seconds",
  }));
  assert.deepEqual(result.transferTimeouts, [
    { direction: 'upload', time: Date.parse('2026-09-08T15:00:00Z'), minimumBytesPerSecond: 1, windowSeconds: 300,
      retryDelayMilliseconds: 341, attempt: 1, maximumAttempts: 5 },
    { direction: 'download', time: Date.parse('2026-09-08T15:00:00Z'), minimumBytesPerSecond: 2, windowSeconds: 60,
      retryDelayMilliseconds: null, attempt: null, maximumAttempts: null },
  ]);
  assert.equal(result.uploadCount, 0);
  assert.equal(result.downloadCount, 0);
  assert.equal(result.transferActiveSeconds, 0);
  assert.equal(result.transferTimeoutCount, 2);
  assert.equal(JSON.stringify(result).includes('private-sentinel'), false);
  assert.equal(summarizeSuite(spec, { status: 'running' }, [], [{ ...result, logStatus: 'ok' }]).transferTimeoutCount, 2);
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

test('HTTP/2 and other failed cache transfers are sanitized incidents, not completed bytes', () => {
  const records = Array.from({ length: 15 }, (_, i) => JSON.stringify({ utc_time: `2026-09-08T15:53:${String(i).padStart(2, '0')}Z`, relative_nanoseconds: i * 1e9,
    log_message: "error: unable to download 'https://cache.example/nar?token=private-sentinel': HTTP error 200 (curl error: Stream error in the HTTP/2 framing layer)" }));
  records.push(JSON.stringify({ utc_time: '2026-09-08T15:53:15Z', relative_nanoseconds: 15e9,
    log_message: "warning: unable to upload 'https://cache.example/private-sentinel': HTTP error 503; retrying in 341 ms (attempt 2/5)" }));
  const parsed = parseLog(records.join('\n'));
  assert.equal(parsed.transferErrorCount, 16);
  assert.equal(parsed.incidents.filter(incident => incident.kind === 'http2-transfer-error').length, 15);
  assert.deepEqual(parsed.transferErrors.at(-1), { kind: 'cache-transfer-error', direction: 'upload', time: Date.parse('2026-09-08T15:53:15Z'),
    httpStatus: 503, retryDelayMilliseconds: 341, attempt: 2, maximumAttempts: 5 });
  assert.equal(parsed.downloadReportedBytes, 0);
  assert.equal(parsed.uploadCount, 0);
  assert.equal(JSON.stringify(parsed).includes('private-sentinel'), false);
  assert.ok(parsed.incidents.every(incident => incident.cause === 'unknown'));
});

test('invalid clocks, same-URL attempts and failed suites cannot become accepted comparisons', () => {
  const config = { type: 'config', status: 'success', checkStartedAt: '2026-09-06T00:00:00Z', checkCompletedAt: '2026-09-06T00:00:01Z' };
  const core = { type: 'test', attribute: attr, url: suiteUrl + '/test', status: 'success', logStatus: 'ok', checkStartedAt: '2026-09-06T00:00:01Z',
    checkCompletedAt: '2026-09-06T00:00:21Z', ...parseLog(overlap) };
  const baseline = summarizeSuite(spec, { status: 'success' }, [{}], [config, core]);
  const backwards = parseLog(overlap.replace('2026-09-06 00:00:05', '2026-09-05 00:00:05'));
  assert.equal(backwards.workerSeconds, null);
  assert.equal(backwards.incidents[0].kind, 'invalid-clock');
  const candidate = summarizeSuite({ ...spec, variant: 'candidate' }, { status: 'success' }, [{}], [config, { ...core, ...backwards }]);
  assert.equal(candidate.evidenceComplete, false);
  assert.equal(comparePairs([baseline, candidate])[0].comparable, false);
  for (const change of [{ checkStartedAt: null }, { checkCompletedAt: '2026-09-05T00:00:00Z' },
    { checkAttempts: [{ id: 1 }, { id: 2 }] }]) {
    const result = summarizeSuite(spec, { status: 'success' }, [{}], [config, { ...core, ...change }]);
    assert.equal(result.evidenceComplete, false);
    assert.ok(result.incidents.length > 0);
  }
  const failed = summarizeSuite({ ...spec, variant: 'candidate' }, { status: 'failure' }, [{}], [config, { ...core, status: 'failure' }]);
  assert.equal(comparePairs([baseline, failed])[0].workerMinutesChangePercent, null);
});

test('final-success tail and submission clocks retain their distinct endpoints', () => {
  const config = { type: 'config', status: 'success', checkStartedAt: '2026-09-06T00:00:10Z', checkCompletedAt: '2026-09-06T00:00:11Z' };
  const core = { type: 'test', attribute: attr, url: suiteUrl + '/test', status: 'success', logStatus: 'ok',
    checkStartedAt: '2026-09-06T00:00:11Z', checkCompletedAt: '2026-09-06T00:00:30Z' };
  const deploy = { type: 'deploy', url: suiteUrl + '/deploy', status: 'success', checkStartedAt: '2026-09-06T00:01:30Z', checkCompletedAt: '2026-09-06T00:01:31Z' };
  const submitted = { ...spec, submittedAt: '2026-09-06T00:00:00Z', submissionCompletedAt: '2026-09-06T00:00:04Z' };
  const result = summarizeSuite(submitted, { status: 'success' }, [{}], [config, core, deploy]);
  assert.equal(result.requiredSeconds, 20);
  assert.equal(result.submissionToRequiredSeconds, 30);
  assert.equal(result.submissionDurationSeconds, 4);
  assert.deepEqual(result.submissionToRequiredBoundsSeconds, { minimum: 26, maximum: 30 });
  assert.equal(result.incidents.find(incident => incident.kind === 'final-success-tail').seconds, 61);
  const atThreshold = summarizeSuite(spec, { status: 'success' }, [{}], [config, core, { ...deploy, checkCompletedAt: '2026-09-06T00:01:30Z' }]);
  assert.equal(atThreshold.incidents.some(incident => incident.kind === 'final-success-tail'), false);
  const impossible = summarizeSuite({ ...submitted, submittedAt: '2026-09-06T00:02:00Z' }, { status: 'success' }, [{}], [config, core, deploy]);
  assert.equal(impossible.evidenceComplete, false);
  assert.equal(impossible.submissionToRequiredSeconds, null);
});

test('queue incidents require same-snapshot declared prerequisites; recovery observations preserve uncertainty', () => {
  const producer = { type: 'build', attribute: 'checks.x86_64-linux.archive', status: 'success', url: suiteUrl + '/producer' };
  const waiting = { type: 'test', attribute: attr, status: 'queued', url: suiteUrl + '/waiting' };
  const evidence = { ...spec, dependencies: { [attr]: [producer.attribute] }, observations: [
    { at: '2026-09-06T00:00:00Z', evidenceFile: 'first.json', suite: { runs: [producer, waiting] } },
    { at: '2026-09-06T00:01:00Z', evidenceFile: 'second.json', suite: { runs: [producer, waiting] } },
  ] };
  const result = summarizeSuite(evidence, { status: 'running' }, [], [producer, waiting]);
  const incident = result.incidents.find(incident => incident.kind === 'queued-after-declared-prerequisites');
  assert.equal(incident.observations.length, 2);
  assert.equal(incident.seconds, undefined);
  assert.equal(incident.cause, 'unknown');
  assert.equal(summarizeSuite({ ...evidence, dependencies: undefined }, { status: 'running' }, [], [producer, waiting]).incidents.some(incident => incident.kind === 'queued-after-declared-prerequisites'), false);
  for (const runs of [[waiting], [{ ...producer, status: 'running' }, waiting], [producer, producer, waiting], [producer, { ...waiting, status: 'running' }]]) {
    const unproven = summarizeSuite({ ...evidence, observations: [{ ...evidence.observations[0], suite: { runs } }] }, { status: 'running' }, [], runs);
    assert.equal(unproven.incidents.some(incident => incident.kind === 'queued-after-declared-prerequisites'), false);
  }
  const recovery = summarizeSuite({ ...evidence, observations: [
    { ...evidence.observations[0], suite: { runs: [producer, { ...waiting, status: 'abandoned' }] } },
    { ...evidence.observations[1], suite: { runs: [producer, { ...waiting, status: 'running' }] } },
  ] }, { status: 'success' }, [], [producer, { ...waiting, status: 'success' }]);
  assert.equal(recovery.interruptedJobs, 1);
  assert.equal(recovery.evidenceComplete, false);
  assert.ok(recovery.incidents.some(incident => incident.kind === 'repeated-attempt' && incident.previousStatus === 'abandoned'));
});

test('saved observations replay with incident exports and reject another revision', async t => {
  const dir = await mkdtemp(join(tmpdir(), 'ci-report-observations-'));
  t.after(() => rm(dir, { recursive: true, force: true }));
  const runs = [
    { type: 'config', status: 'success', url: suiteUrl + '/config' },
    { type: 'build', status: 'cached', attribute: 'checks.x86_64-linux.archive', url: suiteUrl + '/archive' },
    { type: 'test', status: 'success', attribute: attr, url: suiteUrl + '/test' },
  ];
  const checks = runs.map((job, i) => ({ id: i + 1, details_url: job.url, started_at: '2026-09-06T00:00:00Z',
    completed_at: i === 2 ? '2026-09-06T00:00:20Z' : '2026-09-06T00:00:01Z', conclusion: 'success' }));
  const manifest = { repository: 'example/repo', runs: [{ ...spec, submittedAt: '2026-09-05T23:59:59Z',
    observationFiles: ['snapshot.json'], dependencySnapshot: { sha, file: 'config.json' },
    offline: { suite: 'suite.json', checks: 'checks.json', jobs: 'jobs.json', actions: 'actions.json' } }] };
  const snapshot = { at: '2026-09-06T00:00:05Z', suites: [{ spec, suite: { commit: sha, status: 'running',
    runs: runs.map(job => job.type === 'test' ? { ...job, status: 'queued' } : job) } }] };
  for (const [name, value] of Object.entries({ 'manifest.json': manifest, 'snapshot.json': snapshot,
    'suite.json': { commit: sha, status: 'success', runs }, 'checks.json': checks, 'actions.json': [],
    'config.json': { dependencies: { [attr]: [runs[1].attribute] } },
    'jobs.json': [{ url: runs[0].url, file: 'log.ndjson' }, { url: runs[2].url, file: 'log.ndjson' }],
  })) await writeFile(join(dir, name), JSON.stringify(value));
  await writeFile(join(dir, 'log.ndjson'), overlap);
  const result = await createReport(join(dir, 'manifest.json'), join(dir, 'first'));
  assert.equal(result.complete, true);
  assert.equal(result.suites[0].submissionToRequiredSeconds, 21);
  const queue = result.incidents.find(incident => incident.kind === 'queued-after-declared-prerequisites');
  assert.equal(queue.observedAt, new Date(snapshot.at).toISOString());
  for (const file of ['incidents.json', 'incidents.csv', 'report.md'])
    assert.ok((await readFile(join(dir, 'first', file), 'utf8')).includes('queued-after-declared-prerequisites'));
  manifest.runs[0].observationFiles = ['first/raw/1/observations/1.json'];
  manifest.runs[0].dependencySnapshot.file = 'first/raw/1/dependency-snapshot.json';
  await rm(join(dir, 'snapshot.json'));
  await rm(join(dir, 'config.json'));
  await writeFile(join(dir, 'manifest.json'), JSON.stringify(manifest));
  const replay = await createReport(join(dir, 'manifest.json'), join(dir, 'replay'));
  assert.equal(replay.complete, true);
  assert.equal(replay.incidents.filter(incident => incident.kind === queue.kind).length, 1);
  await writeFile(join(dir, 'log.ndjson'), '{"log_message":"clock absent"}');
  const missingClock = await createReport(join(dir, 'manifest.json'), join(dir, 'missing-clock'));
  assert.equal(missingClock.complete, false);
  assert.ok(missingClock.incidents.some(incident => incident.kind === 'missing-clock'));
  assert.ok((await readFile(join(dir, 'missing-clock/raw/1/logs/1.ndjson'), 'utf8')).includes('clock absent'));
  manifest.runs[0].dependencySnapshot.sha = 'b'.repeat(40);
  await writeFile(join(dir, 'manifest.json'), JSON.stringify(manifest));
  await assert.rejects(createReport(join(dir, 'manifest.json'), join(dir, 'wrong-sha')), /same SHA/);
});

test('normal build/test phases are not repeat attempts and only build prerequisites establish readiness', () => {
  const producer = { type: 'build', attribute: 'checks.x86_64-linux.archive', status: 'success', url: suiteUrl + '/producer' };
  const runner = { ...producer, type: 'test', url: suiteUrl + '/runner' };
  const waiting = { type: 'test', attribute: attr, status: 'queued', url: suiteUrl + '/waiting' };
  const result = summarizeSuite(spec, { status: 'running' }, [], [producer, runner]);
  assert.equal(result.incidents.some(incident => incident.kind === 'repeated-attempt'), false);
  const repeated = summarizeSuite(spec, { status: 'running' }, [], [producer, runner, { ...producer, url: suiteUrl + '/producer-repeat' }]);
  assert.equal(repeated.incidents.filter(incident => incident.kind === 'repeated-attempt').length, 1);
  assert.equal(repeated.incidents.find(incident => incident.kind === 'repeated-attempt').attemptCount, 2);
  const observations = [{ at: '2026-09-06T00:00:01Z', evidenceFile: 'snapshot.json', suite: { runs: [producer, { ...runner, status: 'started' }, waiting] } }];
  const queued = summarizeSuite({ ...spec, observations, dependencies: { [attr]: [producer.attribute] } }, { status: 'running' }, [], [producer, runner, waiting]);
  assert.equal(queued.incidents.filter(incident => incident.kind === 'queued-after-declared-prerequisites').length, 1);
  const recovered = summarizeSuite({ ...spec, observations: [
    { ...observations[0], suite: { runs: [{ ...producer, status: 'failed' }] } },
    { ...observations[0], at: '2026-09-06T00:00:02Z', suite: { runs: [{ ...producer, status: 'started' }] } },
  ] }, { status: 'running' }, [], [producer]);
  assert.ok(recovered.incidents.some(incident => incident.kind === 'repeated-attempt' && incident.status === 'started'));
});

test('requested missing snapshot evidence invalidates suite pair eligibility, not only report status', async t => {
  const dir = await mkdtemp(join(tmpdir(), 'ci-report-missing-snapshot-'));
  t.after(() => rm(dir, { recursive: true, force: true }));
  const runs = [ { type: 'config', status: 'cached', url: suiteUrl + '/config' },
    { type: 'test', attribute: attr, status: 'cached', url: suiteUrl + '/test' } ];
  const manifest = { repository: 'example/repo', runs: [{ ...spec, observationFiles: ['missing.json'],
    offline: { suite: 'suite.json', checks: 'checks.json', jobs: 'jobs.json', actions: 'actions.json' } }] };
  for (const [file, data] of Object.entries({ 'manifest.json': manifest, 'suite.json': { commit: sha, status: 'success', runs },
    'jobs.json': [], 'actions.json': [], 'checks.json': runs.map((job, i) => ({ id: i + 1, details_url: job.url,
      started_at: '2026-09-06T00:00:00Z', completed_at: '2026-09-06T00:00:01Z', conclusion: 'success' })) }))
    await writeFile(join(dir, file), JSON.stringify(data));
  const report = await createReport(join(dir, 'manifest.json'), join(dir, 'out'));
  assert.equal(report.complete, false);
  assert.equal(report.suites[0].requiredPassed, true);
  assert.equal(report.suites[0].evidenceComplete, false);
  assert.ok(report.errors.some(error => error.stage === 'observation-evidence'));
  assert.equal(comparePairs([report.suites[0], { ...report.suites[0], variant: 'candidate', evidenceComplete: true }])[0].comparable, false);
});

test('a required runtime group cannot be satisfied by its cached runner package build', () => {
  const expectation = { ...spec, requiredAttributes: [attr] };
  const config = { type: 'config', status: 'success', checkStartedAt: '2026-09-06T00:00:00Z', checkCompletedAt: '2026-09-06T00:00:01Z' };
  const build = { type: 'build', attribute: attr, url: suiteUrl + '/build', status: 'cached',
    checkStartedAt: '2026-09-06T00:00:01Z', checkCompletedAt: '2026-09-06T00:00:02Z' };
  const omitted = summarizeSuite(expectation, { status: 'success' }, [{}], [config, build]);
  assert.deepEqual(omitted.missingRequiredAttributes, [attr]);
  assert.equal(omitted.requiredPassed, false);
  assert.equal(omitted.evidenceComplete, false);
  const present = summarizeSuite(expectation, { status: 'success' }, [{}], [config, build,
    { ...build, type: 'test', url: suiteUrl + '/test', status: 'success', checkCompletedAt: '2026-09-06T00:00:20Z' }]);
  assert.equal(present.requiredPassed, true);
  assert.equal(present.requiredSeconds, 20);
});

test('same-URL started-to-started log replacement survives success and portable offline replay', async t => {
  const dir = await mkdtemp(join(tmpdir(), 'ci-report-log-replacement-'));
  t.after(() => rm(dir, { recursive: true, force: true }));
  const runs = [{ type: 'config', status: 'cached', url: suiteUrl + '/config' },
    { type: 'test', status: 'success', attribute: attr, url: suiteUrl + '/test' }];
  const replacement = overlap.replaceAll('2026-09-06 00:00:', '2026-09-06 00:01:');
  const incident = { jobUrl: runs[1].url, observedAt: '2026-09-06T00:01:15Z', evidenceFile: 'incident.json' };
  const manifest = { repository: 'example/repo', runs: [{ ...spec, observationFiles: ['before.json', 'after.json'],
    offline: { suite: 'suite.json', checks: 'checks.json', jobs: 'jobs.json', actions: 'actions.json' } }] };
  const snapshot = { suites: [{ spec, suite: { commit: sha, status: 'started',
    runs: [runs[0], { ...runs[1], status: 'started' }] } }] };
  const evidence = { ...incident, status: 'started', file: 'prior.ndjson', workerSeconds: 999,
    replacement: { file: 'replacement.ndjson', observedStatus: 'started' } };
  for (const [file, value] of Object.entries({ 'manifest.json': manifest, 'incident.json': evidence,
    'suite.json': { commit: sha, status: 'success', runs }, 'actions.json': [],
    'jobs.json': [{ url: runs[1].url, file: 'current.ndjson' }],
    'checks.json': runs.map((job, i) => ({ id: i + 1, details_url: job.url, started_at: '2026-09-06T00:00:00Z',
      completed_at: i ? '2026-09-06T00:01:20Z' : '2026-09-06T00:00:01Z', conclusion: 'success' })),
    'before.json': { ...snapshot, at: '2026-09-06T00:00:15Z' }, 'after.json': { ...snapshot, at: incident.observedAt },
  })) await writeFile(join(dir, file), JSON.stringify(value));
  for (const [file, log] of Object.entries({ 'prior.ndjson': overlap, 'replacement.ndjson': replacement, 'current.ndjson': replacement }))
    await writeFile(join(dir, file), log);
  const finalOnly = await createReport(join(dir, 'manifest.json'), join(dir, 'final-only'));
  assert.equal(finalOnly.complete, true);
  manifest.runs[0].observedLogReplacements = [incident];
  await writeFile(join(dir, 'manifest.json'), JSON.stringify(manifest));
  const report = await createReport(join(dir, 'manifest.json'), join(dir, 'first'));
  const suite = report.suites[0];
  assert.equal(report.complete, false);
  assert.equal(suite.requiredPassed, true);
  assert.equal(suite.evidenceComplete, false);
  assert.equal(suite.observedResourceLowerBound, true);
  assert.equal(suite.interruptedJobs, 0);
  assert.equal(suite.jobs[1].checkAttempts.length, 1);
  assert.equal(suite.observedWorkerMinutes, finalOnly.suites[0].observedWorkerMinutes);
  assert.equal(suite.downloadReportedBytes, finalOnly.suites[0].downloadReportedBytes);
  assert.equal(suite.observedLogReplacements[0].evidence.workerSeconds, 14);
  assert.equal(suite.observedLogReplacements[0].evidence.replacement.workerSeconds, 14);
  assert.equal(report.incidents.filter(row => row.kind === 'log-replacement').length, 1);
  assert.equal(report.incidents.some(row => ['interrupted-worker', 'repeated-attempt'].includes(row.kind)), false);
  const replacementIncident = report.incidents.find(row => row.kind === 'log-replacement');
  assert.equal(replacementIncident.workerStartChanged, true);
  assert.equal(replacementIncident.preservesPriorPrefix, false);
  assert.equal(replacementIncident.cause, 'unknown');
  assert.equal(replacementIncident.observedAt, new Date(incident.observedAt).toISOString());
  const pair = comparePairs([finalOnly.suites[0], { ...suite, variant: 'candidate' }])[0];
  assert.equal(pair.comparable, false);
  assert.equal(pair.workerMinutesChangePercent, null);
  for (const file of ['incidents.json', 'incidents.csv', 'report.md'])
    assert.ok((await readFile(join(dir, 'first', file), 'utf8')).includes('log-replacement'));
  for (const file of ['prior.ndjson', 'replacement.ndjson', 'incident.json']) await rm(join(dir, file));
  incident.evidenceFile = 'first/raw/1/log-replacements/1.json';
  await writeFile(join(dir, 'manifest.json'), JSON.stringify(manifest));
  const replay = await createReport(join(dir, 'manifest.json'), join(dir, 'replay'));
  assert.equal(replay.errors.length, 0);
  assert.equal(replay.suites[0].evidenceComplete, false);
  assert.equal(replay.incidents.filter(row => row.kind === 'log-replacement').length, 1);
  assert.equal(await readFile(join(dir, 'replay/raw/1/log-replacements/1-prior.ndjson'), 'utf8'), overlap);
  assert.equal(await readFile(join(dir, 'replay/raw/1/log-replacements/1-replacement.ndjson'), 'utf8'), replacement);

  incident.evidenceFile = 'incident.json';
  await writeFile(join(dir, 'manifest.json'), JSON.stringify(manifest));
  await writeFile(join(dir, 'prior.ndjson'), overlap);
  const records = overlap.trim().split('\n');
  for (const [name, log, valid, metadata] of [
    ['identical', overlap, false, evidence],
    ['append-only', overlap + JSON.stringify({ utc_time: '2026-09-06T00:00:15Z', relative_nanoseconds: 15000000000, log_message: 'more output' }) + '\n', false, evidence],
    ['reformatted', records.map(line => JSON.stringify(Object.fromEntries(Object.entries(JSON.parse(line)).reverse()))).join('\n'), false, evidence],
    ['shortened-prefix', records[0] + '\n', true, evidence],
    ['wrong-url', replacement, false, { ...evidence, jobUrl: suiteUrl + '/other' }],
    ['wrong-time', replacement, false, { ...evidence, observedAt: '2026-09-06T00:01:16Z' }],
    ['invalid-log', '{"log_message":"missing clock"}', false, evidence],
  ]) {
    await writeFile(join(dir, 'incident.json'), JSON.stringify(metadata));
    await writeFile(join(dir, 'replacement.ndjson'), log);
    const result = await createReport(join(dir, 'manifest.json'), join(dir, name));
    assert.equal(result.incidents.some(row => row.kind === 'log-replacement'), valid, name);
    assert.equal(result.errors.some(row => row.stage === 'log-replacement-evidence'), !valid, name);
    assert.equal(result.suites[0].evidenceComplete, false, name);
    assert.equal(result.suites[0].observedResourceLowerBound, true, name);
    assert.equal(result.suites[0].interruptedJobs, 0, name);
    if (valid) assert.equal(result.incidents.find(row => row.kind === 'log-replacement').workerStartChanged, false);
  }
});

test('hopeless jobs are terminal failures with required end clocks, never accepted timing samples', () => {
  const config = { type: 'config', status: 'success', url: suiteUrl + '/config',
    checkStartedAt: '2026-09-06T00:00:00Z', checkCompletedAt: '2026-09-06T00:00:01Z' };
  const hopeless = { type: 'test', status: 'hopeless', attribute: attr, url: suiteUrl + '/test',
    checkStartedAt: '2026-09-06T00:00:01Z', checkCompletedAt: '2026-09-06T00:00:20Z', logStatus: 'ok' };
  const queued = { type: 'test', status: 'queued', attribute: attr + '-queued', url: suiteUrl + '/queued' };
  const stopped = summarizeSuite(spec, { status: 'failed' }, [{}], [config, hopeless, queued]);
  assert.equal(stopped.failedJobs, 1);
  assert.equal(stopped.interruptedJobs, 0);
  assert.equal(stopped.suiteSeconds, 20);
  assert.equal(stopped.requiredSeconds, null);
  assert.equal(stopped.requiredPassed, false);
  assert.equal(stopped.evidenceComplete, false);
  assert.equal(stopped.successfulTimingSample, false);
  assert.equal(stopped.incidents.find(row => row.kind === 'job-failure').observation, 'service reports this job as hopeless');
  assert.equal(stopped.incidents.find(row => row.kind === 'job-failure').cause, 'unknown');
  const terminal = summarizeSuite(spec, { status: 'hopeless' }, [{}], [config, hopeless]);
  assert.equal(terminal.suiteSeconds, 20);
  assert.equal(terminal.requiredSeconds, 20);
  assert.equal(terminal.requiredPassed, false);
  assert.equal(terminal.successfulTimingSample, false);
  const pair = comparePairs([terminal, { ...terminal, variant: 'candidate' }])[0];
  assert.equal(pair.comparable, false);
  assert.equal(pair.workerMinutesChangePercent, null);
  const noEnd = summarizeSuite(spec, { status: 'failed' }, [{}], [config, { ...hopeless, checkCompletedAt: null }]);
  assert.equal(noEnd.requiredSeconds, null);
  assert.equal(noEnd.evidenceComplete, false);
  assert.ok(noEnd.incidents.some(row => row.kind === 'missing-clock' && row.jobUrl === hopeless.url));
  for (const status of ['started', 'success']) {
    const recovered = summarizeSuite({ ...spec, observations: [
      { at: '2026-09-06T00:00:21Z', suite: { runs: [hopeless] } },
      { at: '2026-09-06T00:00:22Z', suite: { runs: [{ ...hopeless, status }] } },
    ] }, { status: 'success' }, [{}], [config, { ...hopeless, status: 'success' }]);
    assert.ok(recovered.incidents.some(row => row.kind === 'repeated-attempt' && row.previousStatus === 'hopeless' && row.status === status));
    assert.equal(recovered.evidenceComplete, false);
  }
});

test('failure timestamps use check completion or earliest matching status observation, never the last worker record', () => {
  for (const status of ['failed', 'failure', 'hopeless', 'abandoned']) {
    const job = { type: 'build', attribute: 'packages.x86_64-linux.crate-test-binaries-gammalooprs', url: suiteUrl + '/producer',
      status, checkStartedAt: '2026-09-08T18:00:00Z', workerCompletedAt: '2026-09-08T18:21:41Z' };
    const kind = status === 'abandoned' ? 'interrupted-worker' : 'job-failure';
    const observations = [
      { at: '2026-09-08T18:30:00Z', evidenceFile: 'later.json', suite: { runs: [job] } },
      { at: '2026-09-08T18:26:16.661Z', suite: { runs: [{ ...job, status: 'started' }] } },
      { at: '2026-09-08T18:20:00Z', suite: { runs: [{ ...job, url: suiteUrl + '/other' }] } },
      { at: '2026-09-08T18:27:17.025Z', evidenceFile: 'first.json', suite: { runs: [job] } },
    ];
    const unknown = summarizeSuite(spec, { status: 'failed' }, [], [job]).incidents.find(row => row.kind === kind);
    assert.equal(unknown.time, null);
    assert.equal(unknown.timeSource, null);
    const observed = summarizeSuite({ ...spec, observations }, { status: 'failed' }, [], [job]).incidents.find(row => row.kind === kind);
    assert.equal(observed.time, Date.parse('2026-09-08T18:27:17.025Z'));
    assert.equal(observed.timeSource, 'status-snapshot');
    assert.equal(observed.evidenceFile, 'first.json');
    assert.match(observed.observation, /not a termination time/);
    assert.equal(observed.cause, 'unknown');
    const completed = summarizeSuite({ ...spec, observations }, { status: 'failed' }, [], [
      { ...job, checkCompletedAt: '2026-09-08T18:26:30Z' },
    ]).incidents.find(row => row.kind === kind);
    assert.equal(completed.time, Date.parse('2026-09-08T18:26:30Z'));
    assert.equal(completed.timeSource, 'check-completion');
    assert.equal(completed.evidenceFile, undefined);
  }
});

test('Nextest timed-out tests contribute to failed outcomes and retain a separate timeout count', () => {
  const result = parseLog(JSON.stringify({ utc_time: '2026-09-08T20:00:00Z', relative_nanoseconds: 0,
    log_message: 'TIMEOUT [240.000s] gammaloop_integration_tests se1l_uv\n'
      + 'TIMEOUT [240.000s] gammaloop_integration_tests sunrise_scalar_1_uv\n'
      + 'Summary [240.125s] 104 tests run: 102 passed, 2 timed out, 67 skipped\n',
  }), { type: 'test', status: 'failed', attribute: attr });
  assert.deepEqual(result.testSummaries, [{ runner: 'nextest', executed: 104, passed: 102,
    failed: 2, timedOut: 2, skipped: 67, seconds: 240.125 }]);
  assert.deepEqual(result.tests.map(row => row.status), ['TIMEOUT', 'TIMEOUT']);
  assert.equal(result.testExecution, 'executed');
  const mixed = parseLog(JSON.stringify({ utc_time: '2026-09-08T20:00:00Z', relative_nanoseconds: 0,
    log_message: 'Summary [240.125s] 3 tests run: 0 passed, 1 failed, 2 timed out, 0 skipped\n'
      + 'test result: FAILED. 2 passed; 1 failed; 0 ignored; 0 measured; 0 filtered out; finished in 1.0s\n',
  }));
  assert.deepEqual(mixed.testSummaries.map(row => [row.failed, row.timedOut]), [[3, 2], [1, 0]]);
  assert.deepEqual(parseLog(overlap).testSummaries.map(row => [row.failed, row.timedOut]), [[0, 0], [0, 0]]);
});
