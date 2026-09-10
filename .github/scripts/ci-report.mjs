#!/usr/bin/env node
// Collect existing CI evidence only; this command never dispatches or retries CI.
import { execFile } from 'node:child_process';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { dirname, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
import { isDeepStrictEqual, promisify } from 'node:util';

const exec = promisify(execFile);
const units = { B: 1, KiB: 1024, MiB: 1024 ** 2, GiB: 1024 ** 3, TiB: 1024 ** 4 };
const durationPattern = String.raw`(?:\d+(?:\.\d+)?(?:ms|s|m|h)\s*)+`;
const excluded = new Set(['cancelled', 'skipped', 'abandoned', 'queued', 'pending', 'in_progress', 'running', 'started']);
const sum = (rows, key) => rows.reduce((total, row) => total + (row[key] ?? 0), 0);
const timestamp = value => {
  if (typeof value !== 'string' || !value) return null;
  const normalized = value.replace(' ', 'T').replace(/(\.\d{3})\d+/, '$1');
  const parsed = Date.parse(/(?:Z|[+-]\d\d:\d\d)$/.test(normalized) ? normalized : `${normalized}Z`);
  return Number.isFinite(parsed) ? parsed : null;
};
const secondsBetween = (start, end) => start != null && end != null && end >= start ? (end - start) / 1000 : null;
const duration = value => [...value.matchAll(/([\d.]+)(ms|s|m|h)/g)]
  .reduce((total, match) => total + Number(match[1]) * ({ ms: 0.001, s: 1, m: 60, h: 3600 }[match[2]]), 0);
const min = values => values.length ? Math.min(...values) : null;
const max = values => values.length ? Math.max(...values) : null;

// Each completion event's timer gives an interval ending at its log timestamp.
// Union within one worker; concurrent workers still consume separate resources.
export function intervalSeconds(events) {
  const intervals = events.map(event => [event.time - event.seconds * 1000, event.time])
    .sort((a, b) => a[0] - b[0]);
  let end = -Infinity;
  let total = 0;
  for (const [start, finish] of intervals) {
    total += Math.max(0, finish - Math.max(start, end));
    end = Math.max(end, finish);
  }
  return total / 1000;
}

export function parseLog(ndjson, job = {}) {
  const records = ndjson.split('\n').filter(line => line.trim()).map((line, index) => {
    let record;
    try { record = JSON.parse(line); } catch { throw new Error(`invalid NDJSON at line ${index + 1}`); }
    if (typeof record.log_message !== 'string' || timestamp(record.utc_time) == null)
      throw new Error(`invalid log record at line ${index + 1}`, { cause: timestamp(record.utc_time) == null ? 'missing-clock' : 'invalid-record' });
    return record;
  });
  if (!records.length) throw new Error('empty log');
  const downloads = [], uploads = [], transferTimeouts = [], transferErrors = [], incidents = [], builds = [], cargo = [], compiled = [], checked = [], tests = [], testSummaries = [];
  for (const [index, record] of records.entries()) {
    const text = record.log_message.replace(/\x1b\[[0-?]*[ -/]*[@-~]/g, '');
    const time = timestamp(record.utc_time);
    const previous = records[index - 1];
    if (previous && (time < timestamp(previous.utc_time)
      || (Number.isFinite(record.relative_nanoseconds) && Number.isFinite(previous.relative_nanoseconds)
        && record.relative_nanoseconds < previous.relative_nanoseconds)))
      incidents.push({ kind: 'invalid-clock', time, line: index + 1, cause: 'unknown', observation: 'worker clock moves backwards' });
    // Retry warnings are observations, not completed transfers or measured bytes.
    // Keep numeric diagnostics only; request URLs can contain credentials.
    for (const match of text.matchAll(/unable to (upload|download) '[^'\n]+'[^\n]*?Operation too slow\. Less than ([\d.]+) bytes\/sec transferred the last (\d+) seconds(?:; retrying in (\d+) ms \(attempt (\d+)\/(\d+)\))?/g))
      transferTimeouts.push({ direction: match[1], time, minimumBytesPerSecond: Number(match[2]), windowSeconds: Number(match[3]),
        retryDelayMilliseconds: match[4] == null ? null : Number(match[4]),
        attempt: match[5] == null ? null : Number(match[5]), maximumAttempts: match[6] == null ? null : Number(match[6]) });
    for (const match of text.matchAll(/unable to (upload|download) '[^'\n]+'([^\n]*)/g)) {
      if (/Operation too slow/.test(match[2])) continue;
      const kind = /HTTP\/2|HTTP2|framing layer/i.test(match[2]) ? 'http2-transfer-error' : 'cache-transfer-error';
      transferErrors.push({ kind, direction: match[1], time, httpStatus: Number(match[2].match(/HTTP error (\d{3})/)?.[1]) || null,
        retryDelayMilliseconds: Number(match[2].match(/retrying in (\d+) ms/)?.[1]) || null,
        attempt: Number(match[2].match(/attempt (\d+)\//)?.[1]) || null,
        maximumAttempts: Number(match[2].match(/attempt \d+\/(\d+)/)?.[1]) || null });
    }
    // Ignore "Downloading cached": only the completed line has size and duration.
    for (const match of text.matchAll(new RegExp(`Downloaded cached ([^\\n]+?) \\(([\\d.]+) (B|KiB|MiB|GiB|TiB)\\) in (${durationPattern})`, 'g')))
      downloads.push({ name: match[1], reportedBytes: Number(match[2]) * units[match[3]], seconds: duration(match[4]), time });
    for (const match of text.matchAll(new RegExp(`Uploaded ([^\\n]+?) in (${durationPattern})`, 'g'))) {
      const size = match[1].match(/^(.*?) \(([\d.]+) (B|KiB|MiB|GiB|TiB)\)$/);
      uploads.push({ name: size?.[1] ?? match[1], reportedBytes: size ? Number(size[2]) * units[size[3]] : null, seconds: duration(match[2]), time });
    }
    for (const match of text.matchAll(/(?:Building |building ['"]?)(\/nix\/store\/[0-9a-z]{32}-[^\s'"]+\.drv)/g))
      builds.push({ drv: match[1], time });
    for (const match of text.matchAll(new RegExp(`Finished [^\\n]+? in (${durationPattern})`, 'g')))
      cargo.push({ seconds: duration(match[1]), time });
    for (const match of text.matchAll(/\bCompiling ([^\n]+)/g)) compiled.push(match[1].trim());
    for (const match of text.matchAll(/\bChecking ([\w-]+ v[^\s]+(?: \([^\n]+\))?)/g)) checked.push(match[1]);
    for (const match of text.matchAll(/Summary\s*\[\s*([\d.]+)s\]\s*(\d+) tests? run:([^\n]*)/g)) {
      const count = label => Number(match[3].match(new RegExp(`(\\d+) ${label}`))?.[1] ?? 0);
      const timedOut = count('timed out');
      testSummaries.push({ runner: 'nextest', executed: Number(match[2]), passed: count('passed'),
        failed: count('failed') + timedOut, timedOut, skipped: count('skipped'), seconds: Number(match[1]) });
    }
    for (const match of text.matchAll(/test result: (?:ok|FAILED)\. (\d+) passed; (\d+) failed; (\d+) ignored;[^\n]*?finished in ([\d.]+)s/g))
      testSummaries.push({ runner: 'libtest', executed: Number(match[1]) + Number(match[2]), passed: Number(match[1]), failed: Number(match[2]), timedOut: 0, skipped: Number(match[3]), seconds: Number(match[4]) });
    for (const match of text.matchAll(/\b(PASS|FAIL|SKIP|TIMEOUT)\s+\[[^\]\n]*\]\s+([^\n]+)/g))
      tests.push({ status: match[1], name: match[2].trim() });
    for (const match of text.matchAll(/(?:^|\n|> )test ([^\n]+?) \.\.\. (ok|FAILED|ignored)\b/g))
      tests.push({ status: ({ ok: 'PASS', FAILED: 'FAIL', ignored: 'SKIP' })[match[2]], name: match[1] });
  }
  const first = records[0], last = records.at(-1);
  const relativeSeconds = Number.isFinite(first.relative_nanoseconds) && Number.isFinite(last.relative_nanoseconds)
    ? (last.relative_nanoseconds - first.relative_nanoseconds) / 1e9 : null;
  const workerSeconds = incidents.length ? null : relativeSeconds ?? secondsBetween(timestamp(first.utc_time), timestamp(last.utc_time));
  if (workerSeconds != null && workerSeconds < 0) throw new Error('invalid worker time range', { cause: 'invalid-clock' });
  for (const event of transferTimeouts) incidents.push({ kind: 'transfer-timeout', ...event, cause: 'unknown', observation: 'transfer failed its reported low-speed threshold' });
  for (const event of transferErrors) incidents.push({ ...event, cause: 'unknown', observation: event.kind === 'http2-transfer-error'
    ? 'transfer reported an HTTP/2 stream or framing error' : 'cache transfer reported a failure' });
  const resultName = job.attribute?.split('.').at(-1).replace(/^nix-ci-check-/, '');
  const testExecution = job.type !== 'test' ? 'not-applicable'
    : testSummaries.length || tests.some(test => test.status !== 'SKIP') ? 'executed'
    : job.status === 'cached' || downloads.some(event => event.name === resultName) ? 'reused' : 'unknown';
  return {
    records: records.length, workerStartedAt: first.utc_time, workerCompletedAt: last.utc_time, workerSeconds,
    downloadCount: downloads.length, downloadReportedBytes: sum(downloads, 'reportedBytes'),
    downloadActiveSeconds: intervalSeconds(downloads), uploadCount: uploads.length,
    uploadReportedBytes: uploads.every(event => event.reportedBytes != null) ? sum(uploads, 'reportedBytes') : null,
    uploadActiveSeconds: intervalSeconds(uploads), transferActiveSeconds: intervalSeconds([...downloads, ...uploads]),
    cargoFinishedSeconds: sum(cargo, 'seconds'), cargoActiveSeconds: intervalSeconds(cargo),
    compilationMessages: compiled.length, compiled: [...new Set(compiled)],
    checkingMessages: checked.length, checked: [...new Set(checked)],
    transferTimeoutCount: transferTimeouts.length, transferTimeouts, transferErrorCount: transferErrors.length, transferErrors, incidents,
    testExecution, testSummaries, tests, downloads, uploads, builds, cargo,
  };
}

export function summarizeSuite(spec, suite, checks, jobs) {
  const completed = job => secondsBetween(timestamp(job.checkStartedAt), timestamp(job.checkCompletedAt)) != null
    && !job.checkResultMismatch && !excluded.has(job.status) && !excluded.has(job.checkConclusion);
  const actual = jobs.filter(completed);
  const starts = jobs.filter(job => job.type === 'config').map(job => timestamp(job.checkStartedAt)).filter(value => value != null);
  const started = min(starts);
  const finished = max(actual.map(job => timestamp(job.checkCompletedAt)));
  const wanted = jobs.filter(job => spec.requiredAttributes
    ? spec.requiredAttributes.includes(job.attribute)
      && (!/^packages\.[^.]+\.nix-ci-check-/.test(job.attribute ?? '') || job.type === 'test')
    : job.type === 'test' || /^checks\.[^.]+\.gammaloop-(clippy|fmt|guppy-workspace-graph)$/.test(job.attribute ?? ''));
  // A failed earlier attempt remains in resource totals; latest attempt owns the result.
  const latestWanted = [...new Map(wanted.toSorted((a, b) => (a.attempt ?? 1) - (b.attempt ?? 1)
    || (timestamp(a.checkStartedAt) ?? 0) - (timestamp(b.checkStartedAt) ?? 0))
    .map(job => [job.attribute, job])).values()];
  const missingRequiredAttributes = (spec.requiredAttributes ?? []).filter(attribute => !latestWanted.some(job => job.attribute === attribute));
  const wantedComplete = latestWanted.length > 0 && latestWanted.every(completed) && missingRequiredAttributes.length === 0;
  const wantedEnd = wantedComplete ? max(latestWanted.map(job => timestamp(job.checkCompletedAt))) : null;
  const aggregateEnd = max(actual.filter(job => job.type === 'deploy').map(job => timestamp(job.checkCompletedAt)));
  const builtBy = new Map();
  for (const job of jobs) for (const event of job.builds ?? []) {
    if (!builtBy.has(event.drv)) builtBy.set(event.drv, new Set());
    builtBy.get(event.drv).add(job.url);
  }
  const repeatedDerivations = [...builtBy].filter(([, workers]) => workers.size > 1)
    .map(([drv, workers]) => ({ drv, workerCount: workers.size, jobUrls: [...workers] }));
  const snapshots = spec.observations ?? [];
  const incidents = jobs.flatMap(job => (job.incidents ?? []).map(incident => ({ jobUrl: job.url, logUrl: job.url ? job.url + '/logs' : null, attribute: job.attribute, ...incident })));
  for (const job of jobs) {
    const start = timestamp(job.checkStartedAt), end = timestamp(job.checkCompletedAt);
    if (['failed', 'failure', 'abandoned', 'hopeless'].includes(job.status)) {
      const observed = snapshots.filter(snapshot => timestamp(snapshot.at) != null
        && job.url && snapshot.suite.runs.some(run => run.url === job.url && run.status === job.status))
        .toSorted((a, b) => timestamp(a.at) - timestamp(b.at))[0];
      incidents.push({ kind: job.status === 'abandoned' ? 'interrupted-worker' : 'job-failure',
        time: end ?? timestamp(observed?.at), timeSource: end != null ? 'check-completion' : observed ? 'status-snapshot' : null,
        evidenceFile: end == null ? observed?.evidenceFile : undefined, jobUrl: job.url, logUrl: job.url ? job.url + '/logs' : null,
        attribute: job.attribute, cause: 'unknown', observation: 'service reports this job as ' + job.status
          + (end == null && observed ? '; timestamp is the first saved status observation, not a termination time' : '') });
    }
    const needsStart = !['queued', 'pending', 'skipped', 'cancelled'].includes(job.status);
    const needsEnd = ['success', 'cached', 'failed', 'failure', 'hopeless'].includes(job.status);
    if ((needsStart && start == null) || (needsEnd && end == null) || (start != null && end != null && end < start))
      incidents.push({ kind: start != null && end != null ? 'invalid-clock' : 'missing-clock', time: end ?? start,
        jobUrl: job.url, attribute: job.attribute, cause: 'unknown', observation: 'check timestamps are missing, invalid, or reversed' });
    if ((job.checkAttempts?.length ?? 0) > 1)
      incidents.push({ kind: 'repeated-attempt', time: start, jobUrl: job.url, logUrl: job.url ? job.url + '/logs' : null, attribute: job.attribute,
        checkIds: job.checkAttempts.map(check => check.id), cause: 'unknown', observation: 'multiple check attempts share one job URL; earlier logs may have been replaced' });
  }
  const attempts = new Map();
  for (const job of jobs) {
    const key = JSON.stringify([job.type, job.attribute ?? null]);
    if (!attempts.has(key)) attempts.set(key, new Map());
    if (job.url) attempts.get(key).set(job.url, job);
  }
  const prior = new Map(), queued = new Map();
  for (const snapshot of snapshots) {
    for (const job of snapshot.suite.runs) {
      const previous = prior.get(job.url);
      if (previous && ((['success', 'cached', 'failure', 'failed', 'abandoned', 'cancelled', 'hopeless'].includes(previous.status)
        && ['queued', 'pending', 'running', 'started', 'in_progress'].includes(job.status))
        || (previous.status === 'abandoned' && job.status !== 'abandoned')
        || (['failed', 'failure', 'cancelled', 'hopeless'].includes(previous.status) && ['success', 'cached'].includes(job.status))))
        incidents.push({ kind: 'repeated-attempt', time: timestamp(snapshot.at), jobUrl: job.url, logUrl: job.url ? job.url + '/logs' : null,
          attribute: job.attribute ?? job.type, previousStatus: previous.status, status: job.status, evidenceFile: snapshot.evidenceFile,
          cause: 'unknown', observation: 'saved snapshots show a job URL restarting or recovering; trigger and replaced work are unknown' });
      prior.set(job.url, job);
      const key = JSON.stringify([job.type, job.attribute ?? null]);
      if (!attempts.has(key)) attempts.set(key, new Map());
      if (!attempts.get(key).has(job.url)) attempts.get(key).set(job.url, job);
      const prerequisites = spec.dependencies?.[job.attribute];
      if (job.status !== 'queued' || !prerequisites?.length || !prerequisites.every(attribute => {
        const matches = snapshot.suite.runs.filter(run => run.type === 'build' && run.attribute === attribute);
        return matches.length === 1 && ['success', 'cached'].includes(matches[0].status);
      })) continue;
      if (!queued.has(job.url)) queued.set(job.url, { kind: 'queued-after-declared-prerequisites', time: timestamp(snapshot.at),
        jobUrl: job.url, attribute: job.attribute, prerequisiteAttributes: prerequisites, observations: [], evidenceFile: snapshot.evidenceFile,
        dependencyEvidenceFile: spec.dependencyEvidenceFile ?? null, cause: 'unknown',
        observation: 'queued in a saved snapshot while all declared prerequisites were successful or cached; hidden dependencies and scheduling cause are unknown' });
      queued.get(job.url).observations.push({ at: snapshot.at, evidenceFile: snapshot.evidenceFile });
    }
  }
  incidents.push(...queued.values());
  for (const rows of attempts.values()) if (rows.size > 1) {
    const attempts = [...rows.values()];
    incidents.push({ kind: 'repeated-attempt', time: min(attempts.map(job => timestamp(job.checkStartedAt)).filter(time => time != null)),
      attribute: attempts[0].attribute ?? attempts[0].type, type: attempts[0].type, jobUrl: attempts.at(-1).url, jobUrls: [...rows.keys()], attemptCount: rows.size, cause: 'unknown',
      observation: 'multiple job URLs have the same type and attribute in this suite; dispatch or retry trigger is unknown' });
  }
  for (const event of repeatedDerivations) {
    const starts = jobs.flatMap(job => (job.builds ?? []).filter(build => build.drv === event.drv)
      .map(build => ({ jobUrl: job.url, time: build.time })));
    incidents.push({ kind: 'repeated-derivation', time: min(starts.map(start => start.time)), ...event, buildStarts: starts,
      jobUrl: event.jobUrls[0], logUrl: event.jobUrls[0] + '/logs', cause: 'unknown',
      observation: 'the same full derivation was built under multiple job URLs; substitution or scheduling cause is unknown' });
  }
  const finalAggregateTailSeconds = secondsBetween(wantedEnd, aggregateEnd);
  if (suite.status === 'success' && finalAggregateTailSeconds > 60)
    incidents.push({ kind: 'final-success-tail', time: aggregateEnd, seconds: finalAggregateTailSeconds,
      jobUrl: actual.find(job => job.type === 'deploy' && timestamp(job.checkCompletedAt) === aggregateEnd)?.url,
      cause: 'unknown', observation: 'final success completed more than 60 seconds after the last required check' });
  const submitted = timestamp(spec.submittedAt), submissionCompleted = timestamp(spec.submissionCompletedAt);
  if (submitted != null && ((wantedEnd != null && wantedEnd < submitted)
    || (submissionCompleted != null && submissionCompleted < submitted)))
    incidents.push({ kind: 'invalid-clock', time: submitted, cause: 'unknown', observation: 'submission timestamps disagree with completion order' });
  const observed = jobs.filter(job => job.logStatus === 'ok');
  const missing = jobs.filter(job => job.logStatus === 'missing');
  const observedInterruptions = spec.observedInterruptions ?? [];
  for (const incident of observedInterruptions) incidents.push({ kind: 'interrupted-worker', time: timestamp(incident.observedAt),
    jobUrl: incident.jobUrl, logUrl: incident.jobUrl + '/logs', evidenceFile: incident.evidenceFile,
    cause: 'unknown', observation: 'saved evidence records an abandoned worker; current logs may replace earlier work' });
  const observedLogReplacements = spec.observedLogReplacements ?? [];
  for (const incident of observedLogReplacements.filter(item => item.verified)) incidents.push({ kind: 'log-replacement',
    time: timestamp(incident.observedAt), jobUrl: incident.jobUrl, logUrl: incident.jobUrl + '/logs',
    evidenceFile: incident.evidenceFile, workerStartChanged: incident.evidence.verification.workerStartChanged,
    preservesPriorPrefix: incident.evidence.verification.preservesPriorPrefix, cause: 'unknown',
    observation: 'saved raw logs for the same job URL have a different worker start or replaced record prefix; termination and retry trigger are unknown' });
  const interrupted = new Set([...jobs.filter(job => job.status === 'abandoned').map(job => job.url),
    ...observedInterruptions.map(incident => incident.jobUrl),
    ...snapshots.flatMap(snapshot => snapshot.suite.runs.filter(job => job.status === 'abandoned').map(job => job.url))]);
  const uncertainHistory = observedLogReplacements.length > 0 || jobs.some(job => (job.checkAttempts?.length ?? 0) > 1)
    || [...prior.keys()].some(url => !jobs.some(job => job.url === url))
    || incidents.some(incident => incident.kind === 'repeated-attempt' && incident.previousStatus);
  const invalidChecks = incidents.some(incident => ['invalid-clock', 'missing-clock', 'check-result-mismatch'].includes(incident.kind));
  const completedSuite = ['success', 'failure', 'failed', 'hopeless'].includes(suite.status);
  return {
    layout: spec.layout, variant: spec.variant, scenario: spec.scenario, pair: spec.pair ?? '1', sha: spec.sha, suiteUrl: spec.suiteUrl,
    status: suite.status, successfulTimingSample: suite.status === 'success',
    suiteStartedAt: started == null ? null : new Date(started).toISOString(),
    suiteCompletedAt: completedSuite && finished != null ? new Date(finished).toISOString() : null,
    suiteSeconds: completedSuite ? secondsBetween(started, finished) : null,
    observedCheckSpanSeconds: secondsBetween(started, finished), requiredSeconds: secondsBetween(started, wantedEnd),
    submittedAt: spec.submittedAt ?? null, submissionCompletedAt: spec.submissionCompletedAt ?? null,
    submissionToRequiredSeconds: secondsBetween(submitted, wantedEnd),
    submissionDurationSeconds: secondsBetween(submitted, submissionCompleted),
    submissionToRequiredBoundsSeconds: submitted != null && submissionCompleted != null
      ? { minimum: wantedEnd == null ? null : Math.max(0, (wantedEnd - submissionCompleted) / 1000),
        maximum: secondsBetween(submitted, wantedEnd) } : null,
    requiredPassed: wantedComplete && latestWanted.every(job => ['success', 'cached'].includes(job.status)),
    requiredAttributes: spec.requiredAttributes ?? latestWanted.map(job => job.attribute),
    observedRequiredAttributes: latestWanted.map(job => job.attribute), missingRequiredAttributes,
    finalAggregateSeconds: secondsBetween(started, aggregateEnd), finalAggregateTailSeconds,
    scheduledJobs: jobs.length, cachedJobs: jobs.filter(job => job.status === 'cached').length,
    failedJobs: jobs.filter(job => ['failed', 'failure', 'hopeless'].includes(job.status)).length,
    cancelledOrSkippedJobs: jobs.filter(job => ['cancelled', 'skipped'].includes(job.status)).length,
    observedWorkers: observed.length, missingLogs: missing.length, interruptedJobs: interrupted.size, observedInterruptions, observedLogReplacements,
    observedResourceLowerBound: !completedSuite || missing.length > 0 || interrupted.size > 0 || uncertainHistory || invalidChecks,
    evidenceComplete: completedSuite && wantedComplete && missing.length === 0 && interrupted.size === 0 && checks.length > 0 && !uncertainHistory && !invalidChecks,
    observedWorkerMinutes: sum(observed, 'workerSeconds') / 60,
    downloadReportedBytes: sum(observed, 'downloadReportedBytes'),
    intermediateDownloadReportedBytes: sum(observed.filter(job => job.context === 'artifact-producer'), 'downloadReportedBytes'),
    uploadReportedBytes: observed.every(job => job.uploadReportedBytes != null) ? sum(observed, 'uploadReportedBytes') : null,
    downloadWorkerSeconds: sum(observed, 'downloadActiveSeconds'), uploadWorkerSeconds: sum(observed, 'uploadActiveSeconds'),
    transferWorkerSeconds: sum(observed, 'transferActiveSeconds'), transferTimeoutCount: sum(observed, 'transferTimeoutCount'),
    transferErrorCount: sum(observed, 'transferErrorCount'), incidentCount: incidents.length, incidents,
    cargoFinishedSeconds: sum(observed, 'cargoFinishedSeconds'),
    compilationMessages: sum(observed, 'compilationMessages'), checkingMessages: sum(observed, 'checkingMessages'),
    executedTestJobs: jobs.filter(job => job.testExecution === 'executed').length,
    reusedTestJobs: jobs.filter(job => job.testExecution === 'reused').length,
    unknownTestJobs: jobs.filter(job => job.testExecution === 'unknown').length,
    repeatedDerivationCount: repeatedDerivations.length,
    repeatedBuildWorkers: repeatedDerivations.reduce((total, event) => total + event.workerCount - 1, 0),
    repeatedDerivations,
  };
}

export function summarizeActions(run) {
  const jobs = (run.jobs ?? []).map(job => ({
    id: job.id, name: job.name, status: job.status, conclusion: job.conclusion,
    runnerName: job.runner_name, runnerLabels: job.labels,
    startedAt: job.started_at, completedAt: job.completed_at,
    workerSeconds: job.conclusion === 'skipped' || job.status !== 'completed' ? null
      : secondsBetween(timestamp(job.started_at), timestamp(job.completed_at)),
    steps: (job.steps ?? []).map(step => ({
      name: step.name, conclusion: step.conclusion,
      seconds: excluded.has(step.conclusion) || step.status !== 'completed' ? null
        : secondsBetween(timestamp(step.started_at), timestamp(step.completed_at)),
    })),
  }));
  const done = run.status === 'completed' && !excluded.has(run.conclusion);
  const lastJob = max(jobs.filter(job => job.workerSeconds != null && !excluded.has(job.conclusion)).map(job => timestamp(job.completedAt)));
  return {
    runId: run.id, attempt: run.run_attempt, name: run.name, sha: run.head_sha,
    url: `${run.html_url}/attempts/${run.run_attempt}`, event: run.event,
    status: run.status, conclusion: run.conclusion, successfulTimingSample: done && run.conclusion === 'success',
    createdAt: run.created_at, attemptStartedAt: run.run_started_at, updatedAt: run.updated_at,
    // created_at belongs to the original workflow run, including for reruns.
    createdToUpdatedSeconds: done ? secondsBetween(timestamp(run.created_at), timestamp(run.updated_at)) : null,
    attemptToLastJobSeconds: done ? secondsBetween(timestamp(run.run_started_at), lastJob) : null,
    observedWorkerMinutes: sum(jobs, 'workerSeconds') / 60, jobs,
  };
}

export function comparePairs(suites) {
  const pairs = new Map();
  for (const suite of suites) {
    const key = JSON.stringify([suite.layout, suite.scenario, suite.pair]);
    if (!pairs.has(key)) pairs.set(key, { layout: suite.layout, scenario: suite.scenario, pair: suite.pair });
    pairs.get(key)[suite.variant] = suite;
  }
  return [...pairs.values()].map(({ layout, scenario, pair, baseline, candidate }) => {
    const comparable = Boolean(baseline?.evidenceComplete && candidate?.evidenceComplete
      && baseline.successfulTimingSample && candidate.successfulTimingSample && baseline.requiredPassed && candidate.requiredPassed);
    const change = key => comparable && baseline[key] > 0 && candidate[key] != null
      ? (candidate[key] / baseline[key] - 1) * 100 : null;
    const latest = suite => new Map((suite?.jobs ?? []).filter(job => job.type === 'test')
      .toSorted((a, b) => (a.attempt ?? 1) - (b.attempt ?? 1)
        || (timestamp(a.checkStartedAt) ?? 0) - (timestamp(b.checkStartedAt) ?? 0)).map(job => [job.attribute, job]));
    const beforeGroups = latest(baseline), afterGroups = latest(candidate);
    const groups = [...new Set([...beforeGroups.keys(), ...afterGroups.keys()])].map(attribute => {
      const before = beforeGroups.get(attribute), after = afterGroups.get(attribute);
      const baselineCompletedSeconds = before?.checkSeconds != null
        ? secondsBetween(timestamp(baseline.suiteStartedAt), timestamp(before.checkCompletedAt)) : null;
      const candidateCompletedSeconds = after?.checkSeconds != null
        ? secondsBetween(timestamp(candidate.suiteStartedAt), timestamp(after.checkCompletedAt)) : null;
      return { attribute, baselineCheckSeconds: before?.checkSeconds ?? null, candidateCheckSeconds: after?.checkSeconds ?? null,
        baselineCompletedSeconds, candidateCompletedSeconds,
        completionChangeSeconds: comparable && baselineCompletedSeconds != null && candidateCompletedSeconds != null
          ? candidateCompletedSeconds - baselineCompletedSeconds : null };
    });
    return {
      layout, scenario, pair, baselineUrl: baseline?.suiteUrl ?? null, candidateUrl: candidate?.suiteUrl ?? null,
      comparable, coverageEquivalent: null,
      suiteChangePercent: change('suiteSeconds'), requiredChangePercent: change('requiredSeconds'),
      submissionToRequiredChangePercent: change('submissionToRequiredSeconds'),
      workerMinutesChangePercent: change('observedWorkerMinutes'), restoreChangePercent: change('downloadReportedBytes'),
      intermediateRestoreChangePercent: change('intermediateDownloadReportedBytes'),
      uploadTimeChangePercent: change('uploadWorkerSeconds'), groups,
    };
  });
}

class Collector {
  constructor(manifestPath, outputDir) {
    this.manifestDir = dirname(resolve(manifestPath));
    this.outputDir = resolve(outputDir);
    this.errors = [];
  }

  async json(path, base = this.manifestDir) {
    return JSON.parse(await readFile(resolve(base, path), 'utf8'));
  }

  async save(path, value) {
    const target = resolve(this.outputDir, path);
    await mkdir(dirname(target), { recursive: true });
    await writeFile(target, typeof value === 'string' ? value : `${JSON.stringify(value, null, 2)}\n`, { mode: 0o600 });
    return target;
  }

  async github(endpoint, paginate = false) {
    const args = ['api', endpoint, '-H', 'Accept: application/vnd.github+json'];
    if (paginate) args.push('--paginate', '--slurp');
    const { stdout } = await exec('gh', args, { maxBuffer: 64 * 1024 * 1024, timeout: 120_000 });
    return JSON.parse(stdout);
  }

  async nixci(url, accept) {
    const response = await fetch(url, { headers: { Accept: accept }, signal: AbortSignal.timeout(90_000) });
    if (!response.ok) throw new Error(`HTTP ${response.status} from ${url}`);
    return response.text();
  }

  async collect(spec, index, repository) {
    const prefix = `raw/${index + 1}`;
    const local = spec.offline;
    let suite, checks = [], actionRuns = [], localJobs = [];
    try {
      suite = local ? await this.json(local.suite) : JSON.parse(await this.nixci(spec.suiteUrl, 'application/json'));
      if (suite.commit !== spec.sha) throw new Error(`suite SHA ${suite.commit} does not match manifest SHA ${spec.sha}`);
      if (!Array.isArray(suite.runs)) throw new Error('suite response has no runs array');
      await this.save(`${prefix}/suite.json`, suite);
      if (!local) await this.save(`${prefix}/snapshot.json`, { at: new Date().toISOString(), suites: [{ spec: { sha: spec.sha, suiteUrl: spec.suiteUrl }, suite }] });
    } catch (error) {
      this.errors.push({ suiteUrl: spec.suiteUrl, stage: 'suite', error: error.message });
      return null;
    }
    try {
      const pages = local ? await this.json(local.checks)
        : await this.github(`repos/${repository}/commits/${spec.sha}/check-runs?filter=all&per_page=100`, true);
      checks = (Array.isArray(pages) ? pages : [pages]).flatMap(page => page.check_runs ?? [page]);
      if (!checks.length) throw new Error('no GitHub check timestamps found');
      await this.save(`${prefix}/checks.json`, checks);
    } catch (error) { this.errors.push({ suiteUrl: spec.suiteUrl, stage: 'checks', error: error.message }); }
    if (local?.jobs) localJobs = await this.json(local.jobs);
    try {
      if (local) {
        actionRuns = await this.json(local.actions);
      } else {
        const pages = await this.github(`repos/${repository}/actions/runs?head_sha=${spec.sha}&per_page=100`, true);
        const runs = pages.flatMap(page => page.workflow_runs);
        if (spec.actionRunIds?.some(id => !runs.some(run => run.id === id))) throw new Error('selected Actions run ID not found for this SHA');
        for (const run of runs) {
          if (spec.actionRunIds && !spec.actionRunIds.includes(run.id)) continue;
          for (let attempt = 1; attempt <= run.run_attempt; attempt++) {
            const endpoint = `repos/${repository}/actions/runs/${run.id}/attempts/${attempt}`;
            const metadata = await this.github(endpoint);
            const jobPages = await this.github(`${endpoint}/jobs?per_page=100`, true);
            actionRuns.push({ ...metadata, jobs: jobPages.flatMap(page => page.jobs) });
          }
        }
      }
      if (!Array.isArray(actionRuns) || actionRuns.some(run => !Array.isArray(run.jobs) || !Number.isInteger(run.run_attempt)))
        throw new Error('actions evidence must be an array of run attempts with jobs and run_attempt');
      await this.save(`${prefix}/actions.json`, actionRuns);
    } catch (error) { this.errors.push({ suiteUrl: spec.suiteUrl, stage: 'actions', error: error.message }); }

    const observations = [];
    let dependencies, dependencyEvidenceFile;
    try {
      if (spec.dependencySnapshot) {
        const evidence = await this.json(spec.dependencySnapshot.file);
        const mapping = evidence.dependencies;
        if (!mapping || typeof mapping !== 'object' || Array.isArray(mapping)
          || Object.values(mapping).some(rows => !Array.isArray(rows) || rows.some(attribute => typeof attribute !== 'string')))
          throw new Error('dependency snapshot requires the generated dependencies mapping');
        dependencies = mapping;
        dependencyEvidenceFile = await this.save(prefix + '/dependency-snapshot.json', evidence);
      }
      for (const [index, file] of (spec.observationFiles ?? []).entries()) {
        const snapshot = await this.json(file);
        const entries = snapshot.suites?.filter(entry => entry.spec?.suiteUrl === spec.suiteUrl);
        const entry = entries?.[0];
        if (timestamp(snapshot.at) == null || entries?.length !== 1 || entry.spec.sha !== spec.sha
          || entry.suite?.commit !== spec.sha || !Array.isArray(entry.suite.runs)
          || entry.suite.runs.some(job => typeof job.url !== 'string' || !job.url.startsWith(spec.suiteUrl + '/')
            || !/^[^/?#]+$/.test(job.url.slice(spec.suiteUrl.length + 1))))
          throw new Error('observation requires a timestamp and exactly one matching suite with exact-suite job URLs');
        const evidenceFile = await this.save(prefix + '/observations/' + (index + 1) + '.json', { at: snapshot.at, suites: [entry] });
        observations.push({ at: snapshot.at, suite: entry.suite, evidenceFile });
      }
      observations.sort((a, b) => timestamp(a.at) - timestamp(b.at));
    } catch (error) { this.errors.push({ suiteUrl: spec.suiteUrl, stage: 'observation-evidence', error: error.message }); }

    const observedInterruptions = [];
    for (const [position, incident] of (spec.observedInterruptions ?? []).entries()) {
      const observation = { ...incident };
      observedInterruptions.push(observation);
      try {
        const evidence = await this.json(incident.evidenceFile);
        if (evidence?.jobUrl !== incident.jobUrl || evidence.status !== 'abandoned')
          throw new Error('incident evidence must identify the same abandoned job URL');
        if (evidence.file) {
          const log = await readFile(resolve(dirname(resolve(this.manifestDir, incident.evidenceFile)), evidence.file), 'utf8');
          await this.save(`${prefix}/interruptions/${position + 1}.ndjson`, log);
          evidence.file = `${position + 1}.ndjson`;
        }
        observation.evidenceFile = await this.save(`${prefix}/interruptions/${position + 1}.json`, evidence);
        observation.evidence = evidence;
      } catch (error) {
        this.errors.push({ suiteUrl: spec.suiteUrl, jobUrl: incident.jobUrl, stage: 'interruption-evidence', error: error.message });
      }
      this.errors.push({ suiteUrl: spec.suiteUrl, jobUrl: incident.jobUrl, stage: 'observed-interruption',
        error: 'previously observed abandonment: current logs may replace earlier work; prior observations are kept separate and totals remain lower bounds' });
    }

    const observedLogReplacements = [];
    for (const [position, incident] of (spec.observedLogReplacements ?? []).entries()) {
      const observation = { ...incident, verified: false };
      observedLogReplacements.push(observation);
      try {
        const evidence = await this.json(incident.evidenceFile);
        if (evidence?.jobUrl !== incident.jobUrl || timestamp(evidence.observedAt) !== timestamp(incident.observedAt)
          || typeof evidence.file !== 'string' || !evidence.file
          || typeof evidence.replacement?.file !== 'string' || !evidence.replacement.file)
          throw new Error('log replacement evidence requires the same job URL and observation time, plus both raw log files');
        const directory = dirname(resolve(this.manifestDir, incident.evidenceFile));
        const logs = await Promise.all([evidence.file, evidence.replacement.file].map(file => readFile(resolve(directory, file), 'utf8')));
        const paths = [position + 1 + '-prior.ndjson', position + 1 + '-replacement.ndjson'];
        for (const [index, log] of logs.entries()) await this.save(prefix + '/log-replacements/' + paths[index], log);
        const metrics = logs.map(log => parseLog(log, evidence));
        const records = logs.map(log => log.split('\n').filter(line => line.trim()).map(line => JSON.parse(line)));
        const workerStartChanged = timestamp(metrics[0].workerStartedAt) !== timestamp(metrics[1].workerStartedAt);
        const preservesPriorPrefix = records[0].length <= records[1].length
          && records[0].every((record, index) => {
            const current = records[1][index];
            return isDeepStrictEqual(record, current) || index === records[0].length - 1
              && current.log_message.startsWith(record.log_message)
              && isDeepStrictEqual({ ...record, log_message: current.log_message }, current);
          });
        if (preservesPriorPrefix) throw new Error('raw logs are identical or append-only; they do not prove log replacement');
        observation.evidence = { ...evidence, ...metrics[0], file: paths[0],
          replacement: { ...evidence.replacement, ...metrics[1], file: paths[1] },
          verification: { workerStartChanged, preservesPriorPrefix } };
        observation.evidenceFile = await this.save(prefix + '/log-replacements/' + (position + 1) + '.json', observation.evidence);
        observation.verified = true;
      } catch (error) {
        this.errors.push({ suiteUrl: spec.suiteUrl, jobUrl: incident.jobUrl, stage: 'log-replacement-evidence', error: error.message });
      }
    }

    const runs = [...new Map(suite.runs.map(job => [job.url, job])).values()];
    const retryMetadata = new Map(), retryEvidence = new Map();
    const retryRuns = runs.filter(run => runs.some(other => other.type === run.type && other.attribute === run.attribute
      && (other.attempt ?? 1) > 1));
    let metadataCursor = 0;
    await Promise.all(Array.from({ length: 4 }, async () => {
      while (metadataCursor < retryRuns.length) {
        const run = retryRuns[metadataCursor++];
        try {
          if (!run.url?.startsWith(spec.suiteUrl + '/') || !/^[^/?#]+$/.test(run.url.slice(spec.suiteUrl.length + 1)))
            throw new Error('retry metadata URL does not belong to the exact suite');
          const metadata = local ? localJobs.find(job => job.url === run.url)?.retryMetadata
            : JSON.parse(await this.nixci(run.url, 'application/json'));
          if (!metadata && local) continue;
          if (metadata?.uuid !== run.url.split('/').at(-1) || metadata.attribute !== run.attribute
            || metadata.type !== run.type || metadata.status !== run.status)
            throw new Error('retry metadata identity or status disagrees with the suite');
          retryEvidence.set(run.url, await this.save(`${prefix}/retry-metadata/${runs.indexOf(run) + 1}.json`, metadata));
          retryMetadata.set(run.url, metadata);
        } catch (error) {
          this.errors.push({ suiteUrl: spec.suiteUrl, jobUrl: run.url, stage: 'retry-metadata', error: error.message });
        }
      }
    }));
    const checkJobs = new Map(checks.map(check => [check, check.details_url?.replace(/\/$/, '')]));
    const checkMappings = new Map();
    for (const check of checks) {
      const source = runs.find(run => run.url === checkJobs.get(check));
      if (!source) continue;
      const chain = [source];
      let current = source;
      while (retryMetadata.get(current.url)?.retried_by) {
        const next = runs.find(run => run.url === retryMetadata.get(current.url).retried_by);
        if (!next || chain.includes(next) || retryMetadata.get(next.url)?.retry_of !== current.url
          || next.attribute !== source.attribute || next.type !== source.type
          || (next.attempt ?? 1) <= (current.attempt ?? 1)) break;
        chain.push(next);
        current = next;
      }
      // NixCI can update a check for a retry while retaining the original details URL.
      // Only reciprocal service links can transfer that clock; prior clocks stay unknown.
      if (current === source || retryMetadata.get(current.url)?.retried_by || source.status === current.status
        || check.status !== 'completed' || check.name !== `${current.type} ${current.attribute}`
        || checks.filter(other => other.details_url?.replace(/\/$/, '') === source.url).length !== 1
        || !(check.conclusion === 'success' && ['success', 'cached'].includes(current.status)
          || check.conclusion === 'failure' && ['failed', 'failure', 'hopeless'].includes(current.status))
        || checks.some(other => other !== check && other.details_url?.replace(/\/$/, '') === current.url)) continue;
      checkJobs.set(check, current.url);
      checkMappings.set(check, { source: 'reciprocal-retry-links', fromUrl: source.url,
        urls: chain.map(run => run.url), evidenceFiles: chain.map(run => retryEvidence.get(run.url)) });
    }
    const jobs = new Array(runs.length);
    let cursor = 0;
    // Bound requests and memory use while retaining manifest/API order in outputs.
    await Promise.all(Array.from({ length: 4 }, async () => {
      while (cursor < runs.length) {
        const position = cursor++;
        const run = runs[position];
        const matching = checks.filter(check => checkJobs.get(check) === run.url?.replace(/\/$/, ''));
        const check = matching.toSorted((a, b) => (timestamp(a.started_at) ?? 0) - (timestamp(b.started_at) ?? 0)).at(-1);
        const superseded = [...checkMappings.values()].find(mapping => mapping.fromUrl === run.url)?.urls.at(-1);
        const job = {
          attribute: run.attribute ?? run.type, type: run.type, status: run.status, url: run.url, attempt: run.attempt ?? 1,
          context: /crate-(?:test-|deps-)|cargoArtifacts|Artifacts|prebuild|nextest-binaries|ci-test-inputs/.test(run.attribute ?? '') ? 'artifact-producer'
            : /doctest/.test(run.attribute ?? '') ? 'doctest' : /clippy/.test(run.attribute ?? '') ? 'clippy' : run.type,
          checkAttempts: [...new Map(matching.filter(check => check.id != null).map(check => [check.id, { id: check.id, startedAt: check.started_at, completedAt: check.completed_at }])).values()],
          checkId: check?.id ?? null, checkStartedAt: check?.started_at ?? null,
          checkCompletedAt: check?.completed_at ?? null, checkConclusion: check?.conclusion ?? null,
          checkDetailsUrl: check?.details_url ?? null, checkMapping: checkMappings.get(check) ?? null,
          checkSupersededBy: superseded ?? null,
          retryMetadata: retryMetadata.get(run.url), retryMetadataFile: retryEvidence.get(run.url),
          checkSeconds: excluded.has(run.status) || excluded.has(check?.conclusion) ? null
            : secondsBetween(timestamp(check?.started_at), timestamp(check?.completed_at)),
          logStatus: 'not-run', testExecution: run.type === 'test' ? 'unknown' : 'not-applicable',
        };
        if (!check && !['queued', 'pending', 'skipped', 'cancelled'].includes(run.status))
          this.errors.push({ suiteUrl: spec.suiteUrl, jobUrl: run.url, stage: 'checks', error: superseded
            ? 'GitHub check now describes an explicitly linked retry; original attempt timestamps are unavailable'
            : 'no check matches the exact job URL' });
        if (run.status === 'cached') {
          job.logStatus = 'cached';
          if (run.type === 'test') job.testExecution = 'reused';
        } else if (!['queued', 'pending', 'skipped'].includes(run.status) && !(run.status === 'cancelled' && !check?.started_at)) {
          try {
            if (!run.url?.startsWith(`${spec.suiteUrl}/`)) throw new Error('job URL does not belong to the exact suite URL');
            const stored = localJobs.find(row => row.url === run.url);
            if (local && !stored?.file) throw new Error('no offline log file for this job URL');
            const ndjson = local
              ? await readFile(resolve(dirname(resolve(this.manifestDir, local.jobs)), stored.file), 'utf8')
              : await this.nixci(`${run.url}/logs`, 'application/x-ndjson');
            job.file = await this.save(`${prefix}/logs/${position + 1}.ndjson`, ndjson);
            Object.assign(job, parseLog(ndjson, job), { logStatus: 'ok' });
            job.preWorkerSeconds = secondsBetween(timestamp(job.checkStartedAt), timestamp(job.workerStartedAt));
            job.postWorkerSeconds = job.checkSeconds == null ? null
              : secondsBetween(timestamp(job.workerCompletedAt), timestamp(job.checkCompletedAt));
          } catch (error) {
            job.logStatus = 'missing';
            job.incidents = [{ kind: ['missing-clock', 'invalid-clock'].includes(error.cause) ? error.cause : 'missing-log',
              time: timestamp(job.checkStartedAt), cause: 'unknown', observation: 'worker evidence is missing or malformed; timings are unknown' }];
            this.errors.push({ suiteUrl: spec.suiteUrl, jobUrl: run.url, stage: 'log', error: error.message });
          }
        }
        job.checkResultMismatch = ['failed', 'failure', 'hopeless'].includes(run.status) && check?.conclusion === 'success'
          || ['success', 'cached'].includes(run.status) && check?.conclusion === 'failure';
        if (job.checkResultMismatch) {
          // A retry can update a GitHub check while its details URL still names the failed attempt.
          // Retain both raw outcomes, but do not assign that check's duration to this worker.
          job.checkSeconds = job.preWorkerSeconds = job.postWorkerSeconds = null;
          job.incidents ??= [];
          job.incidents.push({ kind: 'check-result-mismatch', time: timestamp(check.completed_at), timeSource: 'check-metadata',
            status: run.status, checkConclusion: check.conclusion, cause: 'unknown',
            observation: 'NixCI job outcome disagrees with the GitHub check attached to its exact URL; attempt clocks are untrusted' });
        }
        if (run.status === 'abandoned') this.errors.push({ suiteUrl: spec.suiteUrl, jobUrl: run.url, stage: 'interrupted-worker',
          error: 'abandoned worker: observed resource totals are lower bounds, including after a successful retry' });
        jobs[position] = job;
      }
    }));
    await this.save(`${prefix}/jobs.json`, jobs);
    const matchingChecks = checks.filter(check => runs.some(run => check.details_url?.replace(/\/$/, '') === run.url));
    const summary = summarizeSuite({ ...spec, observedInterruptions, observedLogReplacements, observations, dependencies, dependencyEvidenceFile }, suite, matchingChecks, jobs);
    for (const job of jobs) job.completedFromSuiteSeconds = job.checkSeconds == null ? null
      : secondsBetween(timestamp(summary.suiteStartedAt), timestamp(job.checkCompletedAt));
    if (summary.missingRequiredAttributes.length) this.errors.push({ suiteUrl: spec.suiteUrl, stage: 'required-checks',
      error: 'missing required attributes: ' + summary.missingRequiredAttributes.join(', ') });
    summary.evidenceComplete &&= !this.errors.some(error => error.suiteUrl === spec.suiteUrl);
    return { ...summary, jobs, actions: actionRuns.filter(run => run.head_sha === spec.sha).map(summarizeActions) };
  }
}

const csv = rows => {
  const columns = [...new Set(rows.flatMap(Object.keys))];
  const cell = value => `"${String(value == null ? '' : typeof value === 'object' ? JSON.stringify(value) : value).replaceAll('"', '""')}"`;
  return `${columns.map(cell).join(',')}\n${rows.map(row => columns.map(column => cell(row[column])).join(',')).join('\n')}\n`;
};
const display = value => value == null ? 'unknown' : typeof value === 'number' ? value.toFixed(2) : String(value).replaceAll('|', '\\|').replaceAll('\n', ' ');

export async function createReport(manifestPath, outputDir) {
  const collector = new Collector(manifestPath, outputDir);
  const manifest = JSON.parse(await readFile(manifestPath, 'utf8'));
  if (!/^[\w.-]+\/[\w.-]+$/.test(manifest.repository) || !Array.isArray(manifest.runs) || !manifest.runs.length)
    throw new Error('manifest requires repository (owner/repo) and a nonempty runs array');
  const identities = new Set();
  const pairs = new Set();
  for (const spec of manifest.runs) {
    for (const key of ['layout', 'variant', 'scenario', 'sha', 'suiteUrl'])
      if (typeof spec[key] !== 'string' || !spec[key]) throw new Error(`run requires ${key}`);
    if (!/^[0-9a-f]{40}$/.test(spec.sha)) throw new Error('sha must be a full 40-character commit SHA');
    if (!['baseline', 'candidate'].includes(spec.variant)) throw new Error('variant must be baseline or candidate');
    const pair = JSON.stringify([spec.layout, spec.scenario, spec.pair ?? '1', spec.variant]);
    if (pairs.has(pair)) throw new Error('duplicate layout/scenario/pair/variant; give repeated trials different pair values');
    pairs.add(pair);
    if (spec.requiredAttributes && (!Array.isArray(spec.requiredAttributes) || !spec.requiredAttributes.length
      || !spec.requiredAttributes.every(attribute => typeof attribute === 'string' && attribute.length)))
      throw new Error('requiredAttributes must be a nonempty array of attribute names');
    for (const key of ['observedInterruptions', 'observedLogReplacements'])
      if (spec[key] && (!Array.isArray(spec[key])
        || spec[key].some(incident => !incident || typeof incident.jobUrl !== 'string'
          || !incident.jobUrl.startsWith(spec.suiteUrl + '/') || !/^[^/?#]+$/.test(incident.jobUrl.slice(spec.suiteUrl.length + 1))
          || typeof incident.observedAt !== 'string' || timestamp(incident.observedAt) == null || typeof incident.evidenceFile !== 'string' || !incident.evidenceFile)))
        throw new Error(key + ' requires exact-suite jobUrl, observedAt, and evidenceFile');
    for (const key of ['submittedAt', 'submissionCompletedAt']) if (spec[key] != null && timestamp(spec[key]) == null)
      throw new Error(key + ' must be a recorded UTC timestamp');
    if (spec.submissionCompletedAt && !spec.submittedAt) throw new Error('submissionCompletedAt requires submittedAt');
    if (spec.observationFiles && (!Array.isArray(spec.observationFiles) || spec.observationFiles.some(file => typeof file !== 'string' || !file)))
      throw new Error('observationFiles must contain saved snapshot paths');
    if (spec.dependencySnapshot && (spec.dependencySnapshot.sha !== spec.sha || typeof spec.dependencySnapshot.file !== 'string' || !spec.dependencySnapshot.file))
      throw new Error('dependencySnapshot requires the same SHA and a generated config file');
    if (spec.actionRunIds && (!Array.isArray(spec.actionRunIds) || !spec.actionRunIds.every(Number.isSafeInteger)))
      throw new Error('actionRunIds must contain integer run IDs');
    if (spec.offline) for (const key of ['suite', 'checks', 'jobs', 'actions'])
      if (typeof spec.offline[key] !== 'string') throw new Error('offline evidence requires ' + key + ' file');
    const url = new URL(spec.suiteUrl);
    if (url.origin !== 'https://nix-ci.com' || url.search || url.hash
      || !url.pathname.startsWith('/gh:' + manifest.repository.replace('/', ':') + '/') || !url.pathname.endsWith(`/${spec.sha}`))
      throw new Error('suiteUrl must be an exact https://nix-ci.com suite URL ending in the full SHA');
    if (identities.has(spec.suiteUrl)) throw new Error('duplicate suite URL in manifest');
    identities.add(spec.suiteUrl);
  }
  await collector.save('manifest.json', manifest);
  const suites = [];
  for (const [index, spec] of manifest.runs.entries()) {
    const suite = await collector.collect(spec, index, manifest.repository);
    if (suite) suites.push(suite);
  }
  const report = {
    schemaVersion: 1, repository: manifest.repository, collectedAt: new Date().toISOString(),
    complete: collector.errors.length === 0 && suites.every(suite => suite.evidenceComplete),
    limits: [
      'Worker minutes are the sum of observed worker log spans, not billed compute or CHF.',
      'Missing logs and abandoned workers make observed resource totals lower bounds, even after a successful retry; cached and skipped jobs are not missing workers.',
      'Explicit interruption and verified log replacement observations preserve overwritten history separately; their earlier resources are never added to current-log totals because streams may overlap.',
      'Download sizes are log-reported artifact sizes, not established compressed network bytes. Upload bytes may be unknown.',
      'Transfer timers include an unspecified combination of network, decompression and store import; pre-worker gaps are not proven queue time.',
      'Intervals are unioned per worker before summing across workers; download and upload intervals may overlap.',
      'Cargo Finished times are command wall times, not compiler CPU time. Repeated crate names do not establish identical Cargo units.',
      'Repeated exact .drv build starts across job URLs establish repeated build work, not physical worker identity.',
      'Incidents identify observed symptoms; transport errors and repeated attempts do not establish their underlying cause or retry trigger.',
      'Queued-after-prerequisites incidents require a saved suite snapshot and SHA-matched declared dependencies; snapshots do not establish continuous queue duration or hidden prerequisites.',
      'Submission clocks use recorded push-call bounds, never commit author timestamps; submission-to-required time includes the submission call when only its start is known.',
      'Executed/reused/unknown describes visible test evidence. A cached wrapper alone does not establish result reuse; test summaries are not a complete coverage inventory.',
      'GitHub created_at persists across attempts. Compare attempt-to-last-job time for reruns; cancelled/skipped attempts have no completed timing.',
    ],
    suites, pairs: comparePairs(suites), errors: collector.errors,
    incidents: suites.flatMap(suite => suite.incidents.map(incident => ({ layout: suite.layout, variant: suite.variant, scenario: suite.scenario,
      pair: suite.pair, sha: suite.sha, suiteUrl: suite.suiteUrl, ...incident,
      observedAt: incident.time == null ? null : new Date(incident.time).toISOString() })))
      .toSorted((a, b) => (a.time ?? Infinity) - (b.time ?? Infinity)),
  };
  await collector.save('report.json', report);
  await collector.save('incidents.json', report.incidents);
  const suiteRows = [], jobRows = [], actionRows = [], transferRows = [], repeatedRows = [], compilationRows = [];
  for (const suite of suites) {
    const { jobs, actions, repeatedDerivations, incidents, ...summary } = suite;
    suiteRows.push(summary);
    const identity = { layout: suite.layout, variant: suite.variant, scenario: suite.scenario, pair: suite.pair, sha: suite.sha, suiteUrl: suite.suiteUrl };
    for (const job of jobs) {
      const { downloads, uploads, builds, cargo, tests, testSummaries, compiled, checked, incidents, ...metrics } = job;
      jobRows.push({ ...identity, ...metrics });
      for (const [phase, crates] of [['compile', compiled], ['check', checked]])
        for (const crate of crates ?? []) compilationRows.push({ ...identity, jobUrl: job.url, attribute: job.attribute, context: job.context, phase, crate });
      for (const [direction, events] of [['download', downloads], ['upload', uploads]])
        for (const event of events ?? []) transferRows.push({ ...identity, jobUrl: job.url, attribute: job.attribute, direction, ...event });
    }
    for (const { jobs, ...action } of actions) actionRows.push({ ...identity, ...action });
    for (const repeated of repeatedDerivations) repeatedRows.push({ ...identity, ...repeated });
  }
  for (const [name, rows] of Object.entries({ suites: suiteRows, jobs: jobRows, actions: actionRows, transfers: transferRows, compilations: compilationRows,
    pairs: report.pairs.map(({ groups, ...pair }) => pair), 'repeated-derivations': repeatedRows, incidents: report.incidents }))
    await collector.save(`${name}.csv`, csv(rows));
  const lines = [
    '# CI measurements', '',
    `Evidence: **${report.complete ? 'complete' : 'INCOMPLETE'}**. See report.json for every attempt, job, and error.`, '',
    '| Layout / variant / scenario | SHA | Result | Suite s | Required s | Final tail s | Observed worker min | Restore GiB* | Upload worker s | Repeated drv | Tests executed / reused / unknown |',
    '|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|',
    ...suites.map(suite => `| ${[suite.layout, suite.variant, suite.scenario].map(display).join(' / ')} | [${suite.sha.slice(0, 9)}](${suite.suiteUrl}) | ${display(suite.status)} | ${display(suite.suiteSeconds)} | ${display(suite.requiredSeconds)} | ${display(suite.finalAggregateTailSeconds)} | ${display(suite.observedWorkerMinutes)} | ${display(suite.downloadReportedBytes / 1024 ** 3)} | ${display(suite.uploadWorkerSeconds)} | ${suite.repeatedDerivationCount} | ${suite.executedTestJobs} / ${suite.reusedTestJobs} / ${suite.unknownTestJobs} |`),
    '', '*Log-reported artifact sizes; not established network bytes.', '',
    ...suites.flatMap(suite => suite.observedInterruptions.map(incident =>
      `- Prior interruption: [${display(incident.jobUrl.split('/').at(-1))}](${incident.jobUrl}) observed ${display(incident.observedAt)}; ${display(incident.evidence?.workerSeconds)} worker seconds and ${display(incident.evidence?.downloadReportedBytes)} reported download bytes in [saved evidence](${incident.evidenceFile}). Kept separate from current totals; pair is incomplete.`)),
    ...suites.filter(suite => suite.transferTimeoutCount > 0).map(suite =>
      `- [${display(suite.layout)} / ${display(suite.variant)} / ${display(suite.scenario)}](${suite.suiteUrl}): ${suite.transferTimeoutCount} transfer timeout warnings. Unfinished attempts are separate from completed-transfer totals; see per-job numeric diagnostics in report.json.`),
    '',
    '| Incident | Observed UTC | Job / suite | Evidence |',
    '|---|---|---|---|',
    ...report.incidents.map(incident => '| ' + [incident.kind, incident.observedAt].map(display).join(' | ')
      + ' | [' + display(incident.attribute ?? incident.jobUrl?.split('/').at(-1) ?? incident.layout) + '](' + (incident.jobUrl ?? incident.suiteUrl) + ') | '
      + display(incident.observation) + '; cause: ' + display(incident.cause)
      + (incident.logUrl ? ' ([log](' + incident.logUrl + '))' : '')
      + (incident.evidenceFile ? ' ([saved evidence](' + incident.evidenceFile + '))' : '') + ' |'),
    '',
    '| Layout / scenario / pair | Comparable | Worker change % | Intermediate restore change % | Required latency change % |',
    '|---|---|---:|---:|---:|',
    ...report.pairs.map(pair => '| ' + [pair.layout, pair.scenario, pair.pair].map(display).join(' / ') + ' | ' + [pair.comparable, pair.workerMinutesChangePercent, pair.intermediateRestoreChangePercent, pair.requiredChangePercent].map(display).join(' | ') + ' |'),
    '', 'Negative changes indicate improvement. Comparable pairs require successful suites and complete evidence; coverage equivalence requires a separate inventory check.', '',
    '| Layout / variant / scenario | Submitted UTC | Submission → required s | Submission call s |',
    '|---|---|---:|---:|',
    ...suites.map(suite => '| ' + [suite.layout, suite.variant, suite.scenario].map(display).join(' / ') + ' | '
      + [suite.submittedAt, suite.submissionToRequiredSeconds, suite.submissionDurationSeconds].map(display).join(' | ') + ' |'),
    '',
    '| Layout / variant / scenario | Workflow attempt | Result | Created → updated s | Attempt → last job s | Observed worker min |',
    '|---|---|---|---:|---:|---:|',
    ...actionRows.map(run => `| ${[run.layout, run.variant, run.scenario].map(display).join(' / ')} | [${display(run.name)} #${run.runId}/${run.attempt}](${run.url}) | ${display(run.conclusion ?? run.status)} | ${display(run.createdToUpdatedSeconds)} | ${display(run.attemptToLastJobSeconds)} | ${display(run.observedWorkerMinutes)} |`),
    '', ...report.limits.map(limit => `- ${limit}`), '',
    ...report.errors.map(error => `- Missing evidence (${error.stage}): ${display(error.jobUrl ?? error.suiteUrl)} — ${display(error.error)}`), '',
  ];
  await collector.save('report.md', lines.join('\n'));
  return report;
}

if (process.argv[1] && import.meta.url === pathToFileURL(resolve(process.argv[1])).href) {
  if (process.argv.length !== 4) {
    console.error('Usage: node .github/scripts/ci-report.mjs MANIFEST OUTPUT_DIR');
    process.exitCode = 1;
  } else {
    try {
      const report = await createReport(process.argv[2], process.argv[3]);
      console.log(`${report.suites.length} suites; ${report.errors.length} evidence errors; ${resolve(process.argv[3], 'report.md')}`);
      if (!report.complete) process.exitCode = 2;
    } catch (error) {
      console.error(error.message);
      process.exitCode = 1;
    }
  }
}
