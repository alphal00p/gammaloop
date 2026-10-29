#!/usr/bin/env node

import { createHash } from "node:crypto";
import { readFileSync, readdirSync } from "node:fs";
import { basename, join } from "node:path";

const summaryPattern = /test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored;/g;
const runningPattern = /^running (\d+) tests?$/gm;
const testLinePattern = /^test (.+?) \.\.\.(.*)$/gm;
const standaloneOutcomePattern = /^(ok|FAILED)$/gm;

export function account(output, expectedLabel) {
  const summaries = [...output.matchAll(summaryPattern)].map((match) => ({
    status: match[1],
    passed: Number(match[2]),
    failed: Number(match[3]),
    ignored: Number(match[4]),
    selected: Number(match[2]) + Number(match[3]) + Number(match[4]),
    index: match.index,
    end: match.index + match[0].length,
  }));
  const running = [...output.matchAll(runningPattern)].map((match) => Number(match[1]));
  const testLines = [...output.matchAll(testLinePattern)].map((match) => ({
    label: match[1],
    outcome: ["ok", "FAILED"].includes(match[2].trim()) || match[2].trim().startsWith("ignored")
      ? match[2].trim()
      : null,
  }));
  const standaloneOutcomes = [...output.matchAll(standaloneOutcomePattern)].map((match) => match[1]);
  const reasons = [];
  if (summaries.length === 0) reasons.push("missing parent summary");
  const final = summaries.at(-1) ?? { passed: 0, failed: 0, ignored: 0, selected: 0 };
  if (final.selected !== 1) reasons.push(`outer selection is ${final.selected}, expected exactly one`);
  if (running[0] !== 1) reasons.push("missing exact one-test parent start");
  if (testLines.length === 0) reasons.push("missing test identity line");
  if (testLines.some((line) => line.label !== expectedLabel)) reasons.push("test identity mismatch");

  if (summaries.length === 1) {
    if (running.length !== 1 || testLines.length !== 1) reasons.push("ambiguous nested or prefix structure");
    const outcome = testLines[0]?.outcome ?? (standaloneOutcomes.length === 1 ? standaloneOutcomes[0] : null);
    if (final.passed === 1 && outcome !== "ok") reasons.push("pass summary lacks matching completed test line");
    if (final.failed === 1 && outcome !== "FAILED") reasons.push("failure summary lacks matching failed test line");
    if (final.ignored === 1 && !outcome?.startsWith("ignored")) reasons.push("ignored summary lacks matching ignored test line");
  } else if (summaries.length === 3) {
    const nestedShape = running.length === 3
      && running.every((count) => count === 1)
      && testLines.length === 3
      && testLines[0].outcome === null
      && testLines[1].outcome === "ok"
      && testLines[2].outcome === "ok"
      && summaries.every((summary) => summary.selected === 1 && summary.passed === 1 && summary.failed === 0 && summary.ignored === 0)
      && /^ok$/m.test(output.slice(summaries[1].end, summaries[2].index));
    if (!nestedShape) reasons.push("incomplete or ambiguous nested parent/child structure");
  } else {
    reasons.push(`unexpected summary count ${summaries.length}`);
  }
  const parseValid = reasons.length === 0;
  const classification = !parseValid ? "rejected"
    : final.passed === 1 ? "reported-pass"
      : final.failed === 1 ? "reported-failure"
        : final.ignored === 1 ? "reported-ignored"
          : "rejected";
  return {
    selected: final.selected,
    passed: final.passed,
    failed: final.failed,
    ignored: final.ignored,
    summaryCount: summaries.length,
    parseValid,
    classification,
    rejectionReasons: reasons,
    reportedCompletedPass: parseValid && classification === "reported-pass",
  };
}

function expect(output, expectedLabel, expected) {
  const actual = account(output, expectedLabel);
  for (const [key, value] of Object.entries(expected)) {
    if (JSON.stringify(actual[key]) !== JSON.stringify(value)) {
      throw new Error(`accounting mismatch: ${JSON.stringify({ key, actual, expected })}`);
    }
  }
}

function selfTest() {
  expect("running 1 test\ntest color::x ... ok\ntest result: ok. 1 passed; 0 failed; 0 ignored; 0 measured", "color::x", {
    classification: "reported-pass", reportedCompletedPass: true,
  });
  expect("running 1 test\ntest color::x ... FAILED\ntest result: FAILED. 0 passed; 1 failed; 0 ignored; 0 measured", "color::x", {
    classification: "reported-failure", reportedCompletedPass: false,
  });
  expect("running 2 tests\ntest color::x ... ok\ntest color::x_two ... ok\ntest result: ok. 2 passed; 0 failed; 0 ignored; 0 measured", "color::x", {
    classification: "rejected", reportedCompletedPass: false,
  });
  expect("running 1 test\ntest color::x ... ignored, reason\ntest result: ok. 0 passed; 0 failed; 1 ignored; 0 measured", "color::x", {
    classification: "reported-ignored", reportedCompletedPass: false,
  });
  expect("running 0 tests\ntest result: ok. 0 passed; 0 failed; 0 ignored; 0 measured", "color::x", {
    classification: "rejected", reportedCompletedPass: false,
  });
  expect("running 1 test\ntest color::other ... ok\ntest result: ok. 1 passed; 0 failed; 0 ignored; 0 measured", "color::x", {
    classification: "rejected", reportedCompletedPass: false,
  });
  expect("running 1 test\nrunning 1 test\ntest color::x ... ok\ntest result: ok. 1 passed; 0 failed; 0 ignored; 0 measured", "color::x", {
    classification: "rejected", reportedCompletedPass: false,
  });
  expect("running 1 test\ntest color::x ... ok", "color::x", {
    classification: "rejected", reportedCompletedPass: false,
  });
  expect("running 1 test\ntest color::x ... \nrunning 1 test\ntest color::x ... ok\ntest result: ok. 1 passed; 0 failed; 0 ignored; 0 measured\nrunning 1 test\ntest color::x ... ok\ntest result: ok. 1 passed; 0 failed; 0 ignored; 0 measured\nok\ntest result: ok. 1 passed; 0 failed; 0 ignored; 0 measured", "color::x", {
    classification: "reported-pass", reportedCompletedPass: true, summaryCount: 3,
  });
}

function sha256(path) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

function main() {
  if (process.argv[2] === "--self-test") {
    selfTest();
    console.log("accounting guards: PASS");
    return;
  }
  const [logRoot, elfPath] = process.argv.slice(2);
  if (!logRoot || !elfPath) {
    throw new Error("usage: account-color-runs.mjs LOG_ROOT ELF_PATH");
  }
  const elfSha256 = sha256(elfPath);
  const records = readdirSync(logRoot)
    .filter((name) => name.endsWith(".log"))
    .sort()
    .map((name) => {
      const label = basename(name, ".log").replaceAll("__", "::");
      const path = join(logRoot, name);
      return {
        label,
        expectedArgv: [elfPath, label, "--exact", "--test-threads=1", "--nocapture"],
        argvProvenance: "reconstructed expected argv; historical per-selection argv was not retained",
        exitStatus: null,
        exitStatusProvenance: "historical per-selection exit status was not retained",
        elfPath,
        elfSha256,
        logPath: path,
        logSha256: sha256(path),
        ...account(readFileSync(path, "utf8"), label),
      };
    });
  const totals = records.reduce((sum, row) => {
    for (const key of ["selected", "passed", "failed", "ignored"]) sum[key] += row[key];
    if (row.reportedCompletedPass) sum.reportedCompletedPasses += 1;
    return sum;
  }, { selected: 0, passed: 0, failed: 0, ignored: 0, reportedCompletedPasses: 0 });
  console.log(JSON.stringify({ schema: 1, records, totals }, null, 2));
}

main();
