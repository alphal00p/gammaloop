#!/usr/bin/env node

import { createHash } from "node:crypto";
import { readFileSync, readdirSync } from "node:fs";
import { basename, join } from "node:path";

const summaryPattern = /test result: \w+\. (\d+) passed; (\d+) failed; (\d+) ignored;/g;

export function account(output) {
  const matches = [...output.matchAll(summaryPattern)];
  if (matches.length === 0) {
    throw new Error("expected at least one libtest summary");
  }
  const final = matches.at(-1);
  const passed = Number(final[1]);
  const failed = Number(final[2]);
  const ignored = Number(final[3]);
  const selected = passed + failed + ignored;
  return {
    selected,
    passed,
    failed,
    ignored,
    summaryCount: matches.length,
    creditedPass: selected > 0 && passed > 0 && failed === 0 && ignored === 0,
  };
}

function expect(output, expected) {
  const actual = account(output);
  if (JSON.stringify(actual) !== JSON.stringify(expected)) {
    throw new Error(`accounting mismatch: ${JSON.stringify({ actual, expected })}`);
  }
}

function selfTest() {
  expect("test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured", {
    selected: 1, passed: 1, failed: 0, ignored: 0, summaryCount: 1, creditedPass: true,
  });
  expect("test result: FAILED. 0 passed; 1 failed; 0 ignored; 0 measured", {
    selected: 1, passed: 0, failed: 1, ignored: 0, summaryCount: 1, creditedPass: false,
  });
  expect("test result: FAILED. 1 passed; 1 failed; 0 ignored; 0 measured", {
    selected: 2, passed: 1, failed: 1, ignored: 0, summaryCount: 1, creditedPass: false,
  });
  expect("test result: ok. 0 passed; 0 failed; 1 ignored; 0 measured", {
    selected: 1, passed: 0, failed: 0, ignored: 1, summaryCount: 1, creditedPass: false,
  });
  expect("test result: ok. 0 passed; 0 failed; 0 ignored; 0 measured", {
    selected: 0, passed: 0, failed: 0, ignored: 0, summaryCount: 1, creditedPass: false,
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
        argv: [elfPath, label, "--exact", "--test-threads=1", "--nocapture"],
        elfPath,
        elfSha256,
        logPath: path,
        logSha256: sha256(path),
        ...account(readFileSync(path, "utf8")),
      };
    });
  const totals = records.reduce((sum, row) => {
    for (const key of ["selected", "passed", "failed", "ignored"]) sum[key] += row[key];
    if (row.creditedPass) sum.creditedPasses += 1;
    return sum;
  }, { selected: 0, passed: 0, failed: 0, ignored: 0, creditedPasses: 0 });
  console.log(JSON.stringify({ schema: 1, records, totals }, null, 2));
}

main();
