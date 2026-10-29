#!/usr/bin/env node

import { createHash } from "node:crypto";
import { readFileSync } from "node:fs";

const expectedLogSha256 = "bac31cd092deef94f1d7ccf873f7629a6796cf7bb4216c09dc0bd9645e001a5d";
const expectedTargetPasses = [67, 1, 10, 5, 45, 2, 1, 3, 18, 15, 3, 32];
const summaryPattern = /^test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out;/;

function sha256(bytes) {
  return createHash("sha256").update(bytes).digest("hex");
}

function parseSummary(line, lineNumber) {
  const match = line.match(summaryPattern);
  if (!match) return null;
  return {
    line: lineNumber,
    status: match[1],
    passed: Number(match[2]),
    failed: Number(match[3]),
    ignored: Number(match[4]),
    measured: Number(match[5]),
    filtered: Number(match[6]),
  };
}

function classifyHeader(line, lineNumber) {
  const cargo = line.match(/^\s+Running (.+?) \((\/[^)]+)\)$/);
  if (cargo) return { kind: "cargo-target", line: lineNumber, name: cargo[1], executable: cargo[2], summaries: [] };
  const doctest = line.match(/^\s+Doc-tests (\S+)$/);
  if (doctest) return { kind: "doctest", line: lineNumber, name: doctest[1], summaries: [] };
  return null;
}

const [logPath] = process.argv.slice(2);
if (!logPath || process.argv.length !== 3) throw new Error("usage: account-cargo-targets.mjs LOG");
const bytes = readFileSync(logPath);
if (sha256(bytes) !== expectedLogSha256) throw new Error("affected-consumers log hash mismatch");

const sections = [];
let section = null;
for (const [index, line] of bytes.toString("utf8").split(/\r?\n/).entries()) {
  const header = classifyHeader(line, index + 1);
  if (header) {
    section = header;
    sections.push(section);
    continue;
  }
  const summary = parseSummary(line, index + 1);
  if (summary && section) section.summaries.push(summary);
}

const targets = sections.filter(({ kind }) => kind === "cargo-target").map((target) => {
  if (target.summaries.length === 0) throw new Error(`target ${target.name} has no outer summary`);
  const outer = target.summaries.at(-1);
  return { ...target, outer, child_summaries: target.summaries.slice(0, -1) };
});
const doctests = sections.filter(({ kind }) => kind === "doctest").map((target) => {
  if (target.summaries.length !== 1) throw new Error(`doctest ${target.name} has ambiguous summary count`);
  return { name: target.name, line: target.line, summary: target.summaries[0] };
});

if (targets.length !== 12) throw new Error(`expected 12 Cargo targets, found ${targets.length}`);
if (JSON.stringify(targets.map(({ outer }) => outer.passed)) !== JSON.stringify(expectedTargetPasses)) {
  throw new Error("Cargo target pass totals differ from the frozen target-boundary tally");
}
const bootstrap = targets.find(({ name }) => name === "tests/vector_registration_bootstrap.rs");
if (!bootstrap || bootstrap.child_summaries.length !== 1 || bootstrap.child_summaries[0].passed !== 1) {
  throw new Error("expected one separately reported bootstrap child execution");
}

const total = (field) => targets.reduce((sum, target) => sum + target.outer[field], 0);
const result = {
  schema: 1,
  log: { path: logPath, sha256: expectedLogSha256 },
  accounting_boundary: "Each Cargo Running header contributes exactly its final outer test-result summary. Earlier summaries within that boundary are retained as child executions without additional top-level pass credit.",
  cargo_targets: targets,
  totals: { targets: targets.length, passed: total("passed"), failed: total("failed"), ignored: total("ignored") },
  bootstrap_child_execution: bootstrap.child_summaries[0],
  doctests,
};
if (result.totals.passed !== 202 || result.totals.failed !== 0 || result.totals.ignored !== 3) throw new Error("unexpected top-level totals");
if (doctests.find(({ name }) => name === "feynkit_tensor")?.summary.passed !== 17) throw new Error("unexpected tensor doctest total");
if (doctests.find(({ name }) => name === "spenso_hep_lib")?.summary.passed !== 0) throw new Error("unexpected HEP doctest total");
console.log(JSON.stringify(result, null, 2));
