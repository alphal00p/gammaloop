import assert from "node:assert/strict";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const evidence = path.resolve(here, "../..");
const output = process.argv[2] ? path.resolve(process.argv[2]) : fs.mkdtempSync(path.join(os.tmpdir(), "qed-round4-controls-"));
const temporary = !process.argv[2];
fs.mkdirSync(output, { recursive: true });

function run(script, args) {
  const result = spawnSync("node", [script, ...args], { encoding: "utf8" });
  if (result.status !== 0) throw new Error(`${script} failed:\n${result.stdout}\n${result.stderr}`);
  return result.stdout.trim() ? JSON.parse(result.stdout) : {};
}

const prior = run(path.join(here, "prior-controls-all.mjs"), [path.join(output, "prior")]);
const tuples = run(path.join(here, "tuple-controls.mjs"), [path.join(output, "tuples")]);
const semantic = path.join(here, "semantic-validate.mjs");
for (const [script, folder] of [["probes.mjs", "probes"], ["probes2.mjs", "probes2"]]) {
  const result = spawnSync("node", [path.resolve(here, `../review-inputs/${script}`), evidence, semantic, path.join(output, folder)], { encoding: "utf8" });
  if (result.status !== 0) throw new Error(`${script} failed:\n${result.stdout}\n${result.stderr}`);
}
const probes = JSON.parse(fs.readFileSync(path.join(output, "probes", "results.json")));
const probes2 = JSON.parse(fs.readFileSync(path.join(output, "probes2", "results.json")));
const qualifiedProbeIds = new Set(["G18-key-edges-extra-numerator-array", "V1-all-files-key-reversal"]);
const probeDeviations = [];
for (const result of [...probes, ...probes2]) {
  const expected = result.expect === "capture" || result.expect === "accept?" ? "accept" : result.expect;
  if (qualifiedProbeIds.has(result.id)) probeDeviations.push({ id: result.id, expected, actual: result.semantic });
  else assert.equal(result.semantic, expected, result.id);
  assert.equal(result.editError, null, result.id);
}
assert.deepEqual(probeDeviations, [
  { id: "G18-key-edges-extra-numerator-array", expected: "reject", actual: "accept" },
  { id: "V1-all-files-key-reversal", expected: "accept", actual: "reject" },
]);
for (const id of [
  "V3-current-iface-transpose-eedd-edge2",
  "V4-coef-equivalent-i-squared",
  "V5-coef-double-neg",
  "V6-iface-whitespace-all",
  "V7-owner-table-reverse-both",
]) assert.equal(probes.find((result) => result.id === id).semantic, "accept", id);
for (const id of ["T2b-gc-split-power", "T7-neg-pow-precedence-aggregate"]) {
  assert.equal(probes2.find((result) => result.id === id).semantic, "accept", id);
}

const summary = {
  schema: 1,
  prior,
  tuples,
  primary: { total: prior.total + tuples.controls, passed: prior.passed + tuples.passed },
  reviewer_probes: { total: probes.length + probes2.length, matched: probes.length + probes2.length - probeDeviations.length, qualified: probeDeviations },
  prior_expectation_ownership: { reviewer_fixed: 70, author_finalized: 6 },
};
console.log(JSON.stringify(summary, null, 2));
if (temporary) fs.rmSync(output, { recursive: true, force: true });
