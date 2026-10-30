import assert from "node:assert/strict";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const evidence = path.resolve(here, "../..");
const output = process.argv[2] ? path.resolve(process.argv[2]) : fs.mkdtempSync(path.join(os.tmpdir(), "qed-oracle-round3-"));
const temporary = !process.argv[2];
fs.mkdirSync(output, { recursive: true });

function run(args) {
  const result = spawnSync("node", args, { encoding: "utf8" });
  if (result.status !== 0) throw new Error(`${args[0]} failed:\n${result.stdout}\n${result.stderr}`);
}

run([path.join(here, "controls.mjs"), path.join(output, "prior")]);
run([
  path.resolve(here, "../review-inputs/heldouts.mjs"),
  evidence,
  path.join(here, "semantic-validate.mjs"),
  path.resolve(here, "../review-inputs/myfnv.mjs"),
  path.join(output, "heldouts"),
]);

const prior = JSON.parse(fs.readFileSync(path.join(output, "prior", "results.json")));
assert.equal(prior.controls, 49);
assert.equal(prior.passed, 49);
assert.ok(prior.results.every((result) => result.pass));

const heldouts = JSON.parse(fs.readFileSync(path.join(output, "heldouts", "results.json")));
const qualified = new Map([
  ["X14-consistent-coef-mirror", "accept"],
  ["X15-Q-label-mirror", "accept"],
  ["X16-overall-factor-full-mirror", "accept"],
  ["X17-routing-mirror", "accept"],
  ["X23-metadata-unknown-owner-field", "reject"],
  ["X24-mirrored-unknown-payload-field", "reject"],
]);
for (const result of heldouts) {
  const expected = qualified.get(result.id) ?? result.expect;
  assert.ok(expected === "accept" || expected === "reject", `${result.id}: unresolved expectation ${expected}`);
  assert.equal(result.actual, expected, result.id);
  assert.equal(result.editError, null, `${result.id}: edit failed`);
}
for (const id of ["X1-extra-numerator-key-diagram-current", "X2-extra-projector-key-cut-current", "X3-extra-numerator-key-member-historical", "X4-extra-prefactor-key-incidence-current"]) {
  const result = heldouts.find((entry) => entry.id === id);
  assert.equal(result.fnv_changed, false);
  assert.ok(result.failures.some((failure) => failure.includes("schema keys unsupported")));
}
for (const id of ["X5-aggregate-x2-full-mirror", "X6-local-drop-i-full-mirror", "X7-local-sign-full-mirror", "X8-local-precedence-trap-neg-pow", "X9-local-precedence-trap-pow-neg", "X13-gc-multiplicity-aggregate"]) {
  const result = heldouts.find((entry) => entry.id === id);
  assert.ok(result.failures.some((failure) => failure.includes("aggregate")));
}
for (const id of ["X10-valid-power-equivalent", "X11-valid-unary-equivalent"]) {
  assert.equal(heldouts.find((entry) => entry.id === id).actual, "accept");
}

const summary = {
  schema: 1,
  prior: { total: prior.controls, passed: prior.passed },
  heldouts: {
    total: heldouts.length,
    accepted: heldouts.filter((result) => result.actual === "accept").length,
    rejected: heldouts.filter((result) => result.actual === "reject").length,
  },
  total: prior.controls + heldouts.length,
  passed: prior.passed + heldouts.length,
};
console.log(JSON.stringify(summary, null, 2));
if (temporary) fs.rmSync(output, { recursive: true, force: true });
