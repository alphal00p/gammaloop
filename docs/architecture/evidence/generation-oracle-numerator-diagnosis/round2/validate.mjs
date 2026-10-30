import assert from "node:assert/strict";
import crypto from "node:crypto";
import fs from "node:fs";
import path from "node:path";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const evidence = path.resolve(here, "..");
for (const line of fs.readFileSync(path.join(here, "manifest.sha256"), "utf8").trim().split("\n")) {
  const [expected, relative] = line.split(/\s+/, 2);
  const bytes = fs.readFileSync(path.join(here, relative));
  assert.equal(crypto.createHash("sha256").update(bytes).digest("hex"), expected, relative);
}
const inventory = JSON.parse(fs.readFileSync(path.join(here, "changed-inventory.json"), "utf8"));
const workspace = path.resolve(here, "../../../../..");
for (const entry of inventory.files) {
  const bytes = fs.readFileSync(path.join(workspace, entry.path));
  assert.equal(bytes.length, entry.bytes, entry.path);
  assert.equal(crypto.createHash("sha256").update(bytes).digest("hex"), entry.sha256, entry.path);
}

function run(script, args = []) {
  const result = spawnSync("node", [path.join(here, "proof", script), ...args], {
    encoding: "utf8",
  });
  if (result.status !== 0) throw new Error(`${script} failed:\n${result.stdout}\n${result.stderr}`);
  return JSON.parse(result.stdout);
}

const capture = run("capture-validate.mjs");
const semantic = run("semantic-validate.mjs", [evidence]);
const documentation = run("canonical-doc-check.mjs");
const controls = JSON.parse(fs.readFileSync(path.join(here, "proof", "control-results.json")));
assert.equal(controls.controls, 49);
assert.equal(controls.passed, 49);
assert.ok(controls.results.every((result) => result.pass));
assert.deepEqual(semantic.failures, []);
assert.equal(semantic.report.reduce((sum, row) => sum + row.axes, 0), 92);
assert.equal(semantic.report.reduce((sum, row) => sum + row.perms.length, 0), 6);
assert.ok(semantic.report.every((row) => row.cooked === 14));

console.log(JSON.stringify({
  schema: 1,
  capture_rows: capture.rows.length,
  semantic_rows: semantic.report,
  mutation_controls: { total: controls.controls, passed: controls.passed },
  documentation,
}, null, 2));
