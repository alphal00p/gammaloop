import assert from "node:assert/strict";
import crypto from "node:crypto";
import fs from "node:fs";
import path from "node:path";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const evidence = path.resolve(here, "..");
const workspace = path.resolve(here, "../../../../..");
const sha = (bytes) => crypto.createHash("sha256").update(bytes).digest("hex");
for (const line of fs.readFileSync(path.join(here, "manifest.sha256"), "utf8").trim().split("\n")) {
  const [expected, relative] = line.split(/\s+/, 2);
  assert.equal(sha(fs.readFileSync(path.join(here, relative))), expected, relative);
}
const inventory = JSON.parse(fs.readFileSync(path.join(here, "changed-inventory.json"), "utf8"));
for (const entry of inventory.files) {
  const bytes = fs.readFileSync(path.join(workspace, entry.path));
  assert.equal(bytes.length, entry.bytes, entry.path);
  assert.equal(sha(bytes), entry.sha256, entry.path);
}

function requireHash(relative, expected) {
  assert.equal(sha(fs.readFileSync(path.join(workspace, relative))), expected, relative);
}
requireHash("docs/architecture/evidence/generation-oracle-numerator-diagnosis/round2/receipt.json", "9c11b9f748e79509d91f122107a109cb8bc7a826150ee70a8480f1140b07af08");
requireHash("docs/architecture/evidence/generation-oracle-numerator-diagnosis/round2/manifest.sha256", "e24c60731a26698c5904f5e736cc2fedfcf03889c51401a2d56348d13060b2e8");
requireHash("docs/architecture/evidence/generation-oracle-numerator-diagnosis/round2/changed-inventory.json", "6b8ef84ba32f1758ee0a54ca49c6f163df77a70a111608a6f95c40b66fee1477");
requireHash("target/oracle-numerator-diagnosis/adapted-historical/review-rederived-generation_oracle", "31e9ea393246f513f6c7fdbf55fe8a425213ae3276751c95cf8ddbeda47b01ef");
requireHash("target/oracle-numerator-diagnosis/adapted-current/review-rederived-generation_oracle", "1a76aa5134431c7436340c6f231b48e5fba66a67e47e19eaaf376e386b466c50");

function run(script, args = []) {
  const result = spawnSync("node", [path.join(here, "proof", script), ...args], { encoding: "utf8" });
  if (result.status !== 0) throw new Error(`${script} failed:\n${result.stdout}\n${result.stderr}`);
  return JSON.parse(result.stdout);
}
const capture = run("capture-validate.mjs");
const semantic = run("semantic-validate.mjs", [evidence]);
const controls = run("controls-all.mjs");
const documentation = run("canonical-doc-check.mjs");
assert.deepEqual(semantic.failures, []);
assert.equal(semantic.report.reduce((sum, row) => sum + row.axes, 0), 92);
assert.equal(semantic.report.reduce((sum, row) => sum + row.perms.length, 0), 6);
assert.equal(controls.total, 76);
assert.equal(controls.passed, 76);

console.log(JSON.stringify({
  schema: 1,
  capture_rows: capture.rows.length,
  semantic_rows: semantic.report,
  schema_paths: {
    payload: Object.keys(JSON.parse(fs.readFileSync(path.join(here, "proof", "semantic-contract.json"))).schema_signatures.payload).length,
    metadata: Object.keys(JSON.parse(fs.readFileSync(path.join(here, "proof", "semantic-contract.json"))).schema_signatures.metadata).length,
  },
  controls,
  documentation,
  retained_adapted_elfs: 2,
}, null, 2));
