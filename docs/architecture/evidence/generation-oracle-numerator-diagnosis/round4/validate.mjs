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
requireHash("docs/architecture/evidence/generation-oracle-numerator-diagnosis/round3/receipt.json", "c0b5f1a041014420375e7d0789230ecfd79469774ebfe8b4e485640e46d41dc8");
requireHash("docs/architecture/evidence/generation-oracle-numerator-diagnosis/round3/manifest.sha256", "0d6b6a34b96f508f2e7ae9f83a3e89990031e38dc23949ceceb4a4603b533674");
requireHash("docs/architecture/evidence/generation-oracle-numerator-diagnosis/round3/changed-inventory.json", "a6f97a041339fd3414406540549743f172f0f9d8ec66e8c2c9b12af9f0757c33");
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
const contract = JSON.parse(fs.readFileSync(path.join(here, "proof", "semantic-contract.json")));
assert.deepEqual(semantic.failures, []);
assert.equal(semantic.report.reduce((sum, row) => sum + row.axes, 0), 92);
assert.equal(semantic.report.reduce((sum, row) => sum + row.perms.length, 0), 6);
assert.equal(contract.edge_tuple.length, 2);
assert.deepEqual(contract.edge_tuple.incidence_signatures, [["source", "target"]]);
assert.equal(contract.edge_tuple.data_signatures.length, 2);
assert.equal(controls.primary.total, 115);
assert.equal(controls.primary.passed, 115);
assert.equal(controls.reviewer_probes.total, 43);
assert.equal(controls.reviewer_probes.matched, 41);

console.log(JSON.stringify({
  schema: 1,
  capture_rows: capture.rows.length,
  semantic_rows: semantic.report,
  object_key_paths: {
    payload: Object.keys(contract.schema_signatures.payload).length,
    metadata: Object.keys(contract.schema_signatures.metadata).length,
  },
  edge_tuple: contract.edge_tuple,
  controls,
  documentation,
  retained_adapted_elfs: 2,
}, null, 2));
