import assert from "node:assert/strict";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const evidence = path.resolve(here, "../..");
const currentChecker = path.join(here, "semantic-validate.mjs");
const oldChecker = path.resolve(here, "../../round3/proof/semantic-validate.mjs");
const fnv = path.resolve(here, "../review-inputs/myfnv.mjs");
const output = process.argv[2] ? path.resolve(process.argv[2]) : fs.mkdtempSync(path.join(os.tmpdir(), "qed-edge-tuples-"));
const temporary = !process.argv[2];
const folders = ["historical", "current", "metadata-historical", "metadata-current"];
const cases = [
  ["xs_qed_eedd_L1", 0, "0"],
  ["xs_qed_eemm_L1", 4, "-1"],
  ["xs_qed_eeee_1pi_L1", 1, "garbage("],
];

function copyEvidence(destination) {
  fs.mkdirSync(destination, { recursive: true });
  for (const folder of folders) fs.cpSync(path.join(evidence, folder), path.join(destination, folder), { recursive: true });
}

function mutate(directory, folder, name, edit) {
  const file = path.join(directory, folder, `${name}.json`);
  const value = JSON.parse(fs.readFileSync(file, "utf8"));
  edit(value.diagrams[0].payload.edges);
  fs.writeFileSync(file, `${JSON.stringify(value, null, 2)}\n`);
}

function run(checker, directory) {
  const result = spawnSync("node", [checker, directory], { encoding: "utf8" });
  let failures = [];
  try { failures = JSON.parse(result.stdout).failures; } catch {}
  return { status: result.status, actual: result.status === 0 ? "accept" : "reject", failures: failures.slice(0, 4) };
}

const controls = [];
for (const [name, edgeIndex, rogueNumerator] of cases) {
  const variants = [
    ["append-data", (edges) => { const extra = structuredClone(edges[edgeIndex][1]); delete extra.num_interface; extra.numerator = rogueNumerator; edges[edgeIndex].push(extra); }],
    ["append-incidence", (edges) => edges[edgeIndex].push(structuredClone(edges[edgeIndex][0]))],
    ["truncate", (edges) => { edges[edgeIndex] = [edges[edgeIndex][0]]; }],
    ["swap", (edges) => { edges[edgeIndex] = [edges[edgeIndex][1], edges[edgeIndex][0]]; }],
    ["wrong-position", (edges) => { edges[edgeIndex] = [structuredClone(edges[edgeIndex][1]), edges[edgeIndex][1]]; }],
    ["wrong-type", (edges) => { edges[edgeIndex] = { incidence: edges[edgeIndex][0], data: edges[edgeIndex][1] }; }],
  ];
  for (const [variant, edit] of variants) {
    for (const mirrored of [false, true]) {
      controls.push({
        id: `${name}-${mirrored ? "mirrored" : "onesided"}-${variant}`,
        name,
        expected: "reject",
        edit,
        folders: mirrored ? ["historical", "current"] : ["current"],
        original_blocker: mirrored && (variant === "append-data" || variant === "append-incidence"),
      });
    }
  }
  controls.push({
    id: `${name}-valid-component-key-reorder`,
    name,
    expected: "accept",
    folders: ["current"],
    edit: (edges) => {
      edges[edgeIndex] = edges[edgeIndex].map((component) =>
        Object.fromEntries(Object.entries(component).reverse()),
      );
    },
    original_blocker: false,
  });
}

fs.mkdirSync(output, { recursive: true });
const results = [];
for (const control of controls) {
  const directory = path.join(output, control.id);
  copyEvidence(directory);
  for (const folder of control.folders) mutate(directory, folder, control.name, control.edit);
  const old = run(oldChecker, directory), corrected = run(currentChecker, directory);
  const capture = spawnSync("node", [fnv, directory], { encoding: "utf8" });
  const result = {
    id: control.id,
    expected: control.expected,
    old: old.actual,
    corrected: corrected.actual,
    pass: corrected.actual === control.expected,
    original_blocker: control.original_blocker,
    fnv_status: capture.status,
    corrected_failures: corrected.failures,
  };
  if (control.original_blocker) assert.equal(old.actual, "accept", `${control.id}: frozen round3 must reproduce acceptance`);
  results.push(result);
}
assert.ok(results.every((result) => result.pass));
const summary = {
  schema: 1,
  controls: results.length,
  negative: results.filter((result) => result.expected === "reject").length,
  positive: results.filter((result) => result.expected === "accept").length,
  passed: results.filter((result) => result.pass).length,
  old_blocker_acceptances: results.filter((result) => result.original_blocker && result.old === "accept").length,
  results,
};
fs.writeFileSync(path.join(output, "results.json"), `${JSON.stringify(summary, null, 2)}\n`);
console.log(JSON.stringify({
  schema: summary.schema,
  controls: summary.controls,
  negative: summary.negative,
  positive: summary.positive,
  passed: summary.passed,
  old_blocker_acceptances: summary.old_blocker_acceptances,
}, null, 2));
if (temporary) fs.rmSync(output, { recursive: true, force: true });
