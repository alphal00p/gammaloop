// Mutation harness for the corrected reusable semantic entrypoint. Capture hashes
// and fixed digest literals are deliberately outside this path.
import fs from "node:fs";
import path from "node:path";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const base = path.resolve(here, "../..");
const out = path.resolve(process.argv[2] ?? "target/oracle-numerator-diagnosis/round2-controls");
const semantic = path.join(here, "semantic-validate.mjs");
const fnv = path.resolve(here, "../review-inputs/fnv.mjs");
const baselineFnv = spawnSync("node", [fnv, base], { encoding: "utf8" });
if (baselineFnv.status !== 0) throw new Error(baselineFnv.stderr);
const transposeControls = new Set([
  "valid-iface-consistent-reorder",
  "valid-heldout-eeee-interface-transpose",
]);
fs.mkdirSync(out, { recursive: true });
const E = "xs_qed_eedd_L1", PI = "xs_qed_eeee_1pi_L1", MM = "xs_qed_eemm_L1";
const pv = (j) => j.diagrams[0].payload;
const own = (m, o) => m.diagrams[0].owners.find((x) => x.owner === o);
const QC = "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::Q";
const QH = "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q";
const replaceStrings = (value, from, to) => {
  if (typeof value === "string") return value.replaceAll(from, to);
  if (Array.isArray(value)) return value.map((child) => replaceStrings(child, from, to));
  if (value && typeof value === "object") {
    for (const [key, child] of Object.entries(value)) value[key] = replaceStrings(child, from, to);
  }
  return value;
};

// [id, expectReject, description, edits: [[folder, case, fn]]]
const M = [
  ["iface-delete-entry", true, "current vertex:0 num_interface loses one axis", [["current", E, (j) => pv(j).vertices[0].num_interface.pop()]]],
  ["iface-extra-entry", true, "current vertex:0 num_interface gains a foreign axis", [["current", E, (j) => pv(j).vertices[0].num_interface.push("spenso::{spenso::representation,spenso::self_dual}::mink(4,gammalooprs::{spenso::index}::hedge(99,1))")]]],
  ["iface-wrong-variance", true, "current vertex:0 interface: dind(cof) -> cof", [["current", E, (j) => { const a = pv(j).vertices[0].num_interface; a[4] = a[4].replace(/^spenso::\{spenso::lower\}::dind\((.*)\)$/, "$1"); }]]],
  ["iface-wrong-dimension", true, "current edge:2 interface bis(4)->bis(5)", [["current", E, (j) => { const a = pv(j).edges[2][1].num_interface; a[0] = a[0].replace("bis(4,", "bis(5,"); }]]],
  ["iface-swapped-owner", true, "current vertex:0/vertex:1 interfaces swapped", [["current", E, (j) => { const v = pv(j).vertices; [v[0].num_interface, v[1].num_interface] = [v[1].num_interface, v[0].num_interface]; }]]],
  ["iface-deleted-record", true, "current edge:2 num_interface removed (null) while live order is explicit", [["current", E, (j) => { pv(j).edges[2][1].num_interface = null; }]]],
  ["iface-malformed", true, "current vertex:0 num_interface is a garbage string", [["current", E, (j) => { pv(j).vertices[0].num_interface = "garbage("; }]]],
  ["iface-extra-field", true, "unknown extra_interface field on a group member", [["current", MM, (j) => { j.groups[0].members[0].bogus_interface = ["x"]; }]]],
  ["iface-projector-nonnull", true, "current projector_interface becomes non-null wrong axis list", [["current", E, (j) => { pv(j).projector_interface = ["spenso::{spenso::representation,spenso::self_dual}::bis(4,gammalooprs::{spenso::index}::hedge(0,1))"]; }]]],
  ["valid-iface-consistent-reorder", false, "current vertex:0 interface and its live ports reversed together (consistent transpose)", [["current", E, (j) => { pv(j).vertices[0].num_interface.reverse(); }], ["metadata-current", E, (j) => own(j, "vertex:0").ports.reverse()]]],
  ["meta-variance-mirrored", true, "both live port lists: dind(cof) -> cof on vertex:0 (strings still agree)", [["metadata-current", E, (j) => { const p = own(j, "vertex:0").ports; p[4] = p[4].replace(/^spenso::\{spenso::lower\}::dind\((.*)\)$/, "$1"); }], ["metadata-historical", E, (j) => { const p = own(j, "vertex:0").ports; p[0] = p[0].replace(/^spenso::\{spenso::lower\}::dind\((.*)\)$/, "$1"); }]]],
  ["meta-dimension-mirrored", true, "both live port lists: mink(4)->mink(5) vertex:2", [["metadata-current", E, (j) => { const p = own(j, "vertex:2").ports; p[2] = p[2].replace("mink(4,", "mink(5,"); }], ["metadata-historical", E, (j) => { const p = own(j, "vertex:2").ports; p[2] = p[2].replace("mink(4,", "mink(5,"); }]]],
  ["meta-dimension-onesided", true, "current live ports mink(4)->mink(5) vertex:2", [["metadata-current", E, (j) => { const p = own(j, "vertex:2").ports; p[2] = p[2].replace("mink(4,", "mink(5,"); }]]],
  ["meta-owner-duplicate", true, "current owner table duplicates vertex:0", [["metadata-current", E, (j) => j.diagrams[0].owners.push(structuredClone(own(j, "vertex:0")))]]],
  ["meta-owner-missing-current", true, "current owner table drops edge:2", [["metadata-current", E, (j) => { j.diagrams[0].owners = j.diagrams[0].owners.filter((o) => o.owner !== "edge:2"); }]]],
  ["meta-owner-missing-both", true, "both owner tables drop edge:2", [["metadata-current", E, (j) => { j.diagrams[0].owners = j.diagrams[0].owners.filter((o) => o.owner !== "edge:2"); }], ["metadata-historical", E, (j) => { j.diagrams[0].owners = j.diagrams[0].owners.filter((o) => o.owner !== "edge:2"); }]]],
  ["meta-owner-extra-historical", true, "historical owner table gains vertex:9", [["metadata-historical", E, (j) => j.diagrams[0].owners.push({ ...structuredClone(own(j, "vertex:0")), owner: "vertex:9" })]]],
  ["meta-owner-swapped-mirrored", true, "vertex:0<->vertex:1 owner labels swapped in both live tables", [["metadata-current", E, (j) => { own(j, "vertex:0").owner = "tmp"; own(j, "vertex:1").owner = "vertex:0"; own(j, "tmp").owner = "vertex:1"; }], ["metadata-historical", E, (j) => { own(j, "vertex:0").owner = "tmp"; own(j, "vertex:1").owner = "vertex:0"; own(j, "tmp").owner = "vertex:1"; }]]],
  ["meta-expression-current", true, "current live vertex:0 expression -> 0", [["metadata-current", E, (j) => { own(j, "vertex:0").expression = "0"; }]]],
  ["meta-expression-historical", true, "historical live vertex:0 expression -> 2*expr", [["metadata-historical", E, (j) => { own(j, "vertex:0").expression = "2*" + own(j, "vertex:0").expression; }]]],
  ["meta-flags-mirrored", true, "both live tables mark vertex:0 metric=true", [["metadata-current", E, (j) => { own(j, "vertex:0").metric = true; }], ["metadata-historical", E, (j) => { own(j, "vertex:0").metric = true; }]]],
  ["meta-report", true, "current metadata report retained_count changed", [["metadata-current", E, (j) => { j.report = j.report.replace("retained_count: 1", "retained_count: 2"); }]]],
  ["num-trailing-current", true, "current aggregate numerator *2", [["current", E, (j) => { pv(j).numerator += "*2"; }]]],
  ["num-aggregate-mirrored", true, "both aggregate numerators *2 (local fragments untouched)", [["current", E, (j) => { pv(j).numerator += "*2"; }], ["historical", E, (j) => { pv(j).numerator += "*2"; }]]],
  ["num-local-mirrored", true, "both edge:4 photon fragments lose (-i) factor (aggregate untouched)", [["current", E, (j) => { pv(j).edges[4][1].numerator = pv(j).edges[4][1].numerator.replace("(-𝑖)*", ""); }], ["historical", E, (j) => { pv(j).edges[4][1].numerator = pv(j).edges[4][1].numerator.replace("(-𝑖)*", ""); }]]],
  ["prefactor-mirrored", true, "both numerator_prefactor -> 2", [["current", PI, (j) => { pv(j).numerator_prefactor = "2"; }], ["historical", PI, (j) => { pv(j).numerator_prefactor = "2"; }]]],
  ["prefactor-current", true, "current numerator_prefactor -> 2", [["current", PI, (j) => { pv(j).numerator_prefactor = "2"; }]]],
  ["projector-variance-mirrored", true, "both projectors: u bis(hedge0) -> mink (uncontracted)", [["current", E, (j) => { pv(j).projector = pv(j).projector.replace("self_dual}::bis(4,gammalooprs::{spenso::index}::hedge(0,1))", "self_dual}::mink(4,gammalooprs::{spenso::index}::hedge(0,1))"); }], ["historical", E, (j) => { pv(j).projector = pv(j).projector.replace("self_dual}::bis(4,gammalooprs::{spenso::index}::hedge(0,1))", "self_dual}::mink(4,gammalooprs::{spenso::index}::hedge(0,1))"); }]]],
  ["overall-factor-mirrored", true, "both member overall_factor sign flipped (diagram factor untouched)", [["current", MM, (j) => { j.groups[0].members[0].overall_factor = j.groups[0].members[0].overall_factor.replace("InternalFermionLoopSign(-1)", "InternalFermionLoopSign(1)"); }], ["historical", MM, (j) => { j.groups[0].members[0].overall_factor = j.groups[0].members[0].overall_factor.replace("InternalFermionLoopSign(-1)", "InternalFermionLoopSign(1)"); }]]],
  ["ratio-current", true, "current member ratio -> 2", [["current", MM, (j) => { j.groups[0].members[0].ratio = "2"; }]]],
  ["ratio-mirrored", true, "both member ratios -> 2", [["current", MM, (j) => { j.groups[0].members[0].ratio = "2"; }], ["historical", MM, (j) => { j.groups[0].members[0].ratio = "2"; }]]],
  ["master-mirrored", true, "both group masters -> 1 (nonexistent member)", [["current", MM, (j) => { j.groups[0].master = 1; }], ["historical", MM, (j) => { j.groups[0].master = 1; }]]],
  ["member-source-mirrored", true, "both member source_index -> 1", [["current", MM, (j) => { j.groups[0].members[0].source_index = 1; }], ["historical", MM, (j) => { j.groups[0].members[0].source_index = 1; }]]],
  ["routing-current", true, "current loop routing sign flipped", [["current", E, (j) => { const s = pv(j).loop_momentum_basis.edge_signatures; s[2].loops[0] = -s[2].loops[0]; }]]],
  ["ctor-foreign-in-current", true, "one current Q uses the historical (tagged) spelling", [["current", E, (j) => { pv(j).numerator = pv(j).numerator.replace(QC, QH); }]]],
  ["ctor-namespace", true, "current Q namespace gammalooprs -> other", [["current", E, (j) => { pv(j).numerator = pv(j).numerator.replace(QC, QC.replace(/^gammalooprs/, "other")); }]]],
  ["ctor-label", true, "current Q -> K constructor (label mismatch)", [["current", E, (j) => { pv(j).numerator = pv(j).numerator.replace(QC + "(", QC.replace(/::Q$/, "::K") + "("); }]]],
  ["ctor-prefix-extension", true, "current edge index head renamed edgeX (prefix of mapped spelling)", [["current", E, (j) => { pv(j).numerator = pv(j).numerator.replaceAll("gammalooprs::{spenso::index}::edge(", "gammalooprs::{spenso::index}::edgeX("); }], ["historical", E, (j) => { pv(j).numerator = pv(j).numerator.replaceAll("gammalooprs::{spenso::index,spenso::index-label:e}::edge(", "gammalooprs::{spenso::index,spenso::index-label:e}::edgeX("); }]]],
  ["ratio-mirrored-3", true, "both member ratios -> 3", [["current", MM, (j) => { j.groups[0].members[0].ratio = "3"; }], ["historical", MM, (j) => { j.groups[0].members[0].ratio = "3"; }]]],
  ["num-trailing-mirrored-3", true, "both aggregate numerators *3", [["current", MM, (j) => { pv(j).numerator += "*3"; }], ["historical", MM, (j) => { pv(j).numerator += "*3"; }]]],
  // Correction held-outs, independent of the reviewer's original matrix.
  ["heldout-foreign-symbol-mirrored", true, "foreign coupling head mirrored through both payloads and live expressions", [["current", MM, (j) => replaceStrings(j, "UFO::{}::GC_3", "Foreign::{}::GC_3")], ["historical", MM, (j) => replaceStrings(j, "UFO::{}::GC_3", "Foreign::{}::GC_3")], ["metadata-current", MM, (j) => replaceStrings(j, "UFO::{}::GC_3", "Foreign::{}::GC_3")], ["metadata-historical", MM, (j) => replaceStrings(j, "UFO::{}::GC_3", "Foreign::{}::GC_3")]]],
  ["heldout-inferred-order-reorder", true, "missing-interface owner live order changed", [["metadata-current", MM, (j) => own(j, "vertex:0").ports.reverse()]]],
  ["heldout-model-mirrored", true, "captured model parameter changed on both routes", [["metadata-current", MM, (j) => { j.model.parameters[0].name += "_mutated"; }], ["metadata-historical", MM, (j) => { j.model.parameters[0].name += "_mutated"; }]]],
  ["heldout-unknown-historical-interface", true, "unknown historical group interface field", [["historical", MM, (j) => { j.groups[0].members[0].foreign_interface = []; }]]],
  ["heldout-diagram-id-mirrored", true, "diagram id changed consistently across payload and live tables", [["current", MM, (j) => { j.diagrams[0].id = "00000000000000000000000000000000"; j.groups[0].members[0].id = j.diagrams[0].id; }], ["historical", MM, (j) => { j.diagrams[0].id = "00000000000000000000000000000000"; j.groups[0].members[0].id = j.diagrams[0].id; }], ["metadata-current", MM, (j) => { j.diagrams[0].id = "00000000000000000000000000000000"; }], ["metadata-historical", MM, (j) => { j.diagrams[0].id = "00000000000000000000000000000000"; }]]],
  // valid controls (must be accepted)
  ["valid-identity", false, "rewrite with no change", [["current", E, () => {}]]],
  ["valid-key-order", false, "current JSON object keys reordered", [["current", E, (j) => { const p = pv(j); const keys = Object.keys(p).reverse(); const q = Object.fromEntries(keys.map((k) => [k, p[k]])); j.diagrams[0].payload = q; }]]],
  ["valid-owner-order", false, "live owner tables reordered", [["metadata-current", E, (j) => j.diagrams[0].owners.reverse()], ["metadata-historical", PI, (j) => j.diagrams[0].owners.reverse()]]],
  ["valid-heldout-eeee-interface-transpose", false, "second case explicit interface and live ports transposed consistently", [["current", PI, (j) => pv(j).vertices[1].num_interface.reverse()], ["metadata-current", PI, (j) => own(j, "vertex:1").ports.reverse()]]],
];

function copyEvidence(dir) {
  fs.mkdirSync(dir, { recursive: true });
  for (const folder of ["historical", "current", "metadata-historical", "metadata-current"]) {
    fs.cpSync(path.join(base, folder), path.join(dir, folder), { recursive: true });
  }
}

function runSemantic(dir) {
  const result = spawnSync("node", [semantic, dir], { encoding: "utf8" });
  const capture = spawnSync("node", [fnv, dir], { encoding: "utf8" });
  let failures = [];
  if (result.stdout.trim()) failures = JSON.parse(result.stdout).failures.slice(0, 3);
  return {
    status: result.status,
    failures,
    stderr: result.stderr.trim(),
    capture_fnv_changed: capture.stdout !== baselineFnv.stdout,
  };
}
const results = [];
for (const [id, expectReject, desc, edits] of M) {
  const dir = path.join(out, id);
  fs.rmSync(dir, { recursive: true, force: true });
  copyEvidence(dir);
  for (const [folder, name, fn] of edits) {
    const f = path.join(dir, folder, `${name}.json`);
    const j = JSON.parse(fs.readFileSync(f, "utf8")); fn(j);
    fs.writeFileSync(f, JSON.stringify(j, null, 2));
  }
  const result = runSemantic(dir);
  const verdict = (s) => (s === 0 ? "accept" : "reject");
  const expected = expectReject ? "reject" : "accept", actual = verdict(result.status);
  const captureExpectationMet = !transposeControls.has(id) || result.capture_fnv_changed;
  results.push({ id, expected, actual, pass: expected === actual && captureExpectationMet, capture_fnv_changed: result.capture_fnv_changed, desc, failures: result.failures, stderr: result.stderr });
}
const output = { schema: 1, controls: results.length, passed: results.filter((result) => result.pass).length, results };
fs.mkdirSync(out, { recursive: true });
fs.writeFileSync(path.join(out, "results.json"), `${JSON.stringify(output, null, 2)}\n`);
console.log(JSON.stringify(output, null, 2));
process.exit(results.every((result) => result.pass) ? 0 : 1);
