// Reviewer held-out corruptions/valid controls against the frozen round-2 semantic entrypoint.
// Each control copies the four consumed raw folders, edits parsed JSON, and runs ONLY
// proof/semantic-validate.mjs (no manifest, no SHA, no fixed FNV). FNV is recomputed separately.
import fs from "node:fs";
import path from "node:path";
import { spawnSync } from "node:child_process";
const [evidence, semantic, myfnv, out] = process.argv.slice(2);
const E = "xs_qed_eedd_L1", MM = "xs_qed_eemm_L1", PI = "xs_qed_eeee_1pi_L1";
const pv = (j) => j.diagrams[0].payload;
const own = (m, o) => m.diagrams[0].owners.find((x) => x.owner === o);
const deep = (v, f) => typeof v === "string" ? f(v) : Array.isArray(v) ? v.map((x) => deep(x, f)) : v && typeof v === "object" ? (Object.keys(v).forEach((k) => (v[k] = deep(v[k], f))), v) : v;
const must = (s, t) => { if (s === t) throw new Error("edit had no effect"); return t; };
const ALL = ["current", "historical", "metadata-current", "metadata-historical"];
const both = (c, f) => [["current", c, f], ["historical", c, f]];
const live = (c, o, f) => [["metadata-current", c, (j) => { own(j, o).expression = must(own(j, o).expression, f(own(j, o).expression)); }], ["metadata-historical", c, (j) => { own(j, o).expression = must(own(j, o).expression, f(own(j, o).expression)); }]];
const FAC = "feynkit_generator_factor::{}::AutG(1)";
const C = [
  // --- expression-named keys stripped recursively by noExpr (one-sided, should reject) ---
  ["X1-extra-numerator-key-diagram-current", "reject", "one-sided: current diagrams[0].numerator = 'garbage(' (unknown location)", [["current", E, (j) => { j.diagrams[0].numerator = "garbage("; }]]],
  ["X2-extra-projector-key-cut-current", "reject", "one-sided: current payload.cuts[0].projector added", [["current", E, (j) => { pv(j).cuts[0].projector = "2"; }]]],
  ["X3-extra-numerator-key-member-historical", "reject", "one-sided: historical group member gains numerator field", [["historical", MM, (j) => { j.groups[0].members[0].numerator = "0"; }]]],
  ["X4-extra-prefactor-key-incidence-current", "reject", "one-sided: current edges[2][0] incidence object gains numerator_prefactor", [["current", E, (j) => { pv(j).edges[2][0].numerator_prefactor = "-1"; }]]],
  // --- aggregate chain predicate actually firing when live tables are mirrored too ---
  ["X5-aggregate-x2-full-mirror", "reject", "both payload aggregates and both live numerator owners *2 (locals untouched)", [...both(E, (j) => { pv(j).numerator += "*2"; }), ...live(E, "numerator", (s) => s + "*2")]],
  ["X6-local-drop-i-full-mirror", "reject", "edge:4 loses (-i) in both payloads and both live owners (aggregate untouched)", [...both(E, (j) => { const e = pv(j).edges[4][1]; e.numerator = must(e.numerator, e.numerator.replace("(-𝑖)*", "")); }), ...live(E, "edge:4", (s) => s.replace("(-𝑖)*", ""))]],
  ["X7-local-sign-full-mirror", "reject", "edge:4 (-i)->(i) in payloads and live (aggregate untouched): sign only", [...both(E, (j) => { const e = pv(j).edges[4][1]; e.numerator = must(e.numerator, e.numerator.replace("(-𝑖)*", "(𝑖)*")); }), ...live(E, "edge:4", (s) => s.replace("(-𝑖)*", "(𝑖)*"))]],
  ["X8-local-precedence-trap-neg-pow", "reject", "edge:4 (-i) -> -𝑖^3 (= +i) everywhere local", [...both(E, (j) => { const e = pv(j).edges[4][1]; e.numerator = must(e.numerator, e.numerator.replace("(-𝑖)*", "-𝑖^3*")); }), ...live(E, "edge:4", (s) => s.replace("(-𝑖)*", "-𝑖^3*"))]],
  ["X9-local-precedence-trap-pow-neg", "reject", "edge:4 (-i) -> (-𝑖)^3 (= +i) everywhere local", [...both(E, (j) => { const e = pv(j).edges[4][1]; e.numerator = must(e.numerator, e.numerator.replace("(-𝑖)*", "(-𝑖)^3*")); }), ...live(E, "edge:4", (s) => s.replace("(-𝑖)*", "(-𝑖)^3*"))]],
  ["X10-valid-power-equivalent", "accept", "edge:4 (-i) -> 𝑖^3 (= -i) everywhere local: valid equivalent coefficient", [...both(E, (j) => { const e = pv(j).edges[4][1]; e.numerator = must(e.numerator, e.numerator.replace("(-𝑖)*", "𝑖^3*")); }), ...live(E, "edge:4", (s) => s.replace("(-𝑖)*", "𝑖^3*"))]],
  ["X11-valid-unary-equivalent", "accept", "edge:4 (-i)*g -> -𝑖*g everywhere local", [...both(E, (j) => { const e = pv(j).edges[4][1]; e.numerator = must(e.numerator, e.numerator.replace("(-𝑖)*", "-𝑖*")); }), ...live(E, "edge:4", (s) => s.replace("(-𝑖)*", "-𝑖*"))]],
  ["X12-implicit-mult-unsupported", "reject", "edge:4 (-i) -> (-1𝑖) implicit product everywhere local", [...both(E, (j) => { const e = pv(j).edges[4][1]; e.numerator = must(e.numerator, e.numerator.replace("(-𝑖)*", "(-1𝑖)*")); }), ...live(E, "edge:4", (s) => s.replace("(-𝑖)*", "(-1𝑖)*"))]],
  ["X13-gc-multiplicity-aggregate", "reject", "aggregate GC_1^2 -> GC_1^3 in both payloads and live", [...both(E, (j) => { pv(j).numerator = must(pv(j).numerator, pv(j).numerator.replace("GC_1^2", "GC_1^3")); }), ...live(E, "numerator", (s) => s.replace("GC_1^2", "GC_1^3"))]],
  // --- consistent mirrored physics changes (equivalence-preserving; semantic cannot know) ---
  ["X14-consistent-coef-mirror", "accept?", "vertex:0 and aggregate *2 in both payloads and both live owners (consistent)", [...both(E, (j) => { pv(j).vertices[0].numerator = "2*" + pv(j).vertices[0].numerator; pv(j).numerator = "2*" + pv(j).numerator; }), ...live(E, "vertex:0", (s) => "2*" + s), ...live(E, "numerator", (s) => "2*" + s)]],
  ["X15-Q-label-mirror", "reject?", "eedd edge:2 momentum Q(2,..) -> Q(3,..) everywhere (edge owner uses edge 3 momentum)", ALL.map((f) => [f, E, (j) => deep(j, (s) => s.replaceAll("::Q(2,", "::Q(3,"))])],
  ["X16-overall-factor-full-mirror", "reject?", "AutG(1)->AutG(2) in diagram/payload/member factors on both routes; aut field still 1", both(MM, (j) => deep(j, (s) => (s.includes(FAC) ? s.replaceAll(FAC, FAC.replace("(1)", "(2)")) : s)))],
  ["X17-routing-mirror", "reject?", "loop sign of edge 2 flipped in payload LMB on both routes", both(E, (j) => { const s = pv(j).loop_momentum_basis.edge_signatures; s[2].loops[0] = -s[2].loops[0]; })],
  // --- cooked decoding ---
  ["X18-cooked-swap-one-owner", "reject", "historical live vertex:0 swaps two cooked bis symbols only there", [["metadata-historical", E, (j) => { const o = own(j, "vertex:0"); const a = "cooked_e3db38445dbf866e096437d0", b = "cooked_cac23231d644d74d7551011a"; o.expression = must(o.expression, o.expression.replaceAll(a, "TMP").replaceAll(b, a).replaceAll("TMP", b)); }]]],
  ["X19-valid-cooked-global-rename", "accept", "historical live: two cooked symbols swapped consistently in every owner (alpha-equivalent)", [["metadata-historical", E, (j) => { const a = "cooked_e3db38445dbf866e096437d0", b = "cooked_cac23231d644d74d7551011a"; deep(j, (s) => s.replaceAll(a, "TMP").replaceAll(b, a).replaceAll("TMP", b)); }]]],
  ["X20-valid-symmetric-g-swap-live", "accept", "historical live edge:4 g(mink a, mink b) args swapped (symmetric)", [["metadata-historical", E, (j) => { const o = own(j, "edge:4"); const m = /g\((.*),(spenso::\{spenso::representation,spenso::self_dual\}::mink\(4,idenso[^)]*\))\)$/.exec(o.expression); o.expression = must(o.expression, o.expression.replace(m[0], `g(${m[2]},${m[1]})`)); }]]],
  // --- typed axes fully mirrored and self-consistent: only the frozen contract can reject ---
  ["X21-variance-flip-consistent", "reject", "hedge(2) cof variance flipped at vertex:1 and edge:2 in payloads, interfaces, live exprs and ports", ALL.map((f) => [f, E, (j) => { const H = "(3,gammalooprs::{spenso::index}::hedge(2,1))", K = "(3,idenso::{spenso::index}::cooked_0acc2ebe195ecd34e30fcb4c)"; const cof = "spenso::{spenso::representation,spenso::dualizable}::cof"; const low = "spenso::{spenso::lower}::dind("; return deep(j, (s) => { for (const x of [H, K]) { s = s.replaceAll(low + cof + x + ")", "@@").replaceAll(cof + x, low + cof + x + ")").replaceAll("@@", cof + x); } return s; }); }])],
  ["X22-dimension-consistent", "reject", "hedge(2) cof(3)->cof(4) everywhere incl. ports/interfaces/live", ALL.map((f) => [f, E, (j) => deep(j, (s) => s.replaceAll("cof(3,gammalooprs::{spenso::index}::hedge(2,1))", "cof(4,gammalooprs::{spenso::index}::hedge(2,1))").replaceAll("cof(3,idenso::{spenso::index}::cooked_0acc2ebe195ecd34e30fcb4c)", "cof(4,idenso::{spenso::index}::cooked_0acc2ebe195ecd34e30fcb4c)"))])],
  // --- metadata schema ---
  ["X23-metadata-unknown-owner-field", "reject?", "one-sided: current live vertex:0 gains interface:[] record", [["metadata-current", E, (j) => { own(j, "vertex:0").interface = []; }]]],
  ["X24-mirrored-unknown-payload-field", "accept?", "mirrored unknown non-interface payload field", both(E, (j) => { pv(j).bogus = 1; })],
  ["X25-valid-interface-whitespace", "accept", "current vertex:0 interface entry with inserted spaces (same slot)", [["current", E, (j) => { const a = pv(j).vertices[0].num_interface; a[0] = a[0].replace(",gammalooprs", ", gammalooprs"); }]]],
  ["X26-iface-port-hedge-relabel-consistent", "reject", "current vertex:0 interface+live ports+payload use hedge(99) for hedge(6) consistently on current only", [["current", E, (j) => deep(j, (s) => s.replaceAll("hedge(6,1)", "hedge(99,1)"))], ["metadata-current", E, (j) => deep(j, (s) => s.replaceAll("hedge(6,1)", "hedge(99,1)"))]]],
  ["X27-mirrored-hedge-relabel-all", "reject", "hedge(6)<->hedge(7) labels swapped in all four files (vertex/edge incidence violated)", ALL.map((f) => [f, E, (j) => deep(j, (s) => s.replaceAll("hedge(6,1)", "@@").replaceAll("hedge(7,1)", "hedge(6,1)").replaceAll("@@", "hedge(7,1)"))])],
];
const baseF = spawnSync("node", [myfnv, evidence], { encoding: "utf8" }).stdout;
const results = [];
for (const [id, expect, desc, edits] of C) {
  const dir = path.join(out, id);
  fs.mkdirSync(dir, { recursive: true });
  for (const f of ["historical", "current", "metadata-historical", "metadata-current"]) fs.cpSync(path.join(evidence, f), path.join(dir, f), { recursive: true, force: true });
  let editError = null;
  try {
    for (const [folder, name, fn] of edits) {
      const f = path.join(dir, folder, `${name}.json`);
      const before = fs.readFileSync(f, "utf8"); const j = JSON.parse(before); const r = fn(j) ; const v = r && typeof r === "object" && r.diagrams ? r : j;
      const after = JSON.stringify(v, null, 2); if (JSON.stringify(JSON.parse(before)) === JSON.stringify(v)) throw new Error(`no-op edit ${folder}/${name}`);
      fs.writeFileSync(f, after);
    }
  } catch (e) { editError = e.message; }
  const r = spawnSync("node", [semantic, dir], { encoding: "utf8" });
  let failures = []; try { failures = JSON.parse(r.stdout).failures; } catch {}
  const fnvAfter = spawnSync("node", [myfnv, dir], { encoding: "utf8" }).stdout;
  results.push({ id, expect, actual: r.status === 0 ? "accept" : "reject", status: r.status, editError, fnv_changed: fnvAfter !== baseF, desc, failures: failures.slice(0, 3), stderr: r.stderr.slice(0, 200) });
}
fs.writeFileSync(path.join(out, "results.json"), JSON.stringify(results, null, 1) + "\n");
for (const x of results) console.log(`${x.id.padEnd(44)} expect=${x.expect.padEnd(8)} actual=${x.actual} fnvΔ=${x.fnv_changed} ${x.editError ? "EDIT-ERR " + x.editError : ""} ${(x.failures[0] ?? "").slice(0, 120)}`);
